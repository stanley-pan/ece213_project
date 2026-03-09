#include "fastsp.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>

#define CUDA_CHECK(call)                                         \
    do {                                                         \
        cudaError_t _e = (call);                                 \
        if (_e != cudaSuccess)                                   \
            throw std::runtime_error(std::string("CUDA error: ") \
                                     + cudaGetErrorString(_e));  \
    } while (0)


// helpers
namespace {

    inline long long c2(long long n) { return n >= 2 ? n * (n - 1) / 2 : 0; }

    inline char upper(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
    inline bool is_lower(char c) { return std::islower(static_cast<unsigned char>(c)) != 0; }
    inline bool is_alpha(char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0; }

    class Alignment {
        public:
            std::vector<std::string>            names;
            std::vector<std::string>            seqs;
            std::unordered_map<std::string,int> name_to_idx;
            int cols = 0;
            int n    = 0;
    };

    Alignment read_alignment(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot open: " + path);

        Alignment aln;
        std::string line, cur_name, cur_seq;

        auto flush = [&]() {
            if (cur_name.empty()) return;
            aln.name_to_idx[cur_name] = (int)aln.names.size();
            aln.names.push_back(cur_name);
            aln.seqs .push_back(cur_seq);
            cur_name.clear(); cur_seq.clear();
        };

        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                flush();
                std::string nm = line.substr(1);
                auto ws = nm.find_first_of(" \t\r");
                if (ws != std::string::npos) nm.resize(ws);
                cur_name = nm;
            } else {
                for (char c : line)
                    if (c != ' ' && c != '\t' && c != '\r' && c != '\n')
                        cur_seq.push_back(c);
            }
        }
        flush();

        if (aln.names.empty()) throw std::runtime_error("No sequences found in: " + path);
        aln.cols = (int)aln.seqs[0].size();
        aln.n    = (int)aln.names.size();
        for (int i = 0; i < aln.n; i++)
            if ((int)aln.seqs[i].size() != aln.cols)
                throw std::runtime_error("Inconsistent row lengths in: " + path);
        return aln;
    }

    std::unordered_set<char> dash_chars(const std::string& extras) {
        std::unordered_set<char> g = {'-', '?'};
        for (char c : extras) g.insert(c);
        return g;
    }

} // namespace

/**
 * kernel 1: compute prefix sum over reference columns
 * 
 * d_refColInd[i * refColCount + c] = number of non-gap, non-masked chars
 * in sequence i strictly before column c.
 */
__global__ void computeRefColInd(
    const char* d_refSeqs,
    int*        d_refColInd,
    int         n,
    int         refColCount,
    bool        maskLowerRef
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int count = 0;
    for (int c = 0; c < refColCount; c++) {
        d_refColInd[c * n + i] = count; 

        char rc = d_refSeqs[c * n + i];
        char RC = maskLowerRef ? rc : (rc >= 'a' && rc <= 'z') ? (char)(rc - 32) : rc;

        bool isGap       = (RC == '-' || RC == '?');
        bool isAlpha     = (RC >= 'A' && RC <= 'Z') || (RC >= 'a' && RC <= 'z');
        bool isInsertion = maskLowerRef && (rc >= 'a' && rc <= 'z');

        if (!isGap && isAlpha && !isInsertion) count++;
    }
}

__global__ void scoreColumns(
    const int*  d_s,
    const int*  d_refColInd,
    const char* d_refSeqs,
    const int*  d_charsPerEstimatedColumn,
    int         n,
    int         refColCount,
    int         estColCount,
    int         k1,
    bool        maskLowerRef,

    
    long long*  d_sharedHom,
    long long*  d_totalHom,
    int*        d_correctCols,
    int*        d_effectiveRefCols,
    int*        d_singletonRefCols,
    long long*  d_insertRefTotal,
    int*        d_histogram
) {
    int c = blockIdx.x;
    if (c >= refColCount) return;

    int* hist = d_histogram + c * estColCount;

    //  pass 1: scatter into histogram 
    int localRef    = 0;
    int localInsert = 0;

    for (int i = threadIdx.x; i < n; i += blockDim.x) {
        char rc = d_refSeqs[c * n + i];
        char RC = maskLowerRef ? rc : ((rc >= 'a' && rc <= 'z') ? (char)(rc - 32) : rc);

        bool isGap       = (RC == '-' || RC == '?');
        bool isAlpha     = (RC >= 'A' && RC <= 'Z') || (RC >= 'a' && RC <= 'z');
        bool isInsertion = maskLowerRef && (rc >= 'a' && rc <= 'z');

        if (isGap || !isAlpha) continue;
        if (isInsertion) { localInsert++; continue; }

        localRef++;

        int k = d_refColInd[c * n + i];
        int y = d_s[i * k1 + k];
        if (y >= 0) atomicAdd(&hist[y], 1);
    }

    // block-level reduce localRef and localInsert
    __shared__ int sh_refCount;
    __shared__ int sh_insertCount;
    if (threadIdx.x == 0) { sh_refCount = 0; sh_insertCount = 0; }
    __syncthreads();

    atomicAdd(&sh_refCount, localRef);
    atomicAdd(&sh_insertCount, localInsert);
    __syncthreads();

    // pass 2: reduce histogram to homology count 
    // each thread handles a stripe of est columns
    __shared__ unsigned long long sh_sharedHom;
    __shared__ int       sh_correctCol;
    if (threadIdx.x == 0) { sh_sharedHom = 0; sh_correctCol = 0; }
    __syncthreads();

    unsigned long long localHom = 0ULL;
    int       localCorrect = 0;

    for (int j = threadIdx.x; j < estColCount; j += blockDim.x) {
        int cnt = hist[j];
        if (cnt == 0) continue;

        localHom += (unsigned long long)cnt * (cnt - 1) / 2;

        if (cnt == sh_refCount && sh_refCount >= 2 && d_charsPerEstimatedColumn[j] == cnt)
            localCorrect = 1;

        hist[j] = 0;
    }

    atomicAdd(&sh_sharedHom,  localHom);
    atomicAdd(&sh_correctCol, localCorrect);
    __syncthreads();

    if (threadIdx.x == 0) {
        int   rc        = sh_refCount;
        long long total = (long long)rc * (rc - 1) / 2;

        d_sharedHom[c]        = sh_sharedHom;
        d_totalHom[c]         = total;
        d_effectiveRefCols[c] = (rc >= 2) ? 1 : 0;
        d_singletonRefCols[c] = (rc == 1) ? 1 : 0;
        d_insertRefTotal[c]   = sh_insertCount;
        d_correctCols[c]      = sh_correctCol;

        
    }
}

// main function
FastSPResult run_fastsp_cuda(const FastSPOptions& opt) {
    const auto dashChars = dash_chars(opt.gap_chars);
    auto       is_gap    = [&](char c) { return dashChars.count(c) > 0; };

    std::cerr << "Reference alignment: " << opt.reference_path << " ...\n";
    std::cerr << "Estimated alignment: " << opt.estimated_path << " ...\n";

    Alignment ref = read_alignment(opt.reference_path);
    Alignment est = read_alignment(opt.estimated_path);

    bool swapped = false;
    if (est.cols < ref.cols && !opt.mask_lower_est && !opt.mask_lower_ref) {
        swapped = true;
        std::swap(ref, est);
        std::cerr << "(alignments swapped: estimated was shorter)\n";
    }

    const int n            = ref.n;
    const int refColCount  = ref.cols;
    const int estColCount  = est.cols;

    if (est.n != n)
        throw std::runtime_error("Ref and est have different sequence counts.");

    std::vector<std::string> est_ordered(n);
    for (int i = 0; i < n; i++) {
        auto it = est.name_to_idx.find(ref.names[i]);
        if (it == est.name_to_idx.end())
            throw std::runtime_error("Sequence missing from estimated: " + ref.names[i]);
        est_ordered[i] = est.seqs[it->second];
    }

    // k1
    int k1 = 0;
    for (int i = 0; i < n; i++) {
        int nucInd = 0;
        for (char c : est_ordered[i]) {
            char C = opt.mask_lower_est ? c : upper(c);
            if (!is_gap(C) && is_alpha(C)) nucInd++;
        }
        k1 = std::max(k1, nucInd);
    }
    if (k1 <= 0) throw std::runtime_error("Estimated alignment has no letters.");

    // Build s[i][k] and charsPerEstimatedColumn[j]
    std::vector<int> h_s(n * k1, 0); // flat host copy
    std::vector<int> charsPerEstimatedColumn(estColCount, 0);
    long long insertEstTotal = 0;

    for (int i = 0; i < n; i++) {
        int nucInd = 0;
        for (int j = 0; j < estColCount; j++) {
            char c = est_ordered[i][j];
            char C = opt.mask_lower_est ? c : upper(c);
            if (is_gap(C) || !is_alpha(C)) continue;
            if (opt.mask_lower_est && is_lower(c)) {
                h_s[i * k1 + nucInd] = -1;
                insertEstTotal++;
            } else {
                h_s[i * k1 + nucInd] = j;
                charsPerEstimatedColumn[j]++;
            }
            nucInd++;
        }
    }

    // Estimated homology stats
    long long totalHomologiesInEsitmated = 0;
    long long effectiveEstColumns = 0;
    long long singletonEstColumns = 0;
    for (int j = 0; j < estColCount; j++) {
        totalHomologiesInEsitmated += c2(charsPerEstimatedColumn[j]);
        if      (charsPerEstimatedColumn[j] == 1) singletonEstColumns++;
        else if (charsPerEstimatedColumn[j] >= 2) effectiveEstColumns++;
    }

    // Flatten reference sequences: h_refSeqs[i * refColCount + c]
    std::vector<char> h_refSeqs(n * refColCount);
    for (int i = 0; i < n; i++)
        for (int c = 0; c < refColCount; c++)
            h_refSeqs[c * n + i] = ref.seqs[i][c];

    /**
     * MALLOC + MEMCPY to device
     */
    int*  d_s                       = nullptr;
    int*  d_charsPerEstimatedColumn = nullptr;
    char* d_refSeqs                 = nullptr;
    int*  d_refColInd               = nullptr;

    long long* d_sharedHom        = nullptr;
    long long* d_totalHom         = nullptr;
    int*       d_correctCols      = nullptr;
    int*       d_effectiveRefCols = nullptr;
    int*       d_singletonRefCols = nullptr;
    long long* d_insertRefTotal   = nullptr;

    cudaEvent_t gpuStart, gpuStop;
    CUDA_CHECK(cudaEventCreate(&gpuStart));
    CUDA_CHECK(cudaEventCreate(&gpuStop));
    CUDA_CHECK(cudaEventRecord(gpuStart));

    CUDA_CHECK(cudaMalloc(&d_s,                       (size_t)n * k1          * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_charsPerEstimatedColumn, (size_t)estColCount     * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_refSeqs,                 (size_t)n * refColCount * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&d_refColInd,               (size_t)n * refColCount * sizeof(int)));

    CUDA_CHECK(cudaMalloc(&d_sharedHom,        (size_t)refColCount * sizeof(long long)));
    CUDA_CHECK(cudaMalloc(&d_totalHom,         (size_t)refColCount * sizeof(long long)));
    CUDA_CHECK(cudaMalloc(&d_correctCols,      (size_t)refColCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_effectiveRefCols, (size_t)refColCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_singletonRefCols, (size_t)refColCount * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_insertRefTotal,   (size_t)refColCount * sizeof(long long)));

    CUDA_CHECK(cudaMemcpy(d_s, h_s.data(), (size_t)n * k1 * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_charsPerEstimatedColumn, charsPerEstimatedColumn.data(), 
        (size_t)estColCount * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_refSeqs, h_refSeqs.data(), 
        (size_t)n * refColCount * sizeof(char), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemset(d_sharedHom,        0, (size_t)refColCount * sizeof(long long)));
    CUDA_CHECK(cudaMemset(d_totalHom,         0, (size_t)refColCount * sizeof(long long)));
    CUDA_CHECK(cudaMemset(d_correctCols,      0, (size_t)refColCount * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_effectiveRefCols, 0, (size_t)refColCount * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_singletonRefCols, 0, (size_t)refColCount * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_insertRefTotal,   0, (size_t)refColCount * sizeof(long long)));

    /** 
     * KERNEL 1 refColInd
     */
    {
        int blockSize = 128;
        int gridSize  = (n + blockSize - 1) / blockSize;
        computeRefColInd<<<gridSize, blockSize>>>(
            d_refSeqs,
            d_refColInd,
            n,
            refColCount,
            opt.mask_lower_ref);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    /**
     * KERNEL 2 scoreColumns
     */
    int* d_histogram = nullptr;
    CUDA_CHECK(cudaMalloc(&d_histogram, (size_t)refColCount * estColCount * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_histogram, 0, (size_t)refColCount * estColCount * sizeof(int)));
    {
        int blockSize = 512;
        int gridSize  = refColCount;

        scoreColumns<<<gridSize, blockSize>>>(
            d_s, d_refColInd, d_refSeqs, d_charsPerEstimatedColumn,
            n, refColCount, estColCount, k1, opt.mask_lower_ref,
            d_sharedHom, d_totalHom, d_correctCols,
            d_effectiveRefCols, d_singletonRefCols, d_insertRefTotal,
            d_histogram
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }
    cudaFree(d_histogram);

    // reduction via thrust library ==================================================
    auto ll_ptr = [](long long* p) { return thrust::device_ptr<long long>(p); }; // long long
    auto  i_ptr = [](int* p) { return thrust::device_ptr<int>(p); }; // int

    long long sharedHomologies    = thrust::reduce( ll_ptr(d_sharedHom), ll_ptr(d_sharedHom) + refColCount, 0LL);
    long long totalHomologies     = thrust::reduce( ll_ptr(d_totalHom), ll_ptr(d_totalHom) + refColCount, 0LL);
    long long correctColumns      = thrust::reduce( i_ptr(d_correctCols), i_ptr(d_correctCols) + refColCount, 0 );
    long long effectiveRefColumns = thrust::reduce( i_ptr(d_effectiveRefCols), i_ptr(d_effectiveRefCols) + refColCount, 0 );
    long long singletonRefColumns = thrust::reduce( i_ptr(d_singletonRefCols), i_ptr(d_singletonRefCols) + refColCount, 0 );
    long long insertRefTotal      = thrust::reduce( ll_ptr(d_insertRefTotal), ll_ptr(d_insertRefTotal) + refColCount, 0LL );

    CUDA_CHECK(cudaEventRecord(gpuStop));
    CUDA_CHECK(cudaEventSynchronize(gpuStop));

    float gpuMs = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&gpuMs, gpuStart, gpuStop));

    CUDA_CHECK(cudaEventDestroy(gpuStart));
    CUDA_CHECK(cudaEventDestroy(gpuStop));

    // free memory ===================================================================
    cudaFree(d_s);
    cudaFree(d_charsPerEstimatedColumn);
    cudaFree(d_refSeqs);
    cudaFree(d_refColInd);
    cudaFree(d_sharedHom);
    cudaFree(d_totalHom);
    cudaFree(d_correctCols);
    cudaFree(d_effectiveRefCols);
    cudaFree(d_singletonRefCols);
    cudaFree(d_insertRefTotal);

    // swap + results ================================================================
    if (swapped) {
        std::swap(totalHomologies,      totalHomologiesInEsitmated);
        std::swap(effectiveRefColumns,  effectiveEstColumns);
        std::swap(singletonRefColumns,  singletonEstColumns);
        std::swap(insertRefTotal,       insertEstTotal);
    }

    // Diagnostics
    long long cells = (long long)(refColCount + estColCount) * n;
    std::cerr << "MaxLenNoGap= " << k1
              << ", NumSeq= "    << n
              << ", LenRef= "    << (swapped ? estColCount : refColCount)
              << ", LenEst= "    << (swapped ? refColCount : estColCount)
              << ", Cells= "     << cells << "\n";

    FastSPResult r;
    r.sp_score = (double)sharedHomologies / (double)totalHomologies;
    r.modeler  = (double)sharedHomologies / (double)totalHomologiesInEsitmated;
    r.spfn     = 1.0 - r.sp_score;
    r.spfp     = 1.0 - r.modeler;
    r.tc       = (double)correctColumns / (double)effectiveRefColumns;

    r.compression_naive = swapped
        ? (double)refColCount / estColCount
        : (double)estColCount / refColCount;
    r.compression =
        (double)(effectiveEstColumns + insertEstTotal + singletonEstColumns) /
        (double)(effectiveRefColumns + insertRefTotal + singletonRefColumns);

    r.shared_homologies     = sharedHomologies;
    r.total_homologies_ref  = totalHomologies;
    r.total_homologies_est  = totalHomologiesInEsitmated;
    r.correct_columns       = correctColumns;
    r.effective_ref_columns = effectiveRefColumns;
    r.effective_est_columns = effectiveEstColumns;
    r.singleton_ref_columns = singletonRefColumns;
    r.singleton_est_columns = singletonEstColumns;
    r.insert_ref_total      = insertRefTotal;
    r.insert_est_total      = insertEstTotal;
    r.n        = n;
    r.ref_cols = refColCount;
    r.est_cols = estColCount;
    r.k1       = k1;
    r.swapped  = swapped;

    r.gpu_time_ms = gpuMs;

    return r;
}