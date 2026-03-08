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

// CUDA error checking
#define CUDA_CHECK(call)                                         \
    do {                                                         \
        cudaError_t _e = (call);                                 \
        if (_e != cudaSuccess)                                   \
            throw std::runtime_error(std::string("CUDA error: ") \
                                     + cudaGetErrorString(_e));  \
    } while (0)

// HELPERS =====================================================================
namespace {

    inline long long c2(long long n) { return n >= 2 ? n * (n - 1) / 2 : 0; }

    inline char upper(char c) { 
        return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
    inline bool is_lower(char c) { return std::islower(static_cast<unsigned char>(c)) != 0; }
    inline bool is_alpha(char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0; }

    // helper to read alignments
    class Alignment {
        public:
            std::vector<std::string> names;
            std::vector<std::string> seqs;
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
            aln.seqs.push_back(cur_seq);
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
    std::vector<int> h_s(n * k1, 0);
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

    // Flatten reference sequences
    std::vector<char> h_refSeqs(n * refColCount);
    for (int i = 0; i < n; i++)
        for (int c = 0; c < refColCount; c++)
            h_refSeqs[i * refColCount + c] = ref.seqs[i][c];

    // TODO: allocate mem

    // TODO: 1) compute refColInd prefix sums
    // TODO: 2) score every ref column
    // TODO: 3) reduce per-column arrays (use thrust lib?)


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

    // results
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

    return r;
}