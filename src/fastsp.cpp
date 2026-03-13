#include "fastsp.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

inline long long c2(long long n) { return n >= 2 ? n * (n - 1) / 2 : 0; }
inline char upper(char c) { return static_cast<char>(c & ~0x20); }  // branchless ASCII upper
inline bool is_lower(char c) { return (unsigned char)(c - 'a') <= 25u; }
inline bool is_upper(char c) { return (unsigned char)(c - 'A') <= 25u; }
inline bool is_alpha(char c) { return is_upper(upper(c)); }

// ---------------------------------------------------------------------------
// Memory-mapped file reader — avoids getline() allocations and copies
// ---------------------------------------------------------------------------
struct MmapFile {
    const char* data = nullptr;
    size_t      size = 0;
    int         fd   = -1;

    explicit MmapFile(const std::string& path) {
        fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("Cannot open: " + path);
        struct stat sb{};
        if (::fstat(fd, &sb) < 0) { ::close(fd); throw std::runtime_error("fstat: " + path); }
        size = (size_t)sb.st_size;
        if (size == 0) { data = ""; return; }
        void* p = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (p == MAP_FAILED) { ::close(fd); throw std::runtime_error("mmap: " + path); }
        ::madvise(p, size, MADV_SEQUENTIAL);
        data = static_cast<const char*>(p);
    }
    ~MmapFile() {
        if (data && size) ::munmap(const_cast<char*>(data), size);
        if (fd >= 0) ::close(fd);
    }
    MmapFile(const MmapFile&) = delete;
    MmapFile& operator=(const MmapFile&) = delete;
};

// ---------------------------------------------------------------------------
// Flat FASTA alignment: sequences stored as contiguous chars, column-major
// friendly via (seq * cols + col) indexing
// ---------------------------------------------------------------------------
struct Alignment {
    std::vector<char>        seqs;        // flat [n * cols]
    std::vector<std::string> names;
    std::unordered_map<std::string, int> name_to_idx;
    int cols = 0;
    int n    = 0;

    char get(int i, int j) const { return seqs[(size_t)i * cols + j]; }
};

Alignment read_alignment(const std::string& path) {
    MmapFile f(path);
    const char* p   = f.data;
    const char* end = f.data + f.size;

    Alignment aln;
    // --- first pass: count sequences and determine column count ---
    {
        const char* q = p;
        while (q < end) {
            if (*q == '>') { aln.n++; while (q < end && *q != '\n') q++; }
            q++;
        }
    }
    if (aln.n == 0) throw std::runtime_error("No sequences found in: " + path);

    aln.names.reserve(aln.n);

    // --- second pass: parse ---
    // We don't know cols yet; collect each sequence into a tmp buffer then copy.
    // Only one temp vector needed — reused per sequence.
    std::vector<char> tmp;
    tmp.reserve(1 << 20);

    auto flush = [&](std::string& name) {
        if (name.empty()) return;
        int idx = (int)aln.names.size();
        aln.name_to_idx[name] = idx;
        aln.names.push_back(std::move(name));
        name.clear();

        if (idx == 0) {
            aln.cols = (int)tmp.size();
            aln.seqs.resize((size_t)aln.n * aln.cols);
        } else if ((int)tmp.size() != aln.cols) {
            throw std::runtime_error("Inconsistent row lengths in: " + path);
        }
        std::memcpy(aln.seqs.data() + (size_t)idx * aln.cols, tmp.data(), aln.cols);
        tmp.clear();
    };

    std::string cur_name;
    const char* q = p;
    while (q < end) {
        if (*q == '>') {
            flush(cur_name);
            q++;  // skip '>'
            // read name up to whitespace / newline
            const char* ns = q;
            while (q < end && *q != ' ' && *q != '\t' && *q != '\r' && *q != '\n') q++;
            cur_name.assign(ns, q);
            while (q < end && *q != '\n') q++;  // consume rest of header line
            if (q < end) q++;                   // skip '\n'
        } else if (*q == '\n' || *q == '\r' || *q == ' ' || *q == '\t') {
            q++;
        } else {
            tmp.push_back(*q++);
        }
    }
    flush(cur_name);
    return aln;
}

} // namespace

// ---------------------------------------------------------------------------
FastSPResult run_fastsp(const FastSPOptions& opt) {

    // Build gap-char lookup table (256-entry bool array — O(1), branch-free)
    bool is_gap_table[256] = {};
    is_gap_table[(unsigned char)'-'] = true;
    is_gap_table[(unsigned char)'?'] = true;
    for (char c : opt.gap_chars) is_gap_table[(unsigned char)c] = true;
    auto is_gap = [&](char c) __attribute__((always_inline)) {
        return is_gap_table[(unsigned char)c];
    };

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

    const int n           = ref.n;
    const int refColCount = ref.cols;
    const int estColCount = est.cols;

    if (est.n != n)
        throw std::runtime_error("Ref and est have different sequence counts.");

    // Build row-order map: for each ref row i, which est row index?
    // Store a remapping so we can access est rows in ref order without copying.
    std::vector<int> est_row_for_ref(n);
    for (int i = 0; i < n; i++) {
        auto it = est.name_to_idx.find(ref.names[i]);
        if (it == est.name_to_idx.end())
            throw std::runtime_error("Sequence missing from estimated alignment: " + ref.names[i]);
        est_row_for_ref[i] = it->second;
    }

    // ---------------------------------------------------------------------------
    // Count k1 (max non-gap length across all sequences in estimated, in ref order)
    // ---------------------------------------------------------------------------
    int k1 = 0;
    for (int i = 0; i < n; i++) {
        int er = est_row_for_ref[i];
        int cnt = 0;
        const char* row = est.seqs.data() + (size_t)er * estColCount;
        for (int j = 0; j < estColCount; j++) {
            char c = row[j];
            char C = opt.mask_lower_est ? c : upper(c);
            if (!is_gap(C) && is_alpha(C)) cnt++;
        }
        k1 = std::max(k1, cnt);
    }
    if (k1 <= 0) throw std::runtime_error("Estimated alignment has no letters.");

    // ---------------------------------------------------------------------------
    // Build s[i * k1 + k] = estimated column of k-th non-gap char in seq i
    // FLAT array — one allocation, cache-friendly
    // ---------------------------------------------------------------------------
    std::vector<int>  s_flat(static_cast<size_t>(n) * k1, 0);
    std::vector<int>  charsPerEstimatedColumn(estColCount, 0);
    long long insertEstTotal = 0;

    for (int i = 0; i < n; i++) {
        int er = est_row_for_ref[i];
        const char* row = est.seqs.data() + (size_t)er * estColCount;
        int* s_row = s_flat.data() + (size_t)i * k1;
        int nucInd = 0;
        for (int j = 0; j < estColCount; j++) {
            char c = row[j];
            char C = opt.mask_lower_est ? c : upper(c);
            if (is_gap(C) || !is_alpha(C)) continue;
            if (nucInd >= k1) throw std::runtime_error("Non-gap overflow — est parsing bug.");
            if (opt.mask_lower_est && is_lower(c)) {
                s_row[nucInd] = -1;
                insertEstTotal++;
            } else {
                s_row[nucInd] = j;
                charsPerEstimatedColumn[j]++;
            }
            nucInd++;
        }
    }

    long long totalHomologiesInEsitmated = 0;
    long long effectiveEstColumns = 0;
    long long singletonEstColumns = 0;
    for (int j = 0; j < estColCount; j++) {
        int x = charsPerEstimatedColumn[j];
        totalHomologiesInEsitmated += c2(x);
        if      (x == 1) singletonEstColumns++;
        else if (x >= 2) effectiveEstColumns++;
    }

    std::vector<int> refColInd(n, 0);

    std::vector<int>  estSitesCnt(estColCount, 0);
    std::vector<int>  touched;
    touched.reserve(n);

    long long sharedHomologies    = 0;
    long long totalHomologies     = 0;
    long long correctColumns      = 0;
    long long effectiveRefColumns = 0;
    long long singletonRefColumns = 0;
    long long insertRefTotal      = 0;

    for (int c = 0; c < refColCount; c++) {
        long long refCharCount = 0;
        int  singleEstCol = -2;   // used for TC check: -2=uninit, -1=multiple
        bool tcCandidate  = true;

        for (int i = 0; i < n; i++) {
            char rc = ref.get(i, c);
            char RC = opt.mask_lower_ref ? rc : upper(rc);

            if (is_gap(RC) || !is_alpha(RC)) continue;

            if (opt.mask_lower_ref && is_lower(rc)) {
                insertRefTotal++;
                refColInd[i]++;
                continue;
            }

            refCharCount++;
            int y = s_flat[(size_t)i * k1 + refColInd[i]];
            if (y != -1) {
                if (estSitesCnt[y]++ == 0) touched.push_back(y);
                // Track for TC: all residues must land in the same est column
                if (singleEstCol == -2)      singleEstCol = y;
                else if (singleEstCol != y)  singleEstCol = -1;
            } else {
                tcCandidate = false;  // residue unaligned in est → can't be correct col
                singleEstCol = -1;
            }

            refColInd[i]++;
        }

        // Accumulate shared homologies and reset touched slots
        for (int y : touched)
            sharedHomologies += c2(estSitesCnt[y]);

        totalHomologies += c2(refCharCount);

        // TC: all ref chars → same est col, and that est col has no extras
        if (tcCandidate && refCharCount >= 2 && singleEstCol >= 0 &&
            (long long)charsPerEstimatedColumn[singleEstCol] == refCharCount &&
            (long long)estSitesCnt[singleEstCol] == refCharCount)
        {
            correctColumns++;
        }

        if      (refCharCount == 1) singletonRefColumns++;
        else if (refCharCount >= 2) effectiveRefColumns++;

        // Reset only touched slots — O(touched) not O(estColCount)
        for (int y : touched) estSitesCnt[y] = 0;
        touched.clear();
    }

    // ---------------------------------------------------------------------------
    // Swap-back
    // ---------------------------------------------------------------------------
    if (swapped) {
        std::swap(totalHomologies, totalHomologiesInEsitmated);
        std::swap(effectiveRefColumns, effectiveEstColumns);
        std::swap(singletonRefColumns, singletonEstColumns);
        std::swap(insertRefTotal, insertEstTotal);
    }

    if (totalHomologies == 0)
        throw std::runtime_error("Reference has zero homologies.");
    if (totalHomologiesInEsitmated == 0)
        throw std::runtime_error("Estimated has zero homologies.");
    if (effectiveRefColumns == 0)
        throw std::runtime_error("Reference has no columns with >=2 letters.");

    long long cells = (long long)(refColCount + estColCount) * n;
    std::cerr << "MaxLenNoGap= " << k1
              << ", NumSeq= "    << n
              << ", LenRef= "    << (swapped ? estColCount : refColCount)
              << ", LenEst= "    << (swapped ? refColCount : estColCount)
              << ", Cells= "     << cells << "\n";
    std::cerr << "Number of shared homologies:         " << sharedHomologies           << "\n";
    std::cerr << "Number of homologies in reference:   " << totalHomologies            << "\n";
    std::cerr << "Number of homologies in estimated:   " << totalHomologiesInEsitmated << "\n";
    std::cerr << "Number of correctly aligned columns: " << correctColumns             << "\n";
    std::cerr << "Number of aligned columns in ref:    " << effectiveRefColumns        << "\n";
    std::cerr << "Number of aligned columns in est:    " << effectiveEstColumns        << "\n";
    if (opt.mask_lower_ref || singletonRefColumns != 0)
        std::cerr << "Singleton/insertion columns in ref:  "
                  << singletonRefColumns << " " << insertRefTotal << "\n";
    if (opt.mask_lower_est || singletonEstColumns != 0)
        std::cerr << "Singleton/insertion columns in est:  "
                  << singletonEstColumns << " " << insertEstTotal << "\n";

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