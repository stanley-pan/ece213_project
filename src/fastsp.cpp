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

namespace {

    // helpers

    // computes number of pairs that can be formed
    inline long long c2(long long n) { return n >= 2 ? n * (n - 1) / 2 : 0; }

    inline char upper(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
    inline bool is_lower(char c) { return std::islower(static_cast<unsigned char>(c)) != 0; }
    inline bool is_alpha(char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0; }

    // FASTA alignment reader
    class Alignment {
        public:
            std::vector<std::string> names;
            std::vector<std::string> seqs;
            std::unordered_map<std::string,int> name_to_idx;
            int cols = 0;
            int n    = 0;
    };

    // reads sequences
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
            cur_name.clear();
            cur_seq .clear();
        };

        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                flush();
                // name = everything after '>' up to first whitespace
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
                throw std::runtime_error("Inconsistent row lengths in: " + path +
                                        " (sequence: " + aln.names[i] + ")");
        return aln;
    }

    std::unordered_set<char> dash_chars(const std::string& extras) {
        std::unordered_set<char> g = {'-', '?'};
        for (char c : extras) g.insert(c);
        return g;
    }

} // namespace

// main function
FastSPResult run_fastsp(const FastSPOptions& opt) {

    const auto dashChars = dash_chars(opt.gap_chars);
    auto       is_gap    = [&](char c) { return dashChars.count(c) > 0; };

    // Load
    std::cerr << "Reference alignment: " << opt.reference_path << " ...\n";
    std::cerr << "Estimated alignment: " << opt.estimated_path << " ...\n";

    Alignment ref = read_alignment(opt.reference_path);
    Alignment est = read_alignment(opt.estimated_path);

    bool swapped = false;

    // if estimated is shorter, swap and make it reference
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

    // Reorder est rows to match reference order, validate names match
    std::vector<std::string> est_ordered(n);
    for (int i = 0; i < n; i++) {
        auto it = est.name_to_idx.find(ref.names[i]);
        if (it == est.name_to_idx.end())
            throw std::runtime_error("Sequence missing from estimated alignment: " + ref.names[i]);
        est_ordered[i] = est.seqs[it->second];
    }

    // Count max non-gap length (k1) across all sequences ======================
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

    // Build s[i][k] = estimated column of k-th non-gap char in sequence i =====
    std::vector<std::vector<int>> s(n, std::vector<int>(k1, 0));
    std::vector<int> charsPerEstimatedColumn(estColCount, 0);
    long long insertEstTotal = 0;

    for (int i = 0; i < n; i++) {
        int nucInd = 0;
        for (int j = 0; j < estColCount; j++) {
            char c = est_ordered[i][j];
            char C = opt.mask_lower_est ? c : upper(c);

            if (is_gap(C) || !is_alpha(C)) continue;
            if (nucInd >= k1) throw std::runtime_error("Non-gap overflow — est parsing bug.");

            if (opt.mask_lower_est && is_lower(c)) {
                s[i][nucInd] = -1;      // insertion: treated as unaligned
                insertEstTotal++;
            } else {
                s[i][nucInd] = j;
                charsPerEstimatedColumn[j]++;
            }
            nucInd++;
        }
    }

    // Estimated homology stats ================================================
    long long totalHomologiesInEsitmated = 0;
    long long effectiveEstColumns = 0;
    long long singletonEstColumns = 0;

    for (int j = 0; j < estColCount; j++) {
        totalHomologiesInEsitmated += c2(charsPerEstimatedColumn[j]);
        if      (charsPerEstimatedColumn[j] == 1) singletonEstColumns++;
        else if (charsPerEstimatedColumn[j] >= 2) effectiveEstColumns++;
    }

    // Scan reference columns ==================================================
    std::vector<int> refColInd(n, 0);

    long long sharedHomologies    = 0;
    long long totalHomologies     = 0;
    long long correctColumns      = 0;
    long long effectiveRefColumns = 0;
    long long singletonRefColumns = 0;
    long long insertRefTotal      = 0;

    std::unordered_map<int,int> estimatedSitesCount;
    estimatedSitesCount.reserve(n);

    for (int c = 0; c < refColCount; c++) {
        estimatedSitesCount.clear();
        long long refCharCount = 0;

        for (int i = 0; i < n; i++) {
            char rc = ref.seqs[i][c];
            char RC = opt.mask_lower_ref ? rc : upper(rc);

            if (is_gap(RC) || !is_alpha(RC)) continue;

            // Lowercase in reference = insertion site when mask_lower_ref enabled
            if (opt.mask_lower_ref && is_lower(rc)) {
                insertRefTotal++;
                refColInd[i]++;
                continue;
            }

            refCharCount++;
            int y = s[i][refColInd[i]];
            if (y != -1)
                estimatedSitesCount[y]++;

            refColInd[i]++;
        }

        for (const auto& [y, cnt] : estimatedSitesCount)
            sharedHomologies += c2(cnt);

        totalHomologies += c2(refCharCount);

        // TC: correct column = all ref chars map to same est col with no extras
        if (refCharCount >= 2 && estimatedSitesCount.size() == 1) {
            const auto& [estColumn, cnt] = *estimatedSitesCount.begin();
            if (cnt == refCharCount && charsPerEstimatedColumn[estColumn] == (int)refCharCount)
                correctColumns++;
        }

        if      (refCharCount == 1) singletonRefColumns++;
        else if (refCharCount >= 2) effectiveRefColumns++;
    }

    if (swapped) {
        std::swap(totalHomologies, totalHomologiesInEsitmated);
        std::swap(effectiveRefColumns, effectiveEstColumns);
        std::swap(singletonRefColumns, singletonEstColumns);
        std::swap(insertRefTotal, insertEstTotal);
    }

    // SANITY CHECK!!! =========================================================
    if (totalHomologies == 0)
        throw std::runtime_error("Reference has zero homologies — is this really an alignment?");
    if (totalHomologiesInEsitmated == 0)
        throw std::runtime_error("Estimated has zero homologies — is this really an alignment?");
    if (effectiveRefColumns == 0)
        throw std::runtime_error("Reference has no columns with >=2 letters.");

    // Stderr ==================================================================
    long long cells = (long long)(refColCount + estColCount) * n;
    std::cerr << "MaxLenNoGap= " << k1
              << ", NumSeq= "    << n
              << ", LenRef= "    << (swapped ? estColCount : refColCount)
              << ", LenEst= "    << (swapped ? refColCount : estColCount)
              << ", Cells= "     << cells << "\n";
    std::cerr << "Number of shared homologies:         " << sharedHomologies            << "\n";
    std::cerr << "Number of homologies in reference:   " << totalHomologies             << "\n";
    std::cerr << "Number of homologies in estimated:   " << totalHomologiesInEsitmated  << "\n";
    std::cerr << "Number of correctly aligned columns: " << correctColumns              << "\n";
    std::cerr << "Number of aligned columns in ref:    " << effectiveRefColumns         << "\n";
    std::cerr << "Number of aligned columns in est:    " << effectiveEstColumns         << "\n";
    if (opt.mask_lower_ref || singletonRefColumns != 0)
        std::cerr << "Singleton/insertion columns in ref:  "
                  << singletonRefColumns << " " << insertRefTotal << "\n";
    if (opt.mask_lower_est || singletonEstColumns != 0)
        std::cerr << "Singleton/insertion columns in est:  "
                  << singletonEstColumns << " " << insertEstTotal << "\n";

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