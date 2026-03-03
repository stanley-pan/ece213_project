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

// helpers
namespace {

inline long long c2(long long n) { return n >= 2 ? n * (n - 1) / 2 : 0; }

inline char upper(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
inline bool is_lower(char c) { return std::islower(static_cast<unsigned char>(c)) != 0; }
inline bool is_alpha(char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0; }

// ---------------------------------------------------------------------------
// FASTA loader — whole file, multi-line sequences supported
// ---------------------------------------------------------------------------
struct Alignment {
    std::vector<std::string>            names;
    std::vector<std::string>            seqs;
    std::unordered_map<std::string,int> name_to_idx;
    int cols = 0;
    int n    = 0;
};

Alignment read_fasta(const std::string& path) {
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

std::unordered_set<char> make_gap_set(const std::string& extras) {
    std::unordered_set<char> g = {'-', '?'};
    for (char c : extras) g.insert(c);
    return g;
}

} // namespace

FastSPResult run_fastsp(const FastSPOptions& opt) {

    const auto gaps   = make_gap_set(opt.gap_chars);
    auto       is_gap = [&](char c) { return gaps.count(c) > 0; };

    // ---- Load ---------------------------------------------------------------
    std::cerr << "Reference alignment: " << opt.reference_path << " ...\n";
    std::cerr << "Estimated alignment: " << opt.estimated_path << " ...\n";

    Alignment ref = read_fasta(opt.reference_path);
    Alignment est = read_fasta(opt.estimated_path);

    bool swapped = false;
    if (est.cols < ref.cols && !opt.mask_lower_est && !opt.mask_lower_ref) {
        swapped = true;
        std::swap(ref, est);
        std::cerr << "(alignments swapped: estimated was shorter)\n";
    }

    const int n       = ref.n;
    const int refCols = ref.cols;
    const int estCols = est.cols;

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

    // ---- Count max non-gap length (k1) across all sequences -----------------
    int k1 = 0;
    for (int i = 0; i < n; i++) {
        int ng = 0;
        for (char c : est_ordered[i]) {
            char C = opt.mask_lower_est ? c : upper(c);
            if (!is_gap(C) && is_alpha(C)) ng++;
        }
        k1 = std::max(k1, ng);
    }
    if (k1 <= 0) throw std::runtime_error("Estimated alignment has no letters.");

    // ---- Build s[i][k] = estimated column of k-th non-gap char in sequence i -
    std::vector<std::vector<int>> s(n, std::vector<int>(k1, 0));
    std::vector<int> chars_per_est_col(estCols, 0);
    long long insert_est_total = 0;

    for (int i = 0; i < n; i++) {
        int nuc = 0;
        for (int j = 0; j < estCols; j++) {
            char c = est_ordered[i][j];
            char C = opt.mask_lower_est ? c : upper(c);

            if (is_gap(C) || !is_alpha(C)) continue;
            if (nuc >= k1) throw std::runtime_error("Non-gap overflow — est parsing bug.");

            if (opt.mask_lower_est && is_lower(c)) {
                s[i][nuc] = -1;      // insertion: treated as unaligned
                insert_est_total++;
            } else {
                s[i][nuc] = j;
                chars_per_est_col[j]++;
            }
            nuc++;
        }
    }

    // ---- Estimated homology stats -------------------------------------------
    long long total_hom_est      = 0;
    long long effective_est_cols = 0;
    long long singleton_est_cols = 0;

    for (int j = 0; j < estCols; j++) {
        total_hom_est += c2(chars_per_est_col[j]);
        if      (chars_per_est_col[j] == 1) singleton_est_cols++;
        else if (chars_per_est_col[j] >= 2) effective_est_cols++;
    }

    // ---- Scan reference columns ---------------------------------------------
    std::vector<int> ref_col_idx(n, 0);

    long long shared_hom         = 0;
    long long total_hom_ref      = 0;
    long long correct_cols       = 0;
    long long effective_ref_cols = 0;
    long long singleton_ref_cols = 0;
    long long insert_ref_total   = 0;

    std::unordered_map<int,int> est_sites_count;
    est_sites_count.reserve(n);

    for (int c = 0; c < refCols; c++) {
        est_sites_count.clear();
        long long ref_char_count = 0;

        for (int i = 0; i < n; i++) {
            char rc = ref.seqs[i][c];
            char RC = opt.mask_lower_ref ? rc : upper(rc);

            if (is_gap(RC) || !is_alpha(RC)) continue;

            // Lowercase in reference = insertion site when mask_lower_ref enabled
            if (opt.mask_lower_ref && is_lower(rc)) {
                insert_ref_total++;
                ref_col_idx[i]++;
                continue;
            }

            ref_char_count++;
            int est_col = s[i][ref_col_idx[i]];
            if (est_col != -1)
                est_sites_count[est_col]++;

            ref_col_idx[i]++;
        }

        for (const auto& [col, cnt] : est_sites_count)
            shared_hom += c2(cnt);

        total_hom_ref += c2(ref_char_count);

        // TC: correct column = all ref chars map to same est col with no extras
        if (ref_char_count >= 2 && est_sites_count.size() == 1) {
            const auto& [est_col, cnt] = *est_sites_count.begin();
            if (cnt == ref_char_count && chars_per_est_col[est_col] == (int)ref_char_count)
                correct_cols++;
        }

        if      (ref_char_count == 1) singleton_ref_cols++;
        else if (ref_char_count >= 2) effective_ref_cols++;
    }

    if (swapped) {
        std::swap(total_hom_ref,      total_hom_est);
        std::swap(effective_ref_cols, effective_est_cols);
        std::swap(singleton_ref_cols, singleton_est_cols);
        std::swap(insert_ref_total,   insert_est_total);
    }

    // ---- Sanity checks ------------------------------------------------------
    if (total_hom_ref == 0)
        throw std::runtime_error("Reference has zero homologies — is this really an alignment?");
    if (total_hom_est == 0)
        throw std::runtime_error("Estimated has zero homologies — is this really an alignment?");
    if (effective_ref_cols == 0)
        throw std::runtime_error("Reference has no columns with >=2 letters.");

    // ---- Stderr diagnostics -------------------------------------------------
    long long cells = (long long)(refCols + estCols) * n;
    std::cerr << "MaxLenNoGap= " << k1
              << ", NumSeq= "    << n
              << ", LenRef= "    << (swapped ? estCols : refCols)
              << ", LenEst= "    << (swapped ? refCols : estCols)
              << ", Cells= "     << cells << "\n";
    std::cerr << "Number of shared homologies:         " << shared_hom         << "\n";
    std::cerr << "Number of homologies in reference:   " << total_hom_ref      << "\n";
    std::cerr << "Number of homologies in estimated:   " << total_hom_est      << "\n";
    std::cerr << "Number of correctly aligned columns: " << correct_cols       << "\n";
    std::cerr << "Number of aligned columns in ref:    " << effective_ref_cols << "\n";
    std::cerr << "Number of aligned columns in est:    " << effective_est_cols << "\n";
    if (opt.mask_lower_ref || singleton_ref_cols != 0)
        std::cerr << "Singleton/insertion columns in ref:  "
                  << singleton_ref_cols << " " << insert_ref_total << "\n";
    if (opt.mask_lower_est || singleton_est_cols != 0)
        std::cerr << "Singleton/insertion columns in est:  "
                  << singleton_est_cols << " " << insert_est_total << "\n";

    // ---- Build result -------------------------------------------------------
    FastSPResult r;
    r.sp_score = (double)shared_hom   / (double)total_hom_ref;
    r.modeler  = (double)shared_hom   / (double)total_hom_est;
    r.spfn     = 1.0 - r.sp_score;
    r.spfp     = 1.0 - r.modeler;
    r.tc       = (double)correct_cols / (double)effective_ref_cols;

    r.compression_naive = swapped
        ? (double)refCols / estCols
        : (double)estCols / refCols;
    r.compression =
        (double)(effective_est_cols + insert_est_total + singleton_est_cols) /
        (double)(effective_ref_cols + insert_ref_total + singleton_ref_cols);

    r.shared_homologies     = shared_hom;
    r.total_homologies_ref  = total_hom_ref;
    r.total_homologies_est  = total_hom_est;
    r.correct_columns       = correct_cols;
    r.effective_ref_columns = effective_ref_cols;
    r.effective_est_columns = effective_est_cols;
    r.singleton_ref_columns = singleton_ref_cols;
    r.singleton_est_columns = singleton_est_cols;
    r.insert_ref_total      = insert_ref_total;
    r.insert_est_total      = insert_est_total;
    r.n        = n;
    r.ref_cols = refCols;
    r.est_cols = estCols;
    r.k1       = k1;
    r.swapped  = swapped;

    return r;
}