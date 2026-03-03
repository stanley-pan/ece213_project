#pragma once
#include <string>

struct FastSPOptions {
    std::string reference_path;
    std::string estimated_path;
    std::string output_path;
    std::string gap_chars;      // extra gap characters beyond '-' and '?'
    bool use_cuda = false;
    bool mask_lower_est = false;
    bool mask_lower_ref = false;
};

struct FastSPResult {
    double sp_score = 0;
    double modeler  = 0;
    double spfn     = 0;
    double spfp     = 0;
    double tc       = 0;

    // Compression
    double compression_naive = 0;
    double compression       = 0;

    // Diagnostic counts
    long long shared_homologies          = 0;
    long long total_homologies_ref       = 0;
    long long total_homologies_est       = 0;
    long long correct_columns            = 0;
    long long effective_ref_columns      = 0;
    long long effective_est_columns      = 0;
    long long singleton_ref_columns      = 0;
    long long singleton_est_columns      = 0;
    long long insert_ref_total           = 0;
    long long insert_est_total           = 0;

    // Alignment dimensions
    int n        = 0; // number of sequences
    int ref_cols = 0;
    int est_cols = 0;
    int k1       = 0; // max non-gap length across sequences

    bool swapped = false;
};

// CPU version (fastsp.cpp)
FastSPResult run_fastsp(const FastSPOptions& opt);

// GPU version (fastsp_cuda.cu)
FastSPResult run_fastsp_cuda(const FastSPOptions& opt);