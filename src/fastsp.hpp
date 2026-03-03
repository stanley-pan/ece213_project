#pragma once
#include <string>

struct FastSPOptions {
    std::string reference_path;
    std::string estimated_path;
    std::string output_path;   // TODO: set output path
    bool mask_lower_est = false;
    bool mask_lower_ref = false;
    bool use_cuda = false;
};

struct FastSPResult {
    double sp_score = 0;
    double modeler  = 0;
    double spfn     = 0;
    double spfp     = 0;
    double tc       = 0;
};

FastSPResult run_fastsp(const FastSPOptions& opt);