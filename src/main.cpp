#include "fastsp.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <thread>

static void print_usage() {
    std::cerr <<
        "ERROR: ./FastSP -r <ref.fa> -e <est.fa> [-ml] [-mlr]\n";
}

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    FastSPOptions opt; // options for arg

    // parse args
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "-r" && i + 1 < argc) { opt.reference_path = argv[++i]; continue; }
        if (a == "-e" && i + 1 < argc) { opt.estimated_path = argv[++i]; continue; }
        if (a == "-o" && i + 1 < argc) { opt.output_path    = argv[++i]; continue; }
        if (a == "-ml")  { opt.mask_lower_est = true; continue; }
        if (a == "-mlr") { opt.mask_lower_ref = true; continue; }

        std::cerr << "Unknown/invalid arg: " << a << "\n";
        print_usage();
        return 2;
    }

    if (opt.reference_path.empty() || opt.estimated_path.empty()) {
        print_usage();
        return 2;
    }

    // run program
    try {
        FastSPResult res;

        // print results
        std::cout << "SP-Score: " << res.sp_score << "\n";
        std::cout << "Modeler: "  << res.modeler  << "\n";
        std::cout << "SPFN: "     << res.spfn     << "\n";
        std::cout << "SPFP: "     << res.spfp     << "\n";
        std::cout << "TC: "       << res.tc       << "\n";
    } catch (const std::exception& e) {
        std::cerr << "FastSP error: " << e.what() << "\n";
        return 1;
    }


    // check runtime
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Runtime: " << duration.count() << " ms\n";

    return 0;
}