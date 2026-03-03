#include "fastsp.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

static void print_usage() {
    std::cerr <<
        "FastSP C++ — Siavash Mirarab & Tandy Warnow (C++ port)\n\n"
        "Usage: ./FastSP -r <ref.fa> -e <est.fa> [-o <out.txt>] [-c <gap_chars>] [-ml] [-mlr] [-cuda]\n\n"
        "  -r    Reference alignment (FASTA)\n"
        "  -e    Estimated alignment (FASTA)\n"
        "  -o    Output file (default: stdout)\n"
        "  -c    Additional gap characters (e.g. -c '*~')\n"
        "  -ml   Mask lower-case characters in estimated alignment\n"
        "  -mlr  Mask lower-case characters in reference alignment\n"
        "  -cuda Use GPU (CUDA) implementation\n";
}

static std::string fmt_double(double v) {
    std::ostringstream s;
    s << std::setprecision(16) << v;
    std::string str = s.str();
    if (str.find('.') == std::string::npos && str.find('e') == std::string::npos)
        str += ".0";
    return str;
}

static void print_results(std::ostream& out, const FastSPResult& r) {
    out << "SP-Score "            << fmt_double(r.sp_score)          << "\n";
    out << "Modeler "             << fmt_double(r.modeler)           << "\n";
    out << "SPFN "                << fmt_double(r.spfn)              << "\n";
    out << "SPFP "                << fmt_double(r.spfp)              << "\n";
    out << "Compression (naive) " << fmt_double(r.compression_naive) << "\n";
    out << "Compression "         << fmt_double(r.compression)       << "\n";
    out << "TC "                  << fmt_double(r.tc)                << "\n";
}

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    FastSPOptions opt;
    std::string output_path;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "-r"  && i + 1 < argc) { opt.reference_path = argv[++i]; continue; }
        if (a == "-e"  && i + 1 < argc) { opt.estimated_path = argv[++i]; continue; }
        if (a == "-o"  && i + 1 < argc) { output_path        = argv[++i]; continue; }
        if (a == "-c"  && i + 1 < argc) { opt.gap_chars      = argv[++i]; continue; }
        if (a == "-ml")  { opt.mask_lower_est = true; continue; }
        if (a == "-mlr") { opt.mask_lower_ref = true; continue; }
        if (a == "-cuda") { opt.use_cuda = true; continue; }
        if (a == "-h")   { print_usage(); return 0; }

        std::cerr << "Unknown argument: " << a << "\n";
        print_usage();
        return 2;
    }

    if (opt.reference_path.empty() || opt.estimated_path.empty()) {
        print_usage();
        return 2;
    }

    try {
        const FastSPResult r = run_fastsp(opt);

        if (output_path.empty()) {
            print_results(std::cout, r);
        } else {
            std::ofstream f(output_path);
            if (!f) { std::cerr << "Cannot open output: " << output_path << "\n"; return 1; }
            print_results(f, r);
        }
    } catch (const std::exception& e) {
        std::cerr << "FastSP error: " << e.what() << "\n";
        return 1;
    }

    auto ms = std::chrono::duration<double, std::milli>(
                  std::chrono::high_resolution_clock::now() - start).count();
    std::cerr << "Time: " << ms / 1000.0 << " s\n";
    return 0;
}