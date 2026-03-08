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