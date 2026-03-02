#include "fastsp.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

int main(int argc, char** argv) {
    std::cout << "FastSP OK. argc=" << argc << "\n";
    for (int i = 0; i < argc; i++) std::cout << "argv[" << i << "]=" << argv[i] << "\n";
    return 0;
}