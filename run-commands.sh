#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== toolchain ==="
cmake --version
g++ --version
nvcc --version || true

rm -rf build
cmake -S . -B build -DUSE_CUDA=OFF
cmake --build build -j

echo "=== run ==="
# ./build/FastSP -r data/test.true -e data/test.estimated
./build/FastSP -r data/10000/true_aligned.fa -e data/10000/twilight.aln