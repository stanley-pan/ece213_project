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
./build/FastSP -r data/ref.fa -e data/est.fa