cmake_minimum_required(VERSION 3.18)
project(FastSP LANGUAGES CXX)

option(USE_CUDA "Build with CUDA support" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRCS
    src/main.cpp
    src/fastsp.cpp
    src/fastsp_cuda.cu     # GPU version
)

if(USE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    list(APPEND SRCS src/fastsp_cuda.cu)

    # Optional: set arch explicitly (A30 = 80). Comment out to let nvcc decide.
    # set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

add_executable(FastSP ${SRCS})

target_compile_options(FastSP PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -Wall -Wextra>
    $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>
)