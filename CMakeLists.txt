cmake_minimum_required(VERSION 3.20)
project(Stomata VERSION 0.11.0 LANGUAGES CXX CUDA)

# C++17 required
# Use explicit compiler flags for C++ to avoid CMake's feature-mapping
# (workaround for some toolchain/CMake combinations that report
#  "No known features for CXX compiler ...")
add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-std=gnu++17>
    $<$<COMPILE_LANGUAGE:CUDA>:-std=c++17>
)

# CUDA architecture targets.
# Default builds a fat binary covering Volta → Hopper, so the same binary
# runs on V100 (sm_70), GTX 16x0 / RTX 20x0 / T4 (sm_75), A100 (sm_80),
# RTX 30x0 / A40 (sm_86), RTX 40x0 / RTX 6000 Ada (sm_89), and H100 (sm_90).
# Override at configure time for a smaller, faster-linking binary, e.g.:
#   cmake -DCMAKE_CUDA_ARCHITECTURES=89 ...
# (single-arch RTX 6000 Ada build, what the dev machine has used historically).
set(CMAKE_CUDA_STANDARD 17)
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "70;75;80;86;89;90")
endif()

# Find dependencies
find_package(spdlog REQUIRED)
find_package(ZLIB REQUIRED)
find_package(CUDAToolkit REQUIRED)
find_package(Threads REQUIRED)
find_package(Catch2 3 QUIET)

# Include directories
include_directories(include)

# Subdirectories
add_subdirectory(src)
if(Catch2_FOUND)
    enable_testing()
    add_subdirectory(tests)
else()
    message(WARNING "Catch2 not found. Tests will not be built. Install with: conda install -c conda-forge catch2")
endif()

