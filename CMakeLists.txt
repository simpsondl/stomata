cmake_minimum_required(VERSION 3.20)
project(Cystidia VERSION 0.5.0 LANGUAGES CXX CUDA)

# C++17 required
# Use explicit compiler flags for C++ to avoid CMake's feature-mapping
# (workaround for some toolchain/CMake combinations that report
#  "No known features for CXX compiler ...")
add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-std=gnu++17>
    $<$<COMPILE_LANGUAGE:CUDA>:-std=c++17>
)

# CUDA setup
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 86)  # RTX 6000 = compute capability 8.6

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
