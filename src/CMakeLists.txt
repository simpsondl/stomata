# Cystidia library sources

set(CPU_SOURCES
    bio_utils.cpp
    genome_loader.cpp
    myers_cpu.cpp
    shift_add_cpu.cpp
    verifier.cpp
    search_pipeline.cpp
    cfd_scorer.cpp
)

set(CUDA_SOURCES
    gpu_engine.cu
)

add_library(cystidia ${CPU_SOURCES} ${CUDA_SOURCES})
target_link_libraries(cystidia
    PUBLIC CUDA::cudart Threads::Threads
    PRIVATE spdlog::spdlog ZLIB::ZLIB)

add_executable(cystidia_cli main.cpp)
target_link_libraries(cystidia_cli PRIVATE cystidia)
set_target_properties(cystidia_cli PROPERTIES OUTPUT_NAME "cystidia")
