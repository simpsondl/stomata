# Stomata library sources

set(CPU_SOURCES
    annotation_bed.cpp
    bio_utils.cpp
    cfd_scorer.cpp
    genome_loader.cpp
    myers_cpu.cpp
    search_pipeline.cpp
    shift_add_cpu.cpp
    spacer_summary.cpp
    verifier.cpp
)

set(CUDA_SOURCES
    gpu_engine.cu
)

add_library(stomata ${CPU_SOURCES} ${CUDA_SOURCES})
target_link_libraries(stomata
    PUBLIC CUDA::cudart Threads::Threads
    PRIVATE spdlog::spdlog ZLIB::ZLIB)

add_executable(stomata_cli main.cpp)
target_link_libraries(stomata_cli PRIVATE stomata)
set_target_properties(stomata_cli PROPERTIES OUTPUT_NAME "stomata")
