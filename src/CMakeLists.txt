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

if(STOMATA_ENABLE_CUDA)
    set(GPU_SOURCES gpu_engine.cu)
else()
    set(GPU_SOURCES gpu_engine_stub.cpp)
endif()

add_library(stomata ${CPU_SOURCES} ${GPU_SOURCES})
target_link_libraries(stomata
    PUBLIC Threads::Threads
    PRIVATE spdlog::spdlog ZLIB::ZLIB)
if(STOMATA_ENABLE_CUDA)
    target_link_libraries(stomata PUBLIC CUDA::cudart)
endif()

add_executable(stomata_cli main.cpp)
target_link_libraries(stomata_cli PRIVATE stomata)
set_target_properties(stomata_cli PROPERTIES OUTPUT_NAME "stomata")
include(GNUInstallDirs)
install(TARGETS stomata_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
