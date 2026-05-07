# Unit tests

set(UNIT_TESTS
    test_example
    test_genome_loader
    test_myers_cpu
    test_shift_add_cpu
    test_shift_add_gpu
    test_verifier
    test_validation_harness
    test_gpu_engine
    test_search_pipeline
    test_strand_mismatch
    test_bulge_detection
    test_cfd_scorer
)

foreach(T ${UNIT_TESTS})
    add_executable(${T} ${T}.cpp)
    target_compile_definitions(${T} PRIVATE TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data")
    if(${T} STREQUAL "test_example")
        target_link_libraries(${T} PRIVATE Catch2::Catch2WithMain)
    else()
        target_link_libraries(${T} PRIVATE stomata Catch2::Catch2WithMain)
    endif()
    add_test(NAME ${T} COMMAND ${T})
endforeach()
