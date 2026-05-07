# Integration tests

add_executable(test_cli test_cli.cpp)
target_compile_definitions(test_cli PRIVATE
    TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
    STOMATA_EXECUTABLE="$<TARGET_FILE:stomata_cli>"
)
target_link_libraries(test_cli PRIVATE Catch2::Catch2WithMain)
add_test(NAME CliIntegrationTest COMMAND test_cli)
add_dependencies(test_cli stomata_cli)
