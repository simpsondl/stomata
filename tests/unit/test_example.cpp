#include <catch2/catch_test_macros.hpp>

// This is a dummy test to verify Catch2 is working
TEST_CASE("Example test", "[example]") {
    REQUIRE(1 + 1 == 2);
    REQUIRE(true);
}

TEST_CASE("Example test with sections", "[example]") {
    SECTION("Addition") {
        REQUIRE(2 + 2 == 4);
    }
    
    SECTION("Subtraction") {
        REQUIRE(5 - 3 == 2);
    }
}
