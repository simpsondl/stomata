/**
 * Phase 2 Validation Harness
 *
 * Bulk random testing of Myers bit-parallel algorithm against the O(mn) DP verifier.
 * This ensures Myers produces identical results to the ground-truth reference
 * before GPU porting.
 */

#include <catch2/catch_test_macros.hpp>
#include <myers_cpu.hpp>
#include <verifier.hpp>

#include <random>
#include <string>
#include <vector>

// Fixed seed for reproducibility - identical sequences across all runs
static constexpr uint64_t VALIDATION_SEED = 0xDEADBEEF42;

// Generate a random nucleotide string of given length
static std::string random_sequence(std::mt19937_64& rng, size_t length,
                                   double n_probability = 0.02) {
    static const char BASES[] = "ACGT";
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int> base_dist(0, 3);

    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        if (prob(rng) < n_probability) {
            seq += 'N';
        } else {
            seq += BASES[base_dist(rng)];
        }
    }
    return seq;
}

// Generate length in range [min, max] with uniform distribution
static size_t random_length(std::mt19937_64& rng, size_t min, size_t max) {
    std::uniform_int_distribution<size_t> dist(min, max);
    return dist(rng);
}

// Test category definition
struct TestCategory {
    const char* name;
    size_t pattern_min, pattern_max;
    size_t text_min, text_max;
    size_t num_tests;
    double n_probability;
};

// Stratified test categories covering single-word and multi-word patterns
static const TestCategory CATEGORIES[] = {
    // name                pat_min  pat_max  text_min  text_max  count  n_prob
    {"tiny",                    1,       8,        1,       20,   500,  0.05},
    {"short_single_word",       9,      32,       50,      200,   300,  0.02},
    {"full_single_word",       33,      64,      100,      500,   200,  0.02},
    {"boundary_64",            63,      65,      200,      400,   100,  0.02},
    {"multi_word_2",           65,     128,      200,     1000,   100,  0.02},
    {"multi_word_3",          129,     192,      300,     1500,    50,  0.02},
    {"multi_word_large",      193,     255,      500,     2000,    25,  0.01},
};

// =============================================================================
// TEST CASE 1: Bulk random comparison across all categories
// =============================================================================

TEST_CASE("Validation - Myers matches Verifier on random inputs",
          "[validation][random]") {

    for (const auto& cat : CATEGORIES) {
        SECTION(cat.name) {
            std::mt19937_64 rng(VALIDATION_SEED);

            size_t failures = 0;
            size_t first_failure_idx = 0;
            std::string first_failure_pattern, first_failure_text;

            for (size_t i = 0; i < cat.num_tests; ++i) {
                size_t pat_len = random_length(rng, cat.pattern_min, cat.pattern_max);
                size_t txt_len = random_length(rng, cat.text_min, cat.text_max);

                std::string pattern = random_sequence(rng, pat_len, cat.n_probability);
                std::string text = random_sequence(rng, txt_len, cat.n_probability);

                auto myers = myers_edit_distances(pattern, text);
                auto verifier = verify_edit_distances(pattern, text);

                if (myers.distances != verifier.distances) {
                    if (failures == 0) {
                        first_failure_idx = i;
                        first_failure_pattern = pattern;
                        first_failure_text = text;
                    }
                    ++failures;
                }
            }

            INFO("Category: " << cat.name);
            INFO("Tests run: " << cat.num_tests);
            INFO("Failures: " << failures);
            if (failures > 0) {
                INFO("First failure at index: " << first_failure_idx);
                INFO("Pattern (" << first_failure_pattern.size() << " bp): "
                     << first_failure_pattern);
                INFO("Text (" << first_failure_text.size() << " bp): "
                     << first_failure_text);
            }
            REQUIRE(failures == 0);
        }
    }
}

// =============================================================================
// TEST CASE 2: Edge cases
// =============================================================================

TEST_CASE("Validation - Edge cases", "[validation][edge]") {
    std::mt19937_64 rng(VALIDATION_SEED + 1);

    SECTION("pattern equals text length") {
        for (size_t len = 1; len <= 200; len += 7) {
            std::string pattern = random_sequence(rng, len);
            std::string text = random_sequence(rng, len);

            auto myers = myers_edit_distances(pattern, text);
            auto verifier = verify_edit_distances(pattern, text);

            INFO("Length: " << len);
            INFO("Pattern: " << pattern);
            INFO("Text: " << text);
            REQUIRE(myers.distances == verifier.distances);
        }
    }

    SECTION("pattern longer than text") {
        for (size_t i = 0; i < 50; ++i) {
            size_t pat_len = random_length(rng, 10, 100);
            size_t txt_len = random_length(rng, 1, pat_len - 1);

            std::string pattern = random_sequence(rng, pat_len);
            std::string text = random_sequence(rng, txt_len);

            auto myers = myers_edit_distances(pattern, text);
            auto verifier = verify_edit_distances(pattern, text);

            INFO("Pattern len: " << pat_len << ", Text len: " << txt_len);
            INFO("Pattern: " << pattern);
            INFO("Text: " << text);
            REQUIRE(myers.distances == verifier.distances);
        }
    }

    SECTION("all N pattern") {
        for (size_t len = 1; len <= 100; len += 11) {
            std::string pattern(len, 'N');
            std::string text = random_sequence(rng, len * 2);

            auto myers = myers_edit_distances(pattern, text);
            auto verifier = verify_edit_distances(pattern, text);

            INFO("Pattern length: " << len);
            INFO("Text: " << text);
            REQUIRE(myers.distances == verifier.distances);
        }
    }

    SECTION("all N text") {
        for (size_t len = 1; len <= 100; len += 11) {
            std::string pattern = random_sequence(rng, len);
            std::string text(len * 2, 'N');

            auto myers = myers_edit_distances(pattern, text);
            auto verifier = verify_edit_distances(pattern, text);

            INFO("Pattern: " << pattern);
            INFO("Text length: " << len * 2);
            REQUIRE(myers.distances == verifier.distances);
        }
    }

    SECTION("high N density (30%)") {
        for (size_t i = 0; i < 100; ++i) {
            size_t pat_len = random_length(rng, 10, 80);
            size_t txt_len = random_length(rng, 50, 200);

            std::string pattern = random_sequence(rng, pat_len, 0.30);
            std::string text = random_sequence(rng, txt_len, 0.30);

            auto myers = myers_edit_distances(pattern, text);
            auto verifier = verify_edit_distances(pattern, text);

            INFO("Pattern (" << pat_len << " bp): " << pattern);
            INFO("Text (" << txt_len << " bp): " << text);
            REQUIRE(myers.distances == verifier.distances);
        }
    }
}

// =============================================================================
// TEST CASE 3: Invariant properties
// =============================================================================

TEST_CASE("Validation - Invariant properties", "[validation][property]") {
    std::mt19937_64 rng(VALIDATION_SEED + 2);

    SECTION("distances never exceed pattern length") {
        for (size_t i = 0; i < 200; ++i) {
            size_t pat_len = random_length(rng, 1, 150);
            size_t txt_len = random_length(rng, 1, 500);

            std::string pattern = random_sequence(rng, pat_len);
            std::string text = random_sequence(rng, txt_len);

            auto myers = myers_edit_distances(pattern, text);

            for (size_t j = 0; j < myers.distances.size(); ++j) {
                INFO("Test index: " << i);
                INFO("Position: " << j);
                INFO("Distance: " << static_cast<int>(myers.distances[j]));
                INFO("Pattern length: " << pat_len);
                REQUIRE(myers.distances[j] <= pat_len);
            }
        }
    }

    SECTION("exact match embedded has distance 0 (unless pattern has N)") {
        for (size_t i = 0; i < 100; ++i) {
            size_t pat_len = random_length(rng, 5, 80);
            // Generate pattern without N to ensure true exact match
            std::string pattern = random_sequence(rng, pat_len, 0.0);

            // Embed pattern in random text (with no N density)
            size_t prefix_len = random_length(rng, 0, 50);
            size_t suffix_len = random_length(rng, 0, 50);

            std::string text = random_sequence(rng, prefix_len, 0.0) + pattern +
                               random_sequence(rng, suffix_len, 0.0);

            auto myers = myers_edit_distances(pattern, text);

            // Distance at end of embedded pattern should be 0
            size_t embed_end = prefix_len + pat_len - 1;
            INFO("Test index: " << i);
            INFO("Pattern (" << pat_len << " bp): " << pattern);
            INFO("Text (" << text.size() << " bp): " << text);
            INFO("Embed position: " << prefix_len << " to " << embed_end);
            INFO("Distance at embed_end: " << static_cast<int>(myers.distances[embed_end]));
            REQUIRE(myers.distances[embed_end] == 0);
        }
    }
}
