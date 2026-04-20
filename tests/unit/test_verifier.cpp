#include <catch2/catch_test_macros.hpp>
#include <verifier.hpp>
#include <myers_cpu.hpp>

#include <stdexcept>
#include <string>
#include <vector>

// ───────────────────────────────────────────────────────────────────────────
// Full-array checks on small examples (every element hand-traced through DP)
// ───────────────────────────────────────────────────────────────────────────

// AC vs AC
//        ""  A   C
//   ""    0  0   0
//   A     1  0   1
//   C     2  1   0        → distances = [1, 0]
TEST_CASE("Verifier - full array: AC vs AC", "[verifier][full_array]") {
    auto res = verify_edit_distances("AC", "AC");

    REQUIRE(res.pattern_len      == 2);
    REQUIRE(res.text_len         == 2);
    REQUIRE(res.distances.size() == 2);
    REQUIRE(res.distances == std::vector<uint8_t>{1, 0});
}

// AC vs CA
//        ""  C   A
//   ""    0  0   0
//   A     1  1   0
//   C     2  1   1        → distances = [1, 1]
TEST_CASE("Verifier - full array: AC vs CA", "[verifier][full_array]") {
    auto res = verify_edit_distances("AC", "CA");
    REQUIRE(res.distances == std::vector<uint8_t>{1, 1});
}

// A vs CAC   (single-char pattern: 0 where text==A, 1 elsewhere)
TEST_CASE("Verifier - full array: A vs CAC", "[verifier][full_array]") {
    auto res = verify_edit_distances("A", "CAC");
    REQUIRE(res.distances == std::vector<uint8_t>{1, 0, 1});
}

// ACT vs ACGT   (one insertion in text: extra G)
//        ""  A   C   G   T
//   ""    0  0   0   0   0
//   A     1  0   1   1   1
//   C     2  1   0   1   2
//   T     3  2   1   1   1        → distances = [2, 1, 1, 1]
TEST_CASE("Verifier - full array: ACT vs ACGT (insertion in text)", "[verifier][full_array]") {
    auto res = verify_edit_distances("ACT", "ACGT");
    REQUIRE(res.distances == std::vector<uint8_t>{2, 1, 1, 1});
}

// ACGT vs ACT   (one deletion from text: G missing)
//        ""  A   C   T
//   ""    0  0   0   0
//   A     1  0   1   1
//   C     2  1   0   1
//   G     3  2   1   1
//   T     4  3   2   1        → distances = [3, 2, 1]
TEST_CASE("Verifier - full array: ACGT vs ACT (deletion from text)", "[verifier][full_array]") {
    auto res = verify_edit_distances("ACGT", "ACT");
    REQUIRE(res.distances == std::vector<uint8_t>{3, 2, 1});
}

// ───────────────────────────────────────────────────────────────────────────
// Basic correctness (spot-check key positions)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Verifier - identical strings", "[verifier][basic]") {
    auto res = verify_edit_distances("ACGT", "ACGT");
    REQUIRE(res.distances[3] == 0);   // exact match
}

TEST_CASE("Verifier - single substitution", "[verifier][basic]") {
    auto res = verify_edit_distances("ACGT", "ACCT");
    REQUIRE(res.distances[3] == 1);   // G→C at position 2
}

TEST_CASE("Verifier - pattern embedded in longer text", "[verifier][basic]") {
    // ACGT sits at text[2..5]; exact match ends at position 5
    auto res = verify_edit_distances("ACGT", "TTACGTTT");
    REQUIRE(res.distances[5] == 0);
    REQUIRE(res.distances[0] > 0);    // no full match possible this early
}

TEST_CASE("Verifier - lowercase accepted", "[verifier][basic]") {
    auto upper = verify_edit_distances("ACGT", "ACGTACGT");
    auto lower = verify_edit_distances("acgt", "acgtacgt");
    REQUIRE(upper.distances == lower.distances);
}

// ───────────────────────────────────────────────────────────────────────────
// N handling (wildcard semantics)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Verifier - N in pattern matches every nucleotide", "[verifier][n_handling]") {
    auto res = verify_edit_distances("N", "ACGT");
    // Single-char wildcard pattern: distance 0 at every position
    for (size_t i = 0; i < 4; ++i) {
        REQUIRE(res.distances[i] == 0);
    }
}

TEST_CASE("Verifier - N in text is masked (doesn't match pattern)", "[verifier][n_handling]") {
    auto res = verify_edit_distances("ACGT", "NNNN");
    // Each text N is masked and doesn't match pattern nucleotides
    // Distance = 4 (all substitutions)
    REQUIRE(res.distances[3] == 4);
}

TEST_CASE("Verifier - pattern N (wildcard) doesn't match text N (masked)", "[verifier][n_handling]") {
    auto res = verify_edit_distances("N", "N");
    // Pattern wildcard doesn't match masked text N -> distance 1
    REQUIRE(res.distances[0] == 1);
}

TEST_CASE("Verifier - single N embedded in pattern", "[verifier][n_handling]") {
    // ACNT: the N at position 2 is wildcard and matches G at zero cost
    auto res = verify_edit_distances("ACNT", "ACGT");
    REQUIRE(res.distances[3] == 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Edge cases
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Verifier - pattern longer than text", "[verifier][edge]") {
    // 6-bp pattern vs 3-bp text: best = match ACG, delete TAC = 3 dels
    auto res = verify_edit_distances("ACGTAC", "ACG");
    REQUIRE(res.distances[2] == 3);
}

TEST_CASE("Verifier - single character text", "[verifier][edge]") {
    // 4-bp pattern vs 1-bp text: match A, delete CGT = 3 dels
    auto res = verify_edit_distances("ACGT", "A");
    REQUIRE(res.distances[0] == 3);
}

TEST_CASE("Verifier - all mismatches, equal length", "[verifier][edge]") {
    // No nucleotide in common; every position requires a substitution
    auto res = verify_edit_distances("AAAA", "CCCC");
    REQUIRE(res.distances[3] == 4);
}

// ───────────────────────────────────────────────────────────────────────────
// Structural properties (must hold for any valid input)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Verifier - distance never exceeds pattern length", "[verifier][property]") {
    // Semi-global invariant: d[j] ≤ m  (worst case = delete entire pattern)
    std::string pattern = "ACGTACGT";                   // 8 bp
    std::string text    = "TTTTTTTTTTTTTTTTTTTT";       // 20 T's

    auto res = verify_edit_distances(pattern, text);
    REQUIRE(res.distances.size() == 20);
    for (size_t j = 0; j < 20; ++j) {
        REQUIRE(res.distances[j] <= static_cast<uint8_t>(res.pattern_len));
    }
}

TEST_CASE("Verifier - exact-match positions score zero", "[verifier][property]") {
    // Two non-overlapping exact copies at known positions
    std::string text = "ACGTTTTTACGTTTTT";   // copies end at [3] and [11]
    auto res = verify_edit_distances("ACGT", text);

    REQUIRE(res.distances[3]  == 0);
    REQUIRE(res.distances[11] == 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Cross-validation with Myers (the primary purpose of the verifier)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Verifier - matches Myers on diverse inputs", "[verifier][cross_validation]") {
    struct Case { std::string pattern; std::string text; };

    std::vector<Case> cases = {
        {"A",         "A"},
        {"A",         "CACAC"},
        {"AC",        "CA"},
        {"ACG",       "ACGT"},
        {"ACGT",      "ACGT"},
        {"ACGT",      "ACGTACGT"},
        {"ACGT",      "TTTTACGTTTTT"},
        {"ACGT",      "ACCT"},                          // 1 sub
        {"ACT",       "ACGT"},                          // 1 ins in text
        {"ACGT",      "ACT"},                           // 1 del from text
        {"NNN",       "ACGT"},                          // N in pattern (wildcard)
        {"ACGTAC",    "ACG"},                           // pattern > text
        {"ACGT",      "A"},                             // single-char text
        {"AAAA",      "CCCC"},                          // all mismatches
        {"ACGTACGT",  "TTTTTTTTTTTTTTTTTTTT"},          // no match region
        {"ACGT",      "TTACGTTTACCTTTACTTTACGTTTT"},    // multiple near-matches
    };

    for (const auto& c : cases) {
        SECTION(c.pattern + " vs " + c.text) {
            auto ver   = verify_edit_distances(c.pattern, c.text);
            auto myers = myers_edit_distances(c.pattern, c.text);
            REQUIRE(ver.distances == myers.distances);
        }
    }
}

TEST_CASE("Verifier - matches Myers on 80-bp multi-word pattern", "[verifier][cross_validation]") {
    // 80 bp forces Myers into its multi-word path; verifier is word-agnostic
    std::string pattern(80, 'A');
    pattern[40] = 'C';

    std::string text(200, 'T');
    for (size_t i = 0; i < 80; ++i) text[60 + i] = pattern[i];   // exact copy at [60..139]

    auto ver   = verify_edit_distances(pattern, text);
    auto myers = myers_edit_distances(pattern, text);
    REQUIRE(ver.distances == myers.distances);
}

// ───────────────────────────────────────────────────────────────────────────
// Error handling
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Verifier - empty pattern throws", "[verifier][error]") {
    REQUIRE_THROWS_AS(verify_edit_distances("", "ACGT"), std::invalid_argument);
}

TEST_CASE("Verifier - empty text throws", "[verifier][error]") {
    REQUIRE_THROWS_AS(verify_edit_distances("ACGT", ""), std::invalid_argument);
}

TEST_CASE("Verifier - invalid character in pattern throws", "[verifier][error]") {
    REQUIRE_THROWS_AS(verify_edit_distances("ACXGT", "ACGT"), std::invalid_argument);
}

TEST_CASE("Verifier - invalid character in text throws", "[verifier][error]") {
    REQUIRE_THROWS_AS(verify_edit_distances("ACGT", "ACXGT"), std::invalid_argument);
}
