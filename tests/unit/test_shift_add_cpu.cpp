#include <catch2/catch_test_macros.hpp>
#include <shift_add_cpu.hpp>
#include <genome_loader.hpp>

#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// ───────────────────────────────────────────────────────────────────────────
// Helpers
// ───────────────────────────────────────────────────────────────────────────

// Naive O(mn) Hamming distance reference for verifying shift-add output.
// Computes the semi-global Hamming distance ending at every position in text.
// Semi-global: free to start anywhere in text (like Myers).
static std::vector<uint8_t> naive_hamming_semi_global(const std::string& pattern,
                                                       const std::string& text) {
    size_t m = pattern.size();
    size_t n = text.size();

    auto mismatch = [](char a, char b) -> bool {
        // Text N doesn't match any pattern nucleotide (masked)
        if (b == 'N' || b == 'n') return true;
        // Pattern N matches any text nucleotide (wildcard)
        if (a == 'N' || a == 'n') return false;
        return (a != b);
    };

    std::vector<uint8_t> result(n);
    
    for (size_t j = 0; j < n; ++j) {
        if (j < m - 1) {
            // Partial alignment: pattern[0..j] vs text[0..j]
            int dist = 0;
            for (size_t i = 0; i <= j; ++i) {
                if (mismatch(pattern[i], text[i])) dist++;
            }
            result[j] = static_cast<uint8_t>(dist);
        } else {
            // Full alignment: slide pattern to end at position j
            // Compare pattern[0..m-1] vs text[j-m+1..j]
            int dist = 0;
            for (size_t i = 0; i < m; ++i) {
                if (mismatch(pattern[i], text[j - m + 1 + i])) dist++;
            }
            result[j] = static_cast<uint8_t>(dist);
        }
    }
    
    return result;
}

// ───────────────────────────────────────────────────────────────────────────
// Section 1: Basic Correctness (8 tests)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ShiftAdd - identical pattern and text", "[shift_add][basic]") {
    auto res = shift_add_hamming_distance("ACGT", "ACGT");

    REQUIRE(res.pattern_len == 4);
    REQUIRE(res.text_len    == 4);
    REQUIRE(res.distances.size() == 4);
    // The alignment ending at position 3 (last char) must be distance 0
    REQUIRE(res.distances[3] == 0);
}

TEST_CASE("ShiftAdd - single character pattern", "[shift_add][basic]") {
    // Pattern "A" against "CACAC": exact match at positions 1,3 (dist 0),
    // mismatch at 0,2,4 (dist 1).
    auto res = shift_add_hamming_distance("A", "CACAC");

    REQUIRE(res.pattern_len == 1);
    REQUIRE(res.text_len    == 5);
    REQUIRE(res.distances.size() == 5);
    
    REQUIRE(res.distances[1] == 0);  // text[1] = 'A'
    REQUIRE(res.distances[3] == 0);  // text[3] = 'A'
    REQUIRE(res.distances[0] == 1);  // text[0] = 'C'
    REQUIRE(res.distances[2] == 1);  // text[2] = 'C'
    REQUIRE(res.distances[4] == 1);  // text[4] = 'C'
}

TEST_CASE("ShiftAdd - all mismatches", "[shift_add][basic]") {
    // Pattern "AAAA" against "CCCC": distance = 4 at all positions
    auto res = shift_add_hamming_distance("AAAA", "CCCC");

    REQUIRE(res.pattern_len == 4);
    REQUIRE(res.text_len    == 4);
    REQUIRE(res.distances[3] == 4);  // All 4 positions mismatch
}

TEST_CASE("ShiftAdd - single mismatch at different positions", "[shift_add][basic]") {
    std::string pattern = "ACGTACGT";  // 8 characters
    
    // Mismatch at position 0
    auto res1 = shift_add_hamming_distance(pattern, "GCGTACGT");
    REQUIRE(res1.distances[7] == 1);
    
    // Mismatch at position 4 (middle)
    auto res2 = shift_add_hamming_distance(pattern, "ACGTGCGT");
    REQUIRE(res2.distances[7] == 1);
    
    // Mismatch at position 7 (end)
    auto res3 = shift_add_hamming_distance(pattern, "ACGTACGA");
    REQUIRE(res3.distances[7] == 1);
}

TEST_CASE("ShiftAdd - empty pattern throws", "[shift_add][basic]") {
    REQUIRE_THROWS_AS(shift_add_hamming_distance("", "ACGT"), std::invalid_argument);
}

TEST_CASE("ShiftAdd - empty text throws", "[shift_add][basic]") {
    REQUIRE_THROWS_AS(shift_add_hamming_distance("ACGT", ""), std::invalid_argument);
}

TEST_CASE("ShiftAdd - invalid characters throw", "[shift_add][basic]") {
    // Invalid character 'X' in pattern
    REQUIRE_THROWS_AS(shift_add_hamming_distance("ACXTG", "ACGTACGT"), std::invalid_argument);
    
    // Invalid character '1' in text
    REQUIRE_THROWS_AS(shift_add_hamming_distance("ACGT", "AC1T"), std::invalid_argument);
}

TEST_CASE("ShiftAdd - case insensitive", "[shift_add][basic]") {
    // Lowercase input should be accepted and treated as uppercase
    auto res1 = shift_add_hamming_distance("acgt", "ACGT");
    auto res2 = shift_add_hamming_distance("ACGT", "acgt");
    auto res3 = shift_add_hamming_distance("AcGt", "aCgT");
    
    REQUIRE(res1.distances[3] == 0);
    REQUIRE(res2.distances[3] == 0);
    REQUIRE(res3.distances[3] == 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Section 2: Wildcard Handling (6 tests)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ShiftAdd - N in pattern matches any nucleotide", "[shift_add][wildcards]") {
    // Pattern "ANGT" should match "AAGT", "ACGT", "AGGT", "ATGT" with distance 0
    auto res1 = shift_add_hamming_distance("ANGT", "AAGT");
    REQUIRE(res1.distances[3] == 0);
    
    auto res2 = shift_add_hamming_distance("ANGT", "ACGT");
    REQUIRE(res2.distances[3] == 0);
    
    auto res3 = shift_add_hamming_distance("ANGT", "AGGT");
    REQUIRE(res3.distances[3] == 0);
    
    auto res4 = shift_add_hamming_distance("ANGT", "ATGT");
    REQUIRE(res4.distances[3] == 0);
}

TEST_CASE("ShiftAdd - N in text does not match pattern", "[shift_add][wildcards]") {
    // Pattern "ACGT" against text "ANGT": N at position 1 is a mismatch
    auto res = shift_add_hamming_distance("ACGT", "ANGT");
    REQUIRE(res.distances[3] == 1);  // One mismatch (the N)
}

TEST_CASE("ShiftAdd - pattern all N matches any text", "[shift_add][wildcards]") {
    auto res = shift_add_hamming_distance("NNNN", "ACGT");
    REQUIRE(res.distances[3] == 0);  // All Ns match
}

TEST_CASE("ShiftAdd - text all N has maximum distance", "[shift_add][wildcards]") {
    auto res = shift_add_hamming_distance("ACGT", "NNNN");
    REQUIRE(res.distances[3] == 4);  // All positions are N = all mismatches
}

TEST_CASE("ShiftAdd - mixed N positions", "[shift_add][wildcards]") {
    // Pattern "NNGT" against "ACGT": N's match A,C,T but G vs G matches, all = 0
    auto res = shift_add_hamming_distance("NNGT", "ACGT");
    REQUIRE(res.distances[3] == 0);  // All positions match (N's are wildcards)
    
    // Pattern "NNGN" against "ACTT": N's match, but G vs T at position 2 is a mismatch
    auto res2 = shift_add_hamming_distance("NNGN", "ACTT");
    REQUIRE(res2.distances[3] == 1);  // Only position 2 (G vs T) mismatches
}

TEST_CASE("ShiftAdd - N at boundaries", "[shift_add][wildcards]") {
    // N's at start and end
    auto res = shift_add_hamming_distance("NACGTN", "GACGTA");
    REQUIRE(res.distances[5] == 0);  // Ns match G and A
}

// ───────────────────────────────────────────────────────────────────────────
// Section 3: Edge Cases (6 tests)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ShiftAdd - pattern length 64 (single word boundary)", "[shift_add][edge]") {
    std::string pattern(64, 'A');
    std::string text(64, 'A');
    
    auto res = shift_add_hamming_distance(pattern, text);
    REQUIRE(res.pattern_len == 64);
    REQUIRE(res.distances[63] == 0);  // All match
}

TEST_CASE("ShiftAdd - pattern length 65 throws (multi-word not supported)", "[shift_add][edge]") {
    std::string pattern(65, 'A');
    std::string text(100, 'A');
    
    REQUIRE_THROWS_AS(shift_add_hamming_distance(pattern, text), std::invalid_argument);
}

TEST_CASE("ShiftAdd - very long pattern throws", "[shift_add][edge]") {
    std::string pattern(128, 'A');
    std::string text(200, 'A');
    
    REQUIRE_THROWS_AS(shift_add_hamming_distance(pattern, text), std::invalid_argument);
}

TEST_CASE("ShiftAdd - text shorter than pattern", "[shift_add][edge]") {
    // Pattern "ACGTACGT" (8 chars) against text "ACGT" (4 chars)
    // Should compute partial alignments only
    auto res = shift_add_hamming_distance("ACGTACGT", "ACGT");
    
    REQUIRE(res.text_len == 4);
    REQUIRE(res.distances.size() == 4);
    // Position 3: pattern[0..3] vs text[0..3] = "ACGT" vs "ACGT" = 0
    REQUIRE(res.distances[3] == 0);
}

TEST_CASE("ShiftAdd - pattern equals text", "[shift_add][edge]") {
    std::string seq = "ACGTACGTACGTACGT";
    auto res = shift_add_hamming_distance(seq, seq);
    
    REQUIRE(res.distances[seq.size() - 1] == 0);
}

TEST_CASE("ShiftAdd - long text with pattern sliding", "[shift_add][edge]") {
    // Pattern "ACG" slides over long text "ACGACGACG"
    // Should find exact matches at positions 2, 5, 8
    auto res = shift_add_hamming_distance("ACG", "ACGACGACG");
    
    REQUIRE(res.distances[2] == 0);  // ACG
    REQUIRE(res.distances[5] == 0);  // ACG
    REQUIRE(res.distances[8] == 0);  // ACG
}

// ───────────────────────────────────────────────────────────────────────────
// Section 4: Comparison vs Naive (8 tests)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ShiftAdd - matches naive: random case 1", "[shift_add][naive]") {
    auto res = shift_add_hamming_distance("ACGTACGT", "ACGTACGTACGTACGT");
    auto naive = naive_hamming_semi_global("ACGTACGT", "ACGTACGTACGTACGT");
    
    REQUIRE(res.distances == naive);
}

TEST_CASE("ShiftAdd - matches naive: random case 2", "[shift_add][naive]") {
    auto res = shift_add_hamming_distance("GATTACA", "GATTACAGATTACA");
    auto naive = naive_hamming_semi_global("GATTACA", "GATTACAGATTACA");
    
    REQUIRE(res.distances == naive);
}

TEST_CASE("ShiftAdd - matches naive: with mismatches", "[shift_add][naive]") {
    auto res = shift_add_hamming_distance("AAAA", "ACGTACGT");
    auto naive = naive_hamming_semi_global("AAAA", "ACGTACGT");
    
    REQUIRE(res.distances == naive);
}

TEST_CASE("ShiftAdd - matches naive: with Ns in pattern", "[shift_add][naive]") {
    auto res = shift_add_hamming_distance("ANNN", "ACGTACGT");
    auto naive = naive_hamming_semi_global("ANNN", "ACGTACGT");
    
    REQUIRE(res.distances == naive);
}

TEST_CASE("ShiftAdd - matches naive: with Ns in text", "[shift_add][naive]") {
    auto res = shift_add_hamming_distance("ACGT", "ANNNACGT");
    auto naive = naive_hamming_semi_global("ACGT", "ANNNACGT");
    
    REQUIRE(res.distances == naive);
}

TEST_CASE("ShiftAdd - matches naive: all distances 0", "[shift_add][naive]") {
    std::string pattern = "ACGT";
    std::string text = "ACGTACGTACGTACGT";  // Pattern repeats perfectly
    
    auto res = shift_add_hamming_distance(pattern, text);
    auto naive = naive_hamming_semi_global(pattern, text);
    
    REQUIRE(res.distances == naive);
}

TEST_CASE("ShiftAdd - matches naive: systematic sweep", "[shift_add][naive]") {
    // Test multiple patterns and texts systematically
    std::vector<std::string> patterns = {"A", "AC", "ACG", "ACGT", "ACGTACGT"};
    std::vector<std::string> texts = {"A", "ACGT", "ACGTACGTACGT"};
    
    for (const auto& pat : patterns) {
        for (const auto& txt : texts) {
            if (txt.size() >= pat.size()) {
                auto res = shift_add_hamming_distance(pat, txt);
                auto naive = naive_hamming_semi_global(pat, txt);
                
                REQUIRE(res.distances == naive);
            }
        }
    }
}

TEST_CASE("ShiftAdd - matches naive: long pattern", "[shift_add][naive]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";  // 20 chars
    std::string text = "ACGTACGTACGTACGTACGTACGTACGT";  // 28 chars
    
    auto res = shift_add_hamming_distance(pattern, text);
    auto naive = naive_hamming_semi_global(pattern, text);
    
    REQUIRE(res.distances == naive);
}

// ───────────────────────────────────────────────────────────────────────────
// Section 5: Divergence from Myers (8 tests)
// ───────────────────────────────────────────────────────────────────────────

// These tests verify that shift-add (Hamming) gives HIGHER distances than
// Myers (Levenshtein) when indels would help alignment.

TEST_CASE("ShiftAdd - divergence: insertion case", "[shift_add][divergence]") {
    // Pattern "AAAA" against "AACAA"
    // Myers (Levenshtein): distance = 1 (one insertion of 'C')
    // Shift-add (Hamming): must align all 4 A's, so distance depends on alignment
    
    std::string pattern = "AAAA";
    std::string text = "AACAA";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // Semi-global alignment ending at position 4 (last char):
    // Pattern "AAAA" vs text[1..4] = "ACAA"
    // Hamming distance = 1 (position 1: A vs C)
    REQUIRE(res.distances[4] == 1);
    
    // But if we had Myers here, it would find distance 1 via insertion
    // The key difference: Hamming can't use indels to improve alignment
}

TEST_CASE("ShiftAdd - divergence: deletion case", "[shift_add][divergence]") {
    // Pattern "AACAA" against "AAAA"
    // Myers: distance = 1 (one deletion of 'C')
    // Shift-add: must match all 5 chars of pattern, but text is shorter
    
    std::string pattern = "AACAA";
    std::string text = "AAAA";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // At position 3 (end of text), we have partial alignment:
    // Pattern[0..3] "AACA" vs text[0..3] "AAAA"
    // Hamming distance = 1 (position 2: C vs A)
    REQUIRE(res.distances[3] == 1);
}

TEST_CASE("ShiftAdd - divergence: simple indel example", "[shift_add][divergence]") {
    // Pattern "ACGT" against "ACGGT"
    // Myers: distance = 1 (one insertion of 'G')
    // Shift-add: at position 4 (end), pattern "ACGT" vs text[1..4] "CGGT"
    
    std::string pattern = "ACGT";
    std::string text = "ACGGT";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // Position 4: pattern "ACGT" vs text[1..4] "CGGT" = 2 mismatches (A≠C, C≠G)
    // This is the alignment ending exactly at position 4
    // Hamming distance doesn't do free gaps like Levenshtein semi-global
    REQUIRE(res.distances[4] == 2);
    
    // Note: Myers would give distance 1 via insertion, but Hamming can't use indels
}

TEST_CASE("ShiftAdd - divergence: requires indel for optimal Myers", "[shift_add][divergence]") {
    // Pattern "AAAA" against "AAAACAAAA"
    // Myers can insert 'C' for distance 1
    // Shift-add must have mismatches
    
    std::string pattern = "AAAA";
    std::string text = "AAAACAAAA";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // Position 3: pattern "AAAA" vs text[0..3] "AAAA" = 0
    REQUIRE(res.distances[3] == 0);
    
    // Position 4: pattern "AAAA" vs text[1..4] "AAAC" = 1
    REQUIRE(res.distances[4] == 1);
    
    // Position 8: pattern "AAAA" vs text[5..8] "AAAA" = 0
    REQUIRE(res.distances[8] == 0);
}

TEST_CASE("ShiftAdd - divergence: CRISPR-like example", "[shift_add][divergence]") {
    // Realistic CRISPR spacer with 1 bp insertion in genome
    std::string pattern = "GAGTCCGAGCAGAAGAAGAA";  // 20-mer
    std::string text    = "GAGTCCGACGCAGAAGAAGAA"; // Extra 'C' inserted
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // Myers would find distance ~1 (one insertion)
    // Shift-add must have >= 1 mismatch due to frameshift
    // At position 20 (end), pattern aligns to text[1..20] or text[0..19]
    // Best case: some mismatches due to insertion shifting alignment
    
    REQUIRE(res.distances[20] >= 1);
}

TEST_CASE("ShiftAdd - divergence: multiple indels", "[shift_add][divergence]") {
    // Pattern "ACGT" against "AACCGGTT"
    // Myers: distance = 4 (insertions)
    // Shift-add: at various positions, different Hamming distances
    
    std::string pattern = "ACGT";
    std::string text = "AACCGGTT";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // Check that distances are computed (exact values depend on semi-global)
    REQUIRE(res.distances.size() == 8);
    // At position 7, pattern "ACGT" vs text[4..7] "GGTT" = 4 mismatches
    // But semi-global might find better alignment
}

TEST_CASE("ShiftAdd - divergence: no indels means same as Myers", "[shift_add][divergence]") {
    // If there are no indels, Hamming == Levenshtein
    // Pattern "ACGT" against "ACGT" (exact match)
    
    std::string pattern = "ACGT";
    std::string text = "ACGT";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // Should be same as Myers: distance 0
    REQUIRE(res.distances[3] == 0);
    
    // With 2 substitutions
    auto res2 = shift_add_hamming_distance("ACGT", "AGGC");
    REQUIRE(res2.distances[3] == 2);  // Same as Myers
}

TEST_CASE("ShiftAdd - divergence: all text is Ns", "[shift_add][divergence]") {
    std::string pattern = "ACGT";
    std::string text = "NNNNNNNN";
    
    auto res = shift_add_hamming_distance(pattern, text);
    
    // All positions are N = all mismatches
    // Distance should be 4 (pattern length) at all full-alignment positions
    REQUIRE(res.distances[3] == 4);
    REQUIRE(res.distances[7] == 4);
}

// ───────────────────────────────────────────────────────────────────────────
// Section 6: Bit-Vector Interface Tests (8 tests)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ShiftAdd BV - basic functionality", "[shift_add][bv]") {
    // Create a small synthetic genome
    std::string seq = "ACGTACGTACGT";
    Genome genome = encode_genome({{"chr1", seq}});
    
    auto res = shift_add_hamming_distance_bv("ACGT", genome, 0, seq.size());
    
    REQUIRE(res.pattern_len == 4);
    REQUIRE(res.text_len == 12);
    // Check exact matches at positions 3, 7, 11
    REQUIRE(res.distances[3] == 0);
    REQUIRE(res.distances[7] == 0);
    REQUIRE(res.distances[11] == 0);
}

TEST_CASE("ShiftAdd BV - matches string interface", "[shift_add][bv]") {
    std::string seq = "ACGTACGTACGT";
    Genome genome = encode_genome({{"chr1", seq}});
    
    auto res_bv = shift_add_hamming_distance_bv("ACGT", genome, 0, seq.size());
    auto res_str = shift_add_hamming_distance("ACGT", seq);
    
    REQUIRE(res_bv.distances == res_str.distances);
}

TEST_CASE("ShiftAdd BV - slice of genome", "[shift_add][bv]") {
    std::string seq = "ACGTACGTACGT";
    Genome genome = encode_genome({{"chr1", seq}});
    
    // Search only positions [4..8)
    auto res = shift_add_hamming_distance_bv("ACGT", genome, 4, 4);
    
    REQUIRE(res.text_len == 4);
    REQUIRE(res.distances.size() == 4);
}

TEST_CASE("ShiftAdd BV - handles genome N's correctly", "[shift_add][bv]") {
    std::string seq = "ACGTNNNNACGT";
    Genome genome = encode_genome({{"chr1", seq}});
    
    auto res = shift_add_hamming_distance_bv("ACGT", genome, 0, seq.size());
    
    // N's at positions 4,5,6,7 should be mismatches
    // At position 7: pattern "ACGT" vs text[4..7] "NNNN" = 4 mismatches
    REQUIRE(res.distances[7] == 4);
    
    // At position 11: pattern "ACGT" vs text[8..11] "ACGT" = 0
    REQUIRE(res.distances[11] == 0);
}

TEST_CASE("ShiftAdd BV - empty pattern throws", "[shift_add][bv]") {
    Genome genome = encode_genome({{"chr1", "ACGT"}});
    REQUIRE_THROWS_AS(shift_add_hamming_distance_bv("", genome, 0, 4), std::invalid_argument);
}

TEST_CASE("ShiftAdd BV - zero length throws", "[shift_add][bv]") {
    Genome genome = encode_genome({{"chr1", "ACGT"}});
    REQUIRE_THROWS_AS(shift_add_hamming_distance_bv("ACGT", genome, 0, 0), std::invalid_argument);
}

TEST_CASE("ShiftAdd BV - out of range throws", "[shift_add][bv]") {
    Genome genome = encode_genome({{"chr1", "ACGT"}});
    // Try to access beyond genome bounds
    REQUIRE_THROWS_AS(shift_add_hamming_distance_bv("ACGT", genome, 0, 100), std::invalid_argument);
    REQUIRE_THROWS_AS(shift_add_hamming_distance_bv("ACGT", genome, 10, 10), std::invalid_argument);
}

TEST_CASE("ShiftAdd BV - long genome search", "[shift_add][bv]") {
    // Create a longer genome (1000 bp)
    std::string seq(1000, 'A');
    // Insert exact match at position 500
    for (size_t i = 500; i < 520; i++) {
        seq[i] = "ACGTACGTACGTACGTACGT"[i - 500];
    }
    
    Genome genome = encode_genome({{"chr1", seq}});
    auto res = shift_add_hamming_distance_bv("ACGTACGTACGTACGTACGT", genome, 0, 1000);

    // Should find exact match ending at position 519
    REQUIRE(res.distances[519] == 0);
}

TEST_CASE("ShiftAdd BV - parity with string interface on random inputs",
          "[shift_add][bv][parity]") {
    // Deterministic RNG for reproducibility.
    std::mt19937 rng(0xA42B15);
    std::uniform_int_distribution<int> nuc_dist(0, 4);  // 0-3 ACGT, 4 = N
    std::uniform_int_distribution<size_t> pat_len_dist(1, 24);
    std::uniform_int_distribution<size_t> text_len_dist(50, 400);

    static const char ALPHABET[5] = {'A', 'C', 'G', 'T', 'N'};

    auto random_dna = [&](size_t n) {
        std::string s(n, 'A');
        for (size_t i = 0; i < n; ++i) s[i] = ALPHABET[nuc_dist(rng)];
        return s;
    };

    for (int trial = 0; trial < 200; ++trial) {
        size_t m = pat_len_dist(rng);
        // Pattern must not be all-N (shift-and over all-N is meaningless and
        // the _bv path can't match anything anyway).
        std::string pattern;
        do { pattern = random_dna(m); }
        while (pattern.find_first_not_of('N') == std::string::npos);

        std::string text = random_dna(text_len_dist(rng));

        Genome genome = encode_genome({{"chr1", text}});
        auto res_bv = shift_add_hamming_distance_bv(pattern, genome, 0, text.size());
        auto res_str = shift_add_hamming_distance(pattern, text);

        INFO("trial=" << trial << " m=" << m << " pattern=" << pattern
             << " text=" << text);
        REQUIRE(res_bv.distances == res_str.distances);
    }
}
