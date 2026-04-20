#include <catch2/catch_test_macros.hpp>
#include <myers_cpu.hpp>
#include <genome_loader.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef TEST_DATA_DIR
#error "TEST_DATA_DIR not defined. Set via CMake target_compile_definitions."
#endif

// ───────────────────────────────────────────────────────────────────────────
// Helpers
// ───────────────────────────────────────────────────────────────────────────

// Naive O(mn) edit-distance reference for verifying Myers output.
// Computes the semi-global distance ending at every position in text.
static std::vector<uint8_t> naive_semi_global(const std::string& pattern,
                                              const std::string& text) {
    size_t m = pattern.size();
    size_t n = text.size();

    // dp[i][j] = edit distance between pattern[0..i-1] and text[j-k..j-1]
    // for the best starting position k.  Classic semi-global: free gaps at
    // the start of the text (first row = 0).
    std::vector<std::vector<uint8_t>> dp(m + 1, std::vector<uint8_t>(n + 1, 0));

    // First column: aligning pattern[0..i-1] against empty text = i deletions
    for (size_t i = 0; i <= m; ++i) dp[i][0] = static_cast<uint8_t>(i);
    // First row: free to start anywhere in text
    for (size_t j = 0; j <= n; ++j) dp[0][j] = 0;

    auto mismatch = [](char a, char b) -> bool {
        // Text N (b) is masked: does NOT match ANY pattern nucleotide
        if (b == 'N' || b == 'n') return true;
        // Pattern N (a) is wildcard: matches any text nucleotide A/C/G/T
        // (but text N was already handled above as mismatch)
        if (a == 'N' || a == 'n') return false;
        return (a != b);
    };

    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            uint8_t cost = mismatch(pattern[i-1], text[j-1]) ? 1 : 0;
            dp[i][j] = std::min({
                static_cast<uint8_t>(dp[i-1][j]   + 1),   // deletion
                static_cast<uint8_t>(dp[i][j-1]   + 1),   // insertion
                static_cast<uint8_t>(dp[i-1][j-1] + cost) // match/sub
            });
        }
    }

    // Last row = distances ending at each text position
    std::vector<uint8_t> result(n);
    for (size_t j = 0; j < n; ++j) {
        result[j] = dp[m][j + 1];
    }
    return result;
}

// ───────────────────────────────────────────────────────────────────────────
// Basic correctness
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers - identical pattern and text", "[myers_cpu][basic]") {
    auto res = myers_edit_distances("ACGT", "ACGT");

    REQUIRE(res.pattern_len == 4);
    REQUIRE(res.text_len    == 4);
    REQUIRE(res.distances.size() == 4);
    // The alignment ending at position 3 (last char) must be 0
    REQUIRE(res.distances[3] == 0);
}

TEST_CASE("Myers - single character pattern", "[myers_cpu][basic]") {
    // Pattern "A" against "CACAC": exact match at positions 1,3 (dist 0),
    // mismatch at 0,2,4 (dist 1).
    auto res = myers_edit_distances("A", "CACAC");

    REQUIRE(res.distances.size() == 5);
    REQUIRE(res.distances[0] == 1);   // C != A
    REQUIRE(res.distances[1] == 0);   // A == A
    REQUIRE(res.distances[2] == 1);   // C != A
    REQUIRE(res.distances[3] == 0);   // A == A
    REQUIRE(res.distances[4] == 1);   // C != A
}

TEST_CASE("Myers - single substitution", "[myers_cpu][basic]") {
    // One substitution in the middle
    auto res = myers_edit_distances("ACGT", "ACCT");

    REQUIRE(res.distances[3] == 1);   // one sub at position 2
}

TEST_CASE("Myers - single insertion in text", "[myers_cpu][basic]") {
    // Pattern ACGT, text ACGGT — one extra G inserted
    auto res = myers_edit_distances("ACGT", "ACGGT");

    REQUIRE(res.distances.size() == 5);
    // Best alignment ending at position 4 ('T') should be distance 1 (one insertion)
    REQUIRE(res.distances[4] == 1);
}

TEST_CASE("Myers - single deletion from text", "[myers_cpu][basic]") {
    // Pattern ACGT, text ACT — G deleted
    auto res = myers_edit_distances("ACGT", "ACT");

    REQUIRE(res.distances.size() == 3);
    // Best alignment ending at position 2 ('T') should be distance 1 (one deletion)
    REQUIRE(res.distances[2] == 1);
}

// ───────────────────────────────────────────────────────────────────────────
// N handling
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers - N handling: pattern N is wildcard, text N is masked", "[myers_cpu][n_handling]") {
    SECTION("N in pattern matches any text ACGT (but not text N)") {
        auto res = myers_edit_distances("N", "ACGT");
        // Pattern N matches every nucleotide: distance 0 at every position
        for (size_t i = 0; i < 4; ++i) {
            REQUIRE(res.distances[i] == 0);
        }
    }

    SECTION("N in text is masked (does NOT match pattern nucleotides)") {
        auto res = myers_edit_distances("A", "NANNN");
        // Position 0: N (masked) doesn't match A -> distance 1
        REQUIRE(res.distances[0] == 1);
        // Position 1: A matches A -> distance 0
        REQUIRE(res.distances[1] == 0);
        // Position 2: N doesn't match A -> distance 1
        REQUIRE(res.distances[2] == 1);
        // Position 3: N doesn't match A -> distance 1
        REQUIRE(res.distances[3] == 1);
        // Position 4: N doesn't match A -> distance 1
        REQUIRE(res.distances[4] == 1);
    }

    SECTION("Pattern N (wildcard) matches text N (masked) at distance 1") {
        auto res = myers_edit_distances("N", "N");
        // Pattern wildcard N doesn't match masked text N -> distance 1
        REQUIRE(res.distances[0] == 1);
    }
}

TEST_CASE("Myers - lowercase input accepted", "[myers_cpu][basic]") {
    // Should produce same result as uppercase
    auto upper = myers_edit_distances("ACGT", "ACGTACGT");
    auto lower = myers_edit_distances("acgt", "acgtacgt");
    REQUIRE(upper.distances == lower.distances);
}

// ───────────────────────────────────────────────────────────────────────────
// Exhaustive comparison against naive reference
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers - matches naive reference on short strings", "[myers_cpu][reference]") {
    struct Case {
        std::string pattern;
        std::string text;
    };

    std::vector<Case> cases = {
        {"A",      "A"},
        {"A",      "C"},
        {"AC",     "AC"},
        {"AC",     "CA"},
        {"ACG",    "ACGT"},
        {"ACGT",   "ACGT"},
        {"ACGT",   "ACGTACGT"},
        {"ACGT",   "TTTTACGTTTTT"},
        {"AAAA",   "AAAAAAAAAAAA"},
        {"ACGTACGT","ACGTACGT"},
        {"ACG",    "AACG"},
        {"ACG",    "ACGG"},
        {"ACGT",   "AGT"},          // deletion
        {"ACGT",   "AACGT"},        // insertion
        {"NNN",    "ACGT"},         // all-N pattern (wildcards)
    };

    for (const auto& c : cases) {
        SECTION(c.pattern + " vs " + c.text) {
            auto myers  = myers_edit_distances(c.pattern, c.text);
            auto naive  = naive_semi_global(c.pattern, c.text);
            REQUIRE(myers.distances == naive);
        }
    }
}

TEST_CASE("Myers - random-like 20-mer against 200-bp text", "[myers_cpu][reference]") {
    // A fixed pseudo-random pattern and text (deterministic, no RNG needed)
    std::string pattern = "ACGTACGTACGTACGTACGT";   // 20 bp
    std::string text    = "TTTTTTTTTTTTTTTTTTTT"
                          "ACGTACGTACGTACGTACGT"    // exact match at pos 20-39
                          "TTTTTTTTTTTTTTTTTTTT"
                          "ACGTACGTACGTACGTACGA"    // 1 sub at pos 60-79
                          "TTTTTTTTTTTTTTTTTTTT"
                          "ACGTACGTACGTACGTACG"     // 1 del  at pos 99-117 (19 bp)
                          "TTTTTTTTTTTTTTTTTTTTT"
                          "ACGTACGTACGTACGTACGTT"   // 1 ins  at pos 139-159 (21 bp)
                          "TTTTTTTTTTTTTTTTTTTT"
                          "NNNNNNNNNNNNNNNNNNNN";   // all-N at pos 180-199

    auto myers = myers_edit_distances(pattern, text);
    auto naive = naive_semi_global(pattern, text);
    REQUIRE(myers.distances == naive);
}

// ───────────────────────────────────────────────────────────────────────────
// Multi-word pattern (> 64 bp)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers - pattern longer than 64 bp (multi-word)", "[myers_cpu][multiword]") {
    // 80-bp pattern requires ceil(80/64) = 2 words
    std::string pattern(80, 'A');
    pattern[40] = 'C';   // one distinctive position

    // Text: 200 bp, contains an exact copy of pattern at position 60
    std::string text(200, 'T');
    for (size_t i = 0; i < 80; ++i) text[60 + i] = pattern[i];

    auto myers = myers_edit_distances(pattern, text);
    auto naive = naive_semi_global(pattern, text);
    REQUIRE(myers.distances == naive);
}

TEST_CASE("Myers - 128-bp pattern (exactly 2 words)", "[myers_cpu][multiword]") {
    std::string pattern(128, 'G');
    // Sprinkle some variety
    for (size_t i = 0; i < 128; i += 7) pattern[i] = 'A';

    std::string text(300, 'C');
    // Embed pattern at position 100
    for (size_t i = 0; i < 128; ++i) text[100 + i] = pattern[i];
    // Embed a 2-substitution copy at position 250 (fits: 250+128 > 300, so use 50)
    for (size_t i = 0; i < 128 && 50 + i < 300; ++i) text[50 + i] = pattern[i];
    text[55] = 'T';   // sub 1
    text[60] = 'T';   // sub 2

    auto myers = myers_edit_distances(pattern, text);
    auto naive = naive_semi_global(pattern, text);
    REQUIRE(myers.distances == naive);
}

// ───────────────────────────────────────────────────────────────────────────
// Edge cases
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers - pattern longer than text", "[myers_cpu][edge]") {
    // Pattern 6 bp, text 3 bp: minimum distance = 3 (3 deletions needed)
    auto res = myers_edit_distances("ACGTAC", "ACG");

    auto naive = naive_semi_global("ACGTAC", "ACG");
    REQUIRE(res.distances == naive);
    // The best possible ending at pos 2 is 3 deletions
    REQUIRE(res.distances[2] == 3);
}

TEST_CASE("Myers - pattern equals text length", "[myers_cpu][edge]") {
    auto res = myers_edit_distances("ACGT", "TGCA");
    auto naive = naive_semi_global("ACGT", "TGCA");
    REQUIRE(res.distances == naive);
}

TEST_CASE("Myers - single character text", "[myers_cpu][edge]") {
    auto res = myers_edit_distances("ACGT", "A");
    auto naive = naive_semi_global("ACGT", "A");
    REQUIRE(res.distances == naive);
    // Best alignment: match A, delete CGT = distance 3
    REQUIRE(res.distances[0] == 3);
}

// ───────────────────────────────────────────────────────────────────────────
// Error handling
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers - empty pattern throws", "[myers_cpu][error]") {
    REQUIRE_THROWS_AS(myers_edit_distances("", "ACGT"), std::invalid_argument);
}

TEST_CASE("Myers - empty text throws", "[myers_cpu][error]") {
    REQUIRE_THROWS_AS(myers_edit_distances("ACGT", ""), std::invalid_argument);
}

TEST_CASE("Myers - invalid character in pattern throws", "[myers_cpu][error]") {
    REQUIRE_THROWS_AS(myers_edit_distances("ACXGT", "ACGT"), std::invalid_argument);
}

TEST_CASE("Myers - invalid character in text throws", "[myers_cpu][error]") {
    REQUIRE_THROWS_AS(myers_edit_distances("ACGT", "ACXGT"), std::invalid_argument);
}

// ───────────────────────────────────────────────────────────────────────────
// Bit-vector overload — must match string overload on the same data
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Myers bv - matches string overload on synthetic FASTA", "[myers_cpu][bv]") {
    // Load the shared two-chromosome fixture
    static constexpr const char* SMALL_FASTA = TEST_DATA_DIR "/synthetic/test_small.fa";
    Genome genome = load_fasta(SMALL_FASTA);

    std::string pattern = "ACGTACGTACGTACGTACGT";   // 20 bp

    SECTION("sweep chr1") {
        // Reconstruct the plain text for chr1 from the bit-vectors
        const auto& ch = genome.chromosomes[0];
        std::string text;
        text.reserve(ch.length);
        for (size_t i = ch.start; i < ch.start + ch.length; ++i) {
            size_t word = i / 64;
            uint64_t bit = 1ULL << (i % 64);
            // Masked N is encoded as all four bits clear (0b0000).
            bool bA = genome.bv_A[word] & bit;
            bool bC = genome.bv_C[word] & bit;
            bool bG = genome.bv_G[word] & bit;
            bool bT = genome.bv_T[word] & bit;
            if      (!bA && !bC && !bG && !bT) text += 'N';
            else if (bA)                       text += 'A';
            else if (bC)                       text += 'C';
            else if (bG)                       text += 'G';
            else if (bT)                       text += 'T';
        }

        auto from_str = myers_edit_distances(pattern, text);
        auto from_bv  = myers_edit_distances_bv(pattern, genome, ch.start, ch.length);

        REQUIRE(from_bv.distances == from_str.distances);
    }

    SECTION("sweep chr2") {
        const auto& ch = genome.chromosomes[1];
        std::string text;
        text.reserve(ch.length);
        for (size_t i = ch.start; i < ch.start + ch.length; ++i) {
            size_t word = i / 64;
            uint64_t bit = 1ULL << (i % 64);
            // Masked N is encoded as all four bits clear (0b0000).
            bool bA = genome.bv_A[word] & bit;
            bool bC = genome.bv_C[word] & bit;
            bool bG = genome.bv_G[word] & bit;
            bool bT = genome.bv_T[word] & bit;
            if      (!bA && !bC && !bG && !bT) text += 'N';
            else if (bA)                       text += 'A';
            else if (bC)                       text += 'C';
            else if (bG)                       text += 'G';
            else if (bT)                       text += 'T';
        }

        auto from_str = myers_edit_distances(pattern, text);
        auto from_bv  = myers_edit_distances_bv(pattern, genome, ch.start, ch.length);

        REQUIRE(from_bv.distances == from_str.distances);
    }
}

TEST_CASE("Myers bv - out-of-range slice throws", "[myers_cpu][bv][error]") {
    Genome genome = load_fasta(TEST_DATA_DIR "/synthetic/test_small.fa");
    std::string pattern = "ACGT";

    REQUIRE_THROWS_AS(
        myers_edit_distances_bv(pattern, genome, 0, genome.total_bases + 1),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        myers_edit_distances_bv(pattern, genome, genome.total_bases, 1),
        std::invalid_argument);
}

TEST_CASE("Myers bv - empty pattern throws", "[myers_cpu][bv][error]") {
    Genome genome = load_fasta(TEST_DATA_DIR "/synthetic/test_small.fa");
    REQUIRE_THROWS_AS(
        myers_edit_distances_bv("", genome, 0, 10),
        std::invalid_argument);
}

TEST_CASE("Myers bv - N in genome is masked (does not match pattern)", "[myers_cpu][bv][n_handling]") {
    // Create a genome with N in the middle: "AANAA"
    std::vector<FastaEntry> entries = {{"test", "AANAA"}};
    Genome genome = encode_genome(entries);

    SECTION("Pattern A at N position has distance 1") {
        auto res = myers_edit_distances_bv("A", genome, 0, 5);
        // Position 0: A matches -> distance 0
        REQUIRE(res.distances[0] == 0);
        // Position 1: A matches -> distance 0
        REQUIRE(res.distances[1] == 0);
        // Position 2: N (masked) does not match A -> distance 1 (substitution)
        REQUIRE(res.distances[2] == 1);
        // Position 3: A matches -> semi-global can start fresh at position 3, distance 0
        REQUIRE(res.distances[3] == 0);
        // Position 4: A matches -> semi-global best alignment is AA (positions 3-4), distance 0
        REQUIRE(res.distances[4] == 0);
    }

    SECTION("Pattern AAA spanning N has distance 1") {
        auto res = myers_edit_distances_bv("AAA", genome, 0, 5);
        // Best alignment ending at position 2: AA[N] has 1 substitution
        REQUIRE(res.distances[2] == 1);
        // Best alignment ending at position 3: A[N]A has 1 sub, or delete 1st A and match [N]A (1 del), both = 1
        REQUIRE(res.distances[3] == 1);
        // Best alignment ending at position 4: [N]AA has 1 sub, or start fresh at 3 and match AA (1 del of first A)
        REQUIRE(res.distances[4] == 1);
    }

    SECTION("Pattern N (wildcard) matches ACGT but not masked genome N") {
        auto res = myers_edit_distances_bv("N", genome, 0, 5);
        // N in pattern is wildcard (matches ACGT but not masked N)
        // Masked N in genome (0b0000) does not match the wildcard pattern
        REQUIRE(res.distances[0] == 0);  // matches A
        REQUIRE(res.distances[1] == 0);  // matches A
        REQUIRE(res.distances[2] == 1);  // masked N doesn't match wildcard N
        REQUIRE(res.distances[3] == 0);  // semi-global: can start fresh at position 3, matches A
        REQUIRE(res.distances[4] == 0);  // matches A
    }
}
