/**
 * Phase 3 — GPU Engine Tests
 *
 * Every test compares the GPU kernel output against the CPU reference
 * (myers_edit_distances / myers_edit_distances_bv) to verify bit-identical
 * correctness.  All tests skip when no CUDA device is present.
 *
 * Test organisation:
 *   1. Meta          — report GPU availability (never skips)
 *   2. Basic         — small hand-crafted inputs, single chunk
 *   3. Fixture       — parity on test_small.fa via the bv overload
 *   4. Multi-chunk   — large random inputs that span many chunks
 *   5. Boundary      — exact matches placed precisely at chunk boundaries
 *   6. Slice         — non-zero start offsets
 *   7. Wildcard      — N in pattern and/or text
 *   8. Error         — input validation
 */

#include <catch2/catch_test_macros.hpp>
#include <gpu_engine.hpp>
#include <myers_cpu.hpp>
#include <genome_loader.hpp>

#include <random>
#include <string>
#include <vector>

// ─── Skip helper ─────────────────────────────────────────────────────────────
#define GPU_SKIP_IF_UNAVAILABLE()                                \
    if (!gpu_available()) { SKIP("No CUDA device available"); }

// ─── Helpers ─────────────────────────────────────────────────────────────────

static Genome genome_from_string(const std::string& seq) {
    FastaEntry entry{"chr1", seq};
    return encode_genome({entry});
}

static std::string random_sequence(std::mt19937_64& rng, size_t len,
                                   double n_prob = 0.02) {
    static const char BASES[] = "ACGT";
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    std::uniform_int_distribution<int>     base_dist(0, 3);
    std::string seq;
    seq.reserve(len);
    for (size_t i = 0; i < len; ++i)
        seq += (prob(rng) < n_prob) ? 'N' : BASES[base_dist(rng)];
    return seq;
}

// =============================================================================
// TEST 1: Meta — GPU availability report (never skips)
// =============================================================================

TEST_CASE("GPU - availability report", "[gpu][meta]") {
    if (gpu_available()) {
        WARN("GPU is available — all GPU tests will run");
    } else {
        WARN("No GPU detected — GPU tests will be skipped");
    }
    REQUIRE(true);   // always passes; informational only
}

// =============================================================================
// TEST 2: Basic correctness — small inputs, single chunk
// =============================================================================

TEST_CASE("GPU Myers - basic correctness against CPU", "[gpu][basic]") {
    GPU_SKIP_IF_UNAVAILABLE();

    struct Case { const char* name; const char* pattern; const char* text; };
    const Case cases[] = {
        {"identical",       "ACGT",      "ACGT"},
        {"single sub",      "ACGT",      "ACCT"},
        {"single ins",      "ACG",       "ACGT"},
        {"single del",      "ACGT",      "ACG"},
        {"embedded exact",  "ACG",       "TTACGTTT"},
        {"all N pattern",   "NNNN",      "ACGTACGT"},
        {"mismatch heavy",  "AAAA",      "CCCCCCCC"},
        {"pat > text",      "ACGTACGT",  "ACG"},
        {"single char",     "A",         "ACGTACGT"},
        {"lowercase input", "acgt",      "ACGTacgt"},
    };

    for (const auto& c : cases) {
        SECTION(c.name) {
            std::string pat = c.pattern;
            std::string txt = c.text;
            Genome g = genome_from_string(txt);

            auto gpu = myers_gpu(pat, g, 0, txt.size());
            auto cpu = myers_edit_distances(pat, txt);

            INFO("pattern: " << pat);
            INFO("text:    " << txt);
            REQUIRE(gpu.pattern_len == cpu.pattern_len);
            REQUIRE(gpu.text_len   == cpu.text_len);
            REQUIRE(gpu.distances  == cpu.distances);
        }
    }
}

// =============================================================================
// TEST 3: Fixture parity — test_small.fa via the bv overload
// =============================================================================

TEST_CASE("GPU Myers - parity with CPU bv overload on test_small.fa",
          "[gpu][fixture]") {
    GPU_SKIP_IF_UNAVAILABLE();

    Genome genome = load_fasta(TEST_DATA_DIR "/synthetic/test_small.fa");

    const std::string patterns[] = {
        "ACGTACGTACGTACGTACGT",   // 20-mer, all four bases
        "NNNNACGTACGTACGTNNNN",   // N-flanked 20-mer
        "AAAAAAAAAAAAAAAAAAAAA",  // homo-A 21-mer
    };

    for (const auto& pat : patterns) {
        SECTION(pat) {
            for (size_t ci = 0; ci < genome.chromosomes.size(); ++ci) {
                const auto& chr = genome.chromosomes[ci];
                auto gpu = myers_gpu(pat, genome, chr.start, chr.length);
                auto cpu = myers_edit_distances_bv(pat, genome, chr.start, chr.length);

                INFO("chromosome: " << chr.name);
                INFO("pattern: " << pat);
                REQUIRE(gpu.distances == cpu.distances);
            }
        }
    }
}

// =============================================================================
// TEST 4: Multi-chunk random correctness
// =============================================================================

TEST_CASE("GPU Myers - multi-chunk random correctness",
          "[gpu][multichunk]") {
    GPU_SKIP_IF_UNAVAILABLE();

    static constexpr uint64_t SEED = 0xCA11AB1E00000001ULL;
    std::mt19937_64 rng(SEED);

    struct Spec { const char* name; size_t pat_len; size_t text_len; };
    const Spec specs[] = {
        {"200 k text  / 20 bp pat",  20,  200000},
        {"200 k text  / 50 bp pat",  50,  200000},
        {"200 k text  / 64 bp pat",  64,  200000},   // max single-word
        {"200 k text  /  1 bp pat",   1,  200000},
        { "70 k text  / 10 bp pat",  10,   70000},
    };

    for (const auto& spec : specs) {
        SECTION(spec.name) {
            std::string pat  = random_sequence(rng, spec.pat_len);
            std::string text = random_sequence(rng, spec.text_len);
            Genome g = genome_from_string(text);

            auto gpu = myers_gpu(pat, g, 0, text.size());
            auto cpu = myers_edit_distances(pat, text);

            INFO("pattern (" << spec.pat_len  << " bp): "
                 << pat.substr(0, 40) << (pat.size() > 40 ? "..." : ""));
            INFO("text length: " << spec.text_len);
            REQUIRE(gpu.distances == cpu.distances);
        }
    }
}

// =============================================================================
// TEST 5: Boundary — exact matches placed at chunk boundaries
// =============================================================================

TEST_CASE("GPU Myers - exact match crossing chunk boundary",
          "[gpu][boundary]") {
    GPU_SKIP_IF_UNAVAILABLE();

    constexpr size_t STRIDE = GPU_DEFAULT_STRIDE;
    std::string pattern = "ACGTACGTACGTACGTACGT";   // 20 bp
    size_t m = pattern.size();

    // Text slightly larger than one chunk; filled with T so the embedded
    // pattern stands out clearly.
    size_t text_len = STRIDE + 200;
    std::string text(text_len, 'T');

    // Single exact match that straddles the chunk boundary:
    // starts at STRIDE - 10, ends at STRIDE + 9.
    size_t embed_start = STRIDE - 10;
    for (size_t i = 0; i < m; ++i)
        text[embed_start + i] = pattern[i];
    size_t embed_end = embed_start + m - 1;   // STRIDE + 9

    Genome g = genome_from_string(text);
    auto gpu = myers_gpu(pattern, g, 0, text.size());
    auto cpu = myers_edit_distances(pattern, text);

    REQUIRE(gpu.distances == cpu.distances);
    REQUIRE(gpu.distances[embed_end] == 0);   // exact match must be found
}

TEST_CASE("GPU Myers - multiple matches across three chunks",
          "[gpu][boundary]") {
    GPU_SKIP_IF_UNAVAILABLE();

    constexpr size_t STRIDE = GPU_DEFAULT_STRIDE;
    std::string pattern = "ACGTACGT";   // 8 bp — short for easy spacing
    size_t m = pattern.size();

    size_t text_len = STRIDE * 3;
    std::string text(text_len, 'T');

    // Six non-overlapping embed positions spread across chunk boundaries.
    // Each range [pos, pos+m) is at least 100 positions away from the next.
    std::vector<size_t> positions = {
          1000,                        // middle of chunk 0
        STRIDE -   50 - m,            // well before boundary 0→1
        STRIDE -  m / 2,              // crosses boundary 0→1
        STRIDE + 1000,                // middle of chunk 1
        2 * STRIDE - m / 2,           // crosses boundary 1→2
        3 * STRIDE - m,               // last m positions of chunk 2
    };

    for (size_t pos : positions)
        for (size_t i = 0; i < m; ++i)
            text[pos + i] = pattern[i];

    Genome g = genome_from_string(text);
    auto gpu = myers_gpu(pattern, g, 0, text.size());
    auto cpu = myers_edit_distances(pattern, text);

    REQUIRE(gpu.distances == cpu.distances);
    for (size_t pos : positions) {
        size_t end = pos + m - 1;
        INFO("embedded match at [" << pos << ", " << end << "]");
        REQUIRE(gpu.distances[end] == 0);
    }
}

// =============================================================================
// TEST 6: Non-zero start slice
// =============================================================================

TEST_CASE("GPU Myers - non-zero start offset", "[gpu][slice]") {
    GPU_SKIP_IF_UNAVAILABLE();

    // Longer text so the slice is meaningful even with a large offset.
    static constexpr uint64_t SEED = 0xA11CE00000000002ULL;
    std::mt19937_64 rng(SEED);
    std::string full_text = random_sequence(rng, 5000);
    Genome g = genome_from_string(full_text);

    std::string pattern = "ACGTACGT";   // 8 bp

    // Several different slice windows
    struct SliceSpec { size_t start; size_t length; };
    const SliceSpec slices[] = {
        {   0,  100},
        { 100,  200},
        {2500, 1000},
        {4000,  999},
    };

    for (const auto& s : slices) {
        SECTION("slice [" + std::to_string(s.start) + ", "
                          + std::to_string(s.start + s.length) + ")") {
            auto gpu = myers_gpu(pattern, g, s.start, s.length);
            auto cpu = myers_edit_distances_bv(pattern, g, s.start, s.length);

            INFO("slice [" << s.start << ", " << s.start + s.length << ")");
            REQUIRE(gpu.pattern_len == pattern.size());
            REQUIRE(gpu.text_len   == s.length);
            REQUIRE(gpu.distances  == cpu.distances);
        }
    }
}

// =============================================================================
// TEST 7: N-wildcard handling
// =============================================================================

TEST_CASE("GPU Myers - N wildcard in pattern", "[gpu][wildcard]") {
    GPU_SKIP_IF_UNAVAILABLE();

    SECTION("N pattern matches ACGT nucleotides") {
        std::string pattern = "NNNN";
        std::string text    = "ACGTACGT";
        Genome g = genome_from_string(text);

        auto gpu = myers_gpu(pattern, g, 0, text.size());
        auto cpu = myers_edit_distances(pattern, text);

        REQUIRE(gpu.distances == cpu.distances);
        // All-N pattern of length 4: distance 0 at every position ≥ m-1
        for (size_t j = 3; j < text.size(); ++j)
            REQUIRE(gpu.distances[j] == 0);
    }

    SECTION("high N density in pattern") {
        static constexpr uint64_t SEED = 0xA11CE00000000003ULL;
        std::mt19937_64 rng(SEED);
        std::string pat  = random_sequence(rng, 30, 0.25);
        std::string text = random_sequence(rng, 10000, 0.0);  // no N in genome
        Genome g = genome_from_string(text);

        auto gpu = myers_gpu(pat, g, 0, text.size());
        auto cpu = myers_edit_distances(pat, text);

        REQUIRE(gpu.distances == cpu.distances);
    }
}

TEST_CASE("GPU Myers - N in genome is masked (does not match)", "[gpu][masked_n]") {
    GPU_SKIP_IF_UNAVAILABLE();

    SECTION("Single masked N in genome") {
        // Create genome "AANAA" where N is masked (0b0000)
        std::vector<FastaEntry> entries = {{"test", "AANAA"}};
        Genome g = encode_genome(entries);

        auto gpu = myers_gpu("A", g, 0, 5);
        auto cpu = myers_edit_distances_bv("A", g, 0, 5);

        REQUIRE(gpu.distances == cpu.distances);
        // Position 2 has masked N, should not match A
        REQUIRE(gpu.distances[2] == 1);
    }

    SECTION("Pattern AAA spanning masked N") {
        std::vector<FastaEntry> entries = {{"test", "AANAA"}};
        Genome g = encode_genome(entries);

        auto gpu = myers_gpu("AAA", g, 0, 5);
        auto cpu = myers_edit_distances_bv("AAA", g, 0, 5);

        REQUIRE(gpu.distances == cpu.distances);
        // N at position 2 causes mismatch
        REQUIRE(gpu.distances[2] == 1);
        REQUIRE(gpu.distances[3] == 1);
        REQUIRE(gpu.distances[4] == 1);
    }

    SECTION("Wildcard N pattern vs masked N genome") {
        std::vector<FastaEntry> entries = {{"test", "AANAA"}};
        Genome g = encode_genome(entries);

        auto gpu = myers_gpu("N", g, 0, 5);
        auto cpu = myers_edit_distances_bv("N", g, 0, 5);

        REQUIRE(gpu.distances == cpu.distances);
        // Wildcard N in pattern should match ACGT but not masked N
        REQUIRE(gpu.distances[0] == 0);  // matches A
        REQUIRE(gpu.distances[1] == 0);  // matches A
        REQUIRE(gpu.distances[2] == 1);  // masked N doesn't match wildcard
    }
}

// =============================================================================
// TEST 8: Error handling
// =============================================================================

TEST_CASE("GPU Myers - error handling", "[gpu][error]") {
    GPU_SKIP_IF_UNAVAILABLE();

    Genome g = genome_from_string("ACGTACGT");   // 8 bp

    SECTION("empty pattern") {
        REQUIRE_THROWS_AS(myers_gpu("", g, 0, 8), std::invalid_argument);
    }

    SECTION("pattern > 64 bp") {
        std::string long_pat(65, 'A');
        REQUIRE_THROWS_AS(myers_gpu(long_pat, g, 0, 8), std::invalid_argument);
    }

    SECTION("zero length") {
        REQUIRE_THROWS_AS(myers_gpu("ACGT", g, 0, 0), std::invalid_argument);
    }

    SECTION("slice out of range") {
        REQUIRE_THROWS_AS(myers_gpu("ACGT", g, 0, 100), std::invalid_argument);
    }

    SECTION("invalid character in pattern") {
        REQUIRE_THROWS_AS(myers_gpu("AXGT", g, 0, 8), std::invalid_argument);
    }
}
