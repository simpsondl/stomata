// ───────────────────────────────────────────────────────────────────────────
// Unit tests for GPU shift-add Hamming distance implementation
//
// Test strategy:
//   - Validate GPU shift-add against CPU shift-add (ground truth)
//   - Test chunk boundary behavior (overlap warmup correctness)
//   - Test edge cases (exact matches, all mismatches, wildcards)
//   - Test multi-pattern batch modes (sparse and dense output)
// ───────────────────────────────────────────────────────────────────────────

#include <catch2/catch_test_macros.hpp>

#include <gpu_engine.hpp>
#include <shift_add_cpu.hpp>
#include <genome_loader.hpp>

#include <string>
#include <vector>
#include <algorithm>

// Helper: Create synthetic genome for testing
static Genome genome_from_string(const std::string& seq) {
    FastaEntry entry{"chr1", seq};
    return encode_genome({entry});
}

// Helper: Compare GPU vs CPU distance results
static void validate_gpu_vs_cpu(const std::vector<uint8_t>& gpu,
                                const std::vector<uint8_t>& cpu,
                                const std::string& test_name) {
    REQUIRE(gpu.size() == cpu.size());
    for (size_t i = 0; i < gpu.size(); ++i) {
        if (gpu[i] != cpu[i]) {
            FAIL("GPU vs CPU mismatch at position " << i << " in test '" << test_name << "': "
                 << "GPU=" << (int)gpu[i] << " CPU=" << (int)cpu[i]);
        }
    }
}

TEST_CASE("GPU shift-add availability", "[gpu][shift_add]") {
    if (!gpu_available()) {
        WARN("Skipping GPU shift-add tests: no CUDA device available");
        return;
    }
    
    REQUIRE(gpu_available());
}

TEST_CASE("GPU shift-add: exact match", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    std::string pattern = "ACGTACGTACGTACGT";  // 16bp pattern
    std::string genome_seq = "NNNNACGTACGTACGTACGTNNNN";  // Exact match in middle
    
    Genome genome = genome_from_string(genome_seq);
    
    // CPU reference
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    
    // GPU result
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    // Validate
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "exact_match");
    
    // Check expected distances
    REQUIRE(cpu_result.distances.size() == genome_seq.size());
    size_t match_pos = genome_seq.find(pattern);
    REQUIRE(match_pos != std::string::npos);
    REQUIRE(cpu_result.distances[match_pos + pattern.size() - 1] == 0);
}

TEST_CASE("GPU shift-add: single mismatch", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    std::string pattern = "ACGTACGT";
    std::string genome_seq = "NNNNACGTACGANNNN";  // Single mismatch (T->A at position 7)
    
    Genome genome = genome_from_string(genome_seq);
    
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "single_mismatch");
    
    // Check that minimum distance is 1
    uint8_t min_dist = *std::min_element(cpu_result.distances.begin(), cpu_result.distances.end());
    REQUIRE(min_dist == 1);
}

TEST_CASE("GPU shift-add: wildcard N in pattern", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    std::string pattern = "ACGTNACGT";  // N at position 4
    std::string genome_seq = "ACGTAACGT";  // A at position 4 (should match N at zero cost)
    
    Genome genome = genome_from_string(genome_seq);
    
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "wildcard_N");
    
    // N in pattern matches any nucleotide at zero cost
    REQUIRE(cpu_result.distances[pattern.size() - 1] == 0);
}

TEST_CASE("GPU shift-add: chunk boundary alignment", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    // Test pattern placed exactly at chunk boundaries (GPU_DEFAULT_STRIDE = 65536)
    std::string pattern = "ACGTACGTACGTACGT";  // 16bp
    
    // Create genome with exact match at stride boundary - 1
    size_t boundary = 65536;
    size_t match_end = boundary - 1;  // Match ends exactly at boundary
    size_t match_start = match_end - pattern.size() + 1;
    
    std::string genome_seq(boundary + 1000, 'A');
    genome_seq.replace(match_start, pattern.size(), pattern);
    
    Genome genome = genome_from_string(genome_seq);
    
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "chunk_boundary");
    
    // Verify exact match at expected position
    REQUIRE(cpu_result.distances[match_end] == 0);
}

TEST_CASE("GPU shift-add: all mismatches", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    std::string pattern = "AAAAAAAAAAAAAAAA";  // 16 A's
    std::string genome_seq = "TTTTTTTTTTTTTTTTTTTT";  // All T's (complete mismatch)
    
    Genome genome = genome_from_string(genome_seq);
    
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "all_mismatches");
    
    // All positions should have distance = pattern.size()
    for (size_t i = pattern.size() - 1; i < cpu_result.distances.size(); ++i) {
        REQUIRE(cpu_result.distances[i] == pattern.size());
    }
}

TEST_CASE("GPU shift-add: large genome (multi-chunk)", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    std::string pattern = "ACGTACGTACGTACGT";
    
    // Create 200KB genome (requires 3+ GPU chunks)
    size_t genome_len = 200000;
    std::string genome_seq(genome_len, 'N');
    
    // Place exact matches at 3 positions across different chunks
    size_t pos1 = 10000;
    size_t pos2 = 70000;
    size_t pos3 = 130000;
    genome_seq.replace(pos1, pattern.size(), pattern);
    genome_seq.replace(pos2, pattern.size(), pattern);
    genome_seq.replace(pos3, pattern.size(), pattern);
    
    Genome genome = genome_from_string(genome_seq);
    
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "large_genome");
    
    // Verify exact matches at all 3 positions
    REQUIRE(cpu_result.distances[pos1 + pattern.size() - 1] == 0);
    REQUIRE(cpu_result.distances[pos2 + pattern.size() - 1] == 0);
    REQUIRE(cpu_result.distances[pos3 + pattern.size() - 1] == 0);
}

TEST_CASE("GPU shift-add: maximum pattern length (64bp)", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    // Test single-word limit (64 nucleotides)
    std::string pattern = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    REQUIRE(pattern.size() == 64);
    
    std::string genome_seq = pattern + "NNNN" + pattern;
    
    Genome genome = genome_from_string(genome_seq);
    
    auto cpu_result = shift_add_hamming_distance_bv(pattern, make_view(genome), 0, genome.total_bases);
    auto gpu_result = shift_add_gpu(pattern, make_view(genome), 0, genome.total_bases);
    
    validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, "64bp_pattern");
    
    // Two exact matches
    REQUIRE(cpu_result.distances[63] == 0);  // First match
    REQUIRE(cpu_result.distances[63 + 4 + 64] == 0);  // Second match
}

TEST_CASE("GPU shift-add batch: small batch (< 32 patterns)", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    // Create 10 patterns
    std::vector<std::string> patterns;
    for (int i = 0; i < 10; ++i) {
        patterns.push_back("ACGTACGT" + std::string(i % 4, 'A'));  // Vary pattern
    }
    
    // Genome with first pattern embedded
    std::string genome_seq = "NNNN" + patterns[0] + "NNNN";
    Genome genome = genome_from_string(genome_seq);
    
    // Sparse GPU batch
    auto sparse_results = shift_add_gpu_batch_sparse(patterns, genome, 0, genome.total_bases, 2);
    REQUIRE(sparse_results.size() == patterns.size());
    
    // First pattern should have at least one hit (distance 0)
    REQUIRE(!sparse_results[0].hits.empty());
    bool found_exact = false;
    for (const auto& hit : sparse_results[0].hits) {
        if (hit.distance == 0) {
            found_exact = true;
            break;
        }
    }
    REQUIRE(found_exact);
    
    // Dense GPU batch (backward compatibility)
    auto dense_results = shift_add_gpu_batch_threshold(patterns, genome, 0, genome.total_bases, 2);
    REQUIRE(dense_results.size() == patterns.size());
    REQUIRE(dense_results[0].distances.size() == genome.total_bases);
}

TEST_CASE("GPU shift-add batch: large batch (> 32 patterns)", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    // Create 50 patterns (requires 2 kernel launches)
    std::vector<std::string> patterns;
    const std::string suffixes[] = {"", "A", "C", "G", "T", "AA", "AC", "AG", "AT", "CA"};
    for (int i = 0; i < 50; ++i) {
        patterns.push_back("ACGTACGT" + suffixes[i % 10]);
    }
    
    std::string genome_seq = "NNNN" + patterns[0] + "NNNN" + patterns[35] + "NNNN";
    Genome genome = genome_from_string(genome_seq);
    
    auto sparse_results = shift_add_gpu_batch_sparse(patterns, genome, 0, genome.total_bases, 2);
    REQUIRE(sparse_results.size() == 50);
    
    // Patterns 0 and 35 should have hits
    REQUIRE(!sparse_results[0].hits.empty());
    REQUIRE(!sparse_results[35].hits.empty());
}

TEST_CASE("GPU shift-add batch: threshold filtering", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    std::vector<std::string> patterns = {"AAAAAAAA", "AAAAAAAA"};
    
    // Genome with 0, 1, 2, 3 mismatches
    std::string genome_seq = "AAAAAAAANAAAAAAAAANNAAAAAAANNNAAAAAANNNNAAAAA";
    Genome genome = genome_from_string(genome_seq);
    
    // Threshold 0: only exact matches
    auto results_t0 = shift_add_gpu_batch_sparse(patterns, genome, 0, genome.total_bases, 0);
    int hits_t0 = 0;
    for (const auto& hit : results_t0[0].hits) {
        REQUIRE(hit.distance == 0);
        ++hits_t0;
    }
    REQUIRE(hits_t0 >= 1);  // At least one exact match
    
    // Threshold 2: distances 0, 1, 2
    auto results_t2 = shift_add_gpu_batch_sparse(patterns, genome, 0, genome.total_bases, 2);
    for (const auto& hit : results_t2[0].hits) {
        REQUIRE(hit.distance <= 2);
    }
    REQUIRE(results_t2[0].hits.size() >= results_t0[0].hits.size());
}

TEST_CASE("GPU shift-add vs CPU:RandomizerSeq test patterns", "[gpu][shift_add]") {
    if (!gpu_available()) return;

    // Real CRISPR spacers from test data
    std::vector<std::string> patterns = {
        "GTCATCTTAGTCATTACCTG",
        "GGGTGGAGTCTTCTAACAGGG",
        "ACAGGATCGGAAGAGCGTCG"
    };
    
    // Create synthetic genome with embedded patterns (with variations)
    std::string genome_seq = "NNNNNNNN";
    genome_seq += patterns[0];  // Exact match
    genome_seq += "NNNNNNNN";
    genome_seq += "GGGTGGAGTCTTCTAACAGGA";  // patterns[1] with 1 mismatch
    genome_seq += "NNNNNNNN";
    genome_seq += "ACAGGATCGGAAGAGCGCCG";  // patterns[2] with 2 mismatches
    genome_seq += "NNNNNNNN";
    
    Genome genome = genome_from_string(genome_seq);
    
    // Compare GPU vs CPU for all patterns
    for (size_t i = 0; i < patterns.size(); ++i) {
        auto cpu_result = shift_add_hamming_distance_bv(patterns[i], genome, 0, genome.total_bases);
        auto gpu_result = shift_add_gpu(patterns[i], genome, 0, genome.total_bases);
        
        validate_gpu_vs_cpu(gpu_result.distances, cpu_result.distances, 
                           "pattern_" + std::to_string(i));
    }
}
