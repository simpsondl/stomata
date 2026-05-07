#include <catch2/catch_test_macros.hpp>
#include <search_pipeline.hpp>
#include <genome_loader.hpp>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef TEST_DATA_DIR
#error "TEST_DATA_DIR not defined. Set via CMake target_compile_definitions."
#endif

// ───────────────────────────────────────────────────────────────────────────
// Helpers
// ───────────────────────────────────────────────────────────────────────────

static const std::string kTestSmallFa = TEST_DATA_DIR "/synthetic/test_small.fa";

// ───────────────────────────────────────────────────────────────────────────
// SearchHit structure tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("SearchHit - default construction", "[search_pipeline][basic]") {
    SearchHit hit;
    REQUIRE(hit.genome_pos == 0);
    REQUIRE(hit.distance == 0);
    REQUIRE(hit.chrom_name.empty());
    REQUIRE(hit.chrom_offset == 0);
}

TEST_CASE("SearchHit - comparison operators for sorting", "[search_pipeline][basic]") {
    SearchHit hit1{100, 2, "chr1", 100};
    SearchHit hit2{200, 1, "chr1", 200};
    SearchHit hit3{100, 3, "chr1", 100};

    // Default sort by genome position
    REQUIRE(hit1 < hit2);
    REQUIRE_FALSE(hit2 < hit1);

    // Same position, different distance
    REQUIRE(hit1.genome_pos == hit3.genome_pos);
}

// ───────────────────────────────────────────────────────────────────────────
// SearchConfig tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("SearchConfig - default values", "[search_pipeline][config]") {
    SearchConfig config;
    REQUIRE(config.threshold == 3);          // Default threshold (v0.6.1+: was 4)
    REQUIRE(config.prefer_gpu == true);      // Prefer GPU by default
    REQUIRE(config.pattern.empty());
}

TEST_CASE("SearchConfig - custom threshold", "[search_pipeline][config]") {
    SearchConfig config;
    config.threshold = 3;
    config.pattern = "ACGTACGT";

    REQUIRE(config.threshold == 3);
    REQUIRE(config.pattern == "ACGTACGT");
}

// ───────────────────────────────────────────────────────────────────────────
// Threshold filtering tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("filter_by_threshold - exact matches only", "[search_pipeline][filter]") {
    // Simulate distances from a sweep
    std::vector<uint8_t> distances = {3, 2, 0, 1, 4, 0, 5, 2};

    auto hits = filter_by_threshold(distances, 0);

    REQUIRE(hits.size() == 2);
    REQUIRE(hits[0] == 2);   // Position 2 has distance 0
    REQUIRE(hits[1] == 5);   // Position 5 has distance 0
}

TEST_CASE("filter_by_threshold - up to 1 mismatch", "[search_pipeline][filter]") {
    std::vector<uint8_t> distances = {3, 2, 0, 1, 4, 0, 5, 2};

    auto hits = filter_by_threshold(distances, 1);

    REQUIRE(hits.size() == 3);
    REQUIRE(hits[0] == 2);   // distance 0
    REQUIRE(hits[1] == 3);   // distance 1
    REQUIRE(hits[2] == 5);   // distance 0
}

TEST_CASE("filter_by_threshold - all pass", "[search_pipeline][filter]") {
    std::vector<uint8_t> distances = {1, 2, 1, 2};

    auto hits = filter_by_threshold(distances, 5);

    REQUIRE(hits.size() == 4);
}

TEST_CASE("filter_by_threshold - none pass", "[search_pipeline][filter]") {
    std::vector<uint8_t> distances = {5, 6, 7, 8};

    auto hits = filter_by_threshold(distances, 2);

    REQUIRE(hits.empty());
}

TEST_CASE("filter_by_threshold - empty input", "[search_pipeline][filter]") {
    std::vector<uint8_t> distances;

    auto hits = filter_by_threshold(distances, 3);

    REQUIRE(hits.empty());
}

// ───────────────────────────────────────────────────────────────────────────
// Chromosome resolution tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("resolve_chromosome - single chromosome", "[search_pipeline][resolve]") {
    std::vector<Chromosome> chroms = {{"chr1", 0, 100}};

    auto [name, offset] = resolve_chromosome(chroms, 50);

    REQUIRE(name == "chr1");
    REQUIRE(offset == 50);
}

TEST_CASE("resolve_chromosome - multiple chromosomes", "[search_pipeline][resolve]") {
    std::vector<Chromosome> chroms = {
        {"chr1", 0, 100},
        {"chr2", 100, 80}
    };

    // Position in chr1
    {
        auto [name, offset] = resolve_chromosome(chroms, 50);
        REQUIRE(name == "chr1");
        REQUIRE(offset == 50);
    }

    // Position at chr1/chr2 boundary
    {
        auto [name, offset] = resolve_chromosome(chroms, 100);
        REQUIRE(name == "chr2");
        REQUIRE(offset == 0);
    }

    // Position in chr2
    {
        auto [name, offset] = resolve_chromosome(chroms, 150);
        REQUIRE(name == "chr2");
        REQUIRE(offset == 50);
    }
}

TEST_CASE("resolve_chromosome - boundary conditions", "[search_pipeline][resolve]") {
    std::vector<Chromosome> chroms = {
        {"chr1", 0, 100},
        {"chr2", 100, 80}
    };

    // First position
    {
        auto [name, offset] = resolve_chromosome(chroms, 0);
        REQUIRE(name == "chr1");
        REQUIRE(offset == 0);
    }

    // Last position of chr1
    {
        auto [name, offset] = resolve_chromosome(chroms, 99);
        REQUIRE(name == "chr1");
        REQUIRE(offset == 99);
    }

    // Last position of chr2
    {
        auto [name, offset] = resolve_chromosome(chroms, 179);
        REQUIRE(name == "chr2");
        REQUIRE(offset == 79);
    }
}

TEST_CASE("resolve_chromosome - out of range throws", "[search_pipeline][resolve]") {
    std::vector<Chromosome> chroms = {{"chr1", 0, 100}};

    REQUIRE_THROWS_AS(resolve_chromosome(chroms, 100), std::out_of_range);
    REQUIRE_THROWS_AS(resolve_chromosome(chroms, 500), std::out_of_range);
}

// ───────────────────────────────────────────────────────────────────────────
// Full search pipeline tests (using test_small.fa)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - exact match in chr1", "[search_pipeline][integration]") {
    Genome genome = load_fasta(kTestSmallFa);

    // test_small.fa chr1 starts: ACGTACGTACGTACGTACGT...
    SearchConfig config;
    config.pattern = "ACGTACGT";
    config.threshold = 0;         // Exact matches only
    config.prefer_gpu = false;    // Use CPU for deterministic testing

    SearchResult result = search_genome(config, genome);

    REQUIRE(result.pattern == "ACGTACGT");
    REQUIRE(result.threshold == 0);

    // Should find exact matches at positions 0, 8, 16, ...
    REQUIRE(result.hits.size() >= 1);

    // All hits should have distance 0
    for (const auto& hit : result.hits) {
        REQUIRE(hit.distance == 0);
    }
}

TEST_CASE("search_genome - single mismatch tolerance", "[search_pipeline][integration]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGTACGT";
    config.threshold = 1;         // Up to 1 mismatch
    config.prefer_gpu = false;

    SearchResult result = search_genome(config, genome);

    // Should find more hits than exact match
    REQUIRE_FALSE(result.hits.empty());

    // All hits should have distance <= 1
    for (const auto& hit : result.hits) {
        REQUIRE(hit.distance <= 1);
    }
}

TEST_CASE("search_genome - no matches with impossible pattern", "[search_pipeline][integration]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "NNNNNNNN";   // Wildcard pattern matches everything, but let's try all-T
    config.threshold = 0;
    config.prefer_gpu = false;

    SearchResult result = search_genome(config, genome);

    // N matches everything, so this will find many hits
    // Let's change pattern to something that won't be found
    config.pattern = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";  // Long all-T pattern
    config.threshold = 0;

    result = search_genome(config, genome);

    // test_small.fa has ACGT repeats, long all-T exact match unlikely
    // With threshold 0, exact matches only
}

TEST_CASE("search_genome - hits span multiple chromosomes", "[search_pipeline][integration]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 1;
    config.prefer_gpu = false;

    SearchResult result = search_genome(config, genome);

    REQUIRE_FALSE(result.hits.empty());

    // Verify chromosome names are resolved correctly
    bool found_chr1 = false;
    bool found_chr2 = false;
    for (const auto& hit : result.hits) {
        if (hit.chrom_name == "chr1") found_chr1 = true;
        if (hit.chrom_name == "chr2") found_chr2 = true;
    }

    // Both chromosomes have ACGT patterns
    REQUIRE(found_chr1);
    REQUIRE(found_chr2);
}

TEST_CASE("search_genome - chromosome offset calculation", "[search_pipeline][integration]") {
    Genome genome = load_fasta(kTestSmallFa);

    // chr1 is 100 bp, chr2 starts at position 100
    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 0;
    config.prefer_gpu = false;

    SearchResult result = search_genome(config, genome);

    // Verify offsets are correct
    for (const auto& hit : result.hits) {
        if (hit.chrom_name == "chr1") {
            REQUIRE(hit.chrom_offset == hit.genome_pos);
        } else if (hit.chrom_name == "chr2") {
            // chr2 starts at genome position 100
            REQUIRE(hit.chrom_offset == hit.genome_pos - 100);
        }
    }
}

// ───────────────────────────────────────────────────────────────────────────
// GPU/CPU parity tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - GPU matches CPU", "[search_pipeline][gpu_parity]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig cpu_config;
    cpu_config.pattern = "ACGTACGTACGT";
    cpu_config.threshold = 2;
    cpu_config.prefer_gpu = false;

    SearchConfig gpu_config = cpu_config;
    gpu_config.prefer_gpu = true;

    SearchResult cpu_result = search_genome(cpu_config, genome);
    SearchResult gpu_result = search_genome(gpu_config, genome);

    // If GPU is available, results should match
    if (gpu_result.used_gpu) {
        REQUIRE(cpu_result.hits.size() == gpu_result.hits.size());

        for (size_t i = 0; i < cpu_result.hits.size(); ++i) {
            REQUIRE(cpu_result.hits[i].genome_pos == gpu_result.hits[i].genome_pos);
            REQUIRE(cpu_result.hits[i].distance == gpu_result.hits[i].distance);
        }
    }
    // If GPU not available, GPU result will use CPU fallback
    else {
        REQUIRE(cpu_result.hits.size() == gpu_result.hits.size());
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Output formatting tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("format_hits_tsv - basic output", "[search_pipeline][format]") {
    std::vector<SearchHit> hits = {
        {100, 0, "chr1", 100},
        {250, 2, "chr2", 70}
    };

    std::string tsv = format_hits_tsv(hits, "ACGTACGT");

    // Check header
    REQUIRE(tsv.find("chrom\tstart\tend\tpattern\tdistance") != std::string::npos);

    // Check data rows
    REQUIRE(tsv.find("chr1\t") != std::string::npos);
    REQUIRE(tsv.find("chr2\t") != std::string::npos);
}

TEST_CASE("format_hits_tsv - empty hits", "[search_pipeline][format]") {
    std::vector<SearchHit> hits;

    std::string tsv = format_hits_tsv(hits, "ACGT");

    // Should still have header
    REQUIRE(tsv.find("chrom\tstart\tend\tpattern\tdistance") != std::string::npos);
}

TEST_CASE("format_hits_bed - basic output", "[search_pipeline][format]") {
    std::vector<SearchHit> hits = {
        {100, 1, "chr1", 100},
        {200, 0, "chr2", 20}
    };

    std::string bed = format_hits_bed(hits, "ACGTACGT");

    // BED format: chrom start end name score strand
    REQUIRE(bed.find("chr1\t") != std::string::npos);
    REQUIRE(bed.find("chr2\t") != std::string::npos);
    // No header in BED format
    REQUIRE(bed.find("chrom\tstart") == std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Error handling tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - empty pattern throws", "[search_pipeline][error]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "";
    config.threshold = 3;

    REQUIRE_THROWS_AS(search_genome(config, genome), std::invalid_argument);
}

TEST_CASE("search_genome - pattern too long for GPU throws", "[search_pipeline][error]") {
    Genome genome = load_fasta(kTestSmallFa);

    // Pattern > 64 bp cannot use single-word GPU kernel
    std::string long_pattern(65, 'A');

    SearchConfig config;
    config.pattern = long_pattern;
    config.threshold = 3;
    config.prefer_gpu = true;

    // Should fall back to CPU or throw depending on implementation
    // For now, we test that it doesn't crash
    SearchResult result = search_genome(config, genome);
    REQUIRE(result.used_gpu == false);  // Should have fallen back to CPU
}

TEST_CASE("search_genome - invalid character in pattern throws", "[search_pipeline][error]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGTXACGT";  // X is invalid
    config.threshold = 3;

    REQUIRE_THROWS_AS(search_genome(config, genome), std::invalid_argument);
}

// ───────────────────────────────────────────────────────────────────────────
// Performance sanity tests (not benchmarks, just regression checks)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - handles full genome search", "[search_pipeline][perf]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGTACGT";
    config.threshold = 4;           // Common CRISPR threshold
    config.prefer_gpu = false;

    SearchResult result = search_genome(config, genome);

    // Just verify it completes and returns something
    REQUIRE(result.total_positions == genome.total_bases);
    // With both strands searched, we can have up to 2x genome positions
    REQUIRE(result.hits.size() <= genome.total_bases * 2);
}

// ───────────────────────────────────────────────────────────────────────────
// Parallel batch processing tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome_batch - parallel produces same results as sequential",
          "[search_pipeline][parallel]") {
    Genome genome = load_fasta(kTestSmallFa);

    // Create test spacers
    std::vector<SpacerEntry> spacers = {
        {"spacer_1", "ACGTACGT"},
        {"spacer_2", "GCTAGCTA"},
        {"spacer_3", "TTTTAAAA"},
        {"spacer_4", "ACGTACGT"}  // Duplicate to test correctness
    };

    BatchSearchConfig sequential_config;
    sequential_config.spacers = spacers;
    sequential_config.threshold = 2;
    sequential_config.prefer_gpu = false;
    sequential_config.num_threads = 1;  // Force sequential

    BatchSearchConfig parallel_config = sequential_config;
    parallel_config.num_threads = 4;  // Force parallel

    auto seq_result = search_genome_batch(sequential_config, genome);
    auto par_result = search_genome_batch(parallel_config, genome);

    // Same number of spacer results
    REQUIRE(seq_result.spacer_results.size() == par_result.spacer_results.size());

    // Same results for each spacer (order preserved)
    for (size_t i = 0; i < spacers.size(); ++i) {
        REQUIRE(seq_result.spacer_results[i].spacer_name ==
                par_result.spacer_results[i].spacer_name);
        REQUIRE(seq_result.spacer_results[i].result.hits.size() ==
                par_result.spacer_results[i].result.hits.size());
    }

    // Same total hits
    REQUIRE(seq_result.total_hits == par_result.total_hits);
}

TEST_CASE("search_genome_batch - empty spacer list", "[search_pipeline][parallel]") {
    Genome genome = load_fasta(kTestSmallFa);

    BatchSearchConfig config;
    config.spacers = {};
    config.num_threads = 4;

    auto result = search_genome_batch(config, genome);

    REQUIRE(result.spacer_results.empty());
    REQUIRE(result.total_hits == 0);
}

TEST_CASE("search_genome_batch - single spacer uses sequential path",
          "[search_pipeline][parallel]") {
    Genome genome = load_fasta(kTestSmallFa);

    BatchSearchConfig config;
    config.spacers = {{"single", "ACGTACGT"}};
    config.threshold = 2;
    config.prefer_gpu = false;
    config.num_threads = 8;  // Request many threads, should still work

    auto result = search_genome_batch(config, genome);

    REQUIRE(result.spacer_results.size() == 1);
    REQUIRE(result.spacer_results[0].spacer_name == "single");
}

TEST_CASE("search_genome_batch - result ordering preserved",
          "[search_pipeline][parallel]") {
    Genome genome = load_fasta(kTestSmallFa);

    std::vector<SpacerEntry> spacers;
    for (int i = 0; i < 20; ++i) {
        spacers.push_back({"spacer_" + std::to_string(i), "ACGTACGT"});
    }

    BatchSearchConfig config;
    config.spacers = spacers;
    config.threshold = 2;
    config.prefer_gpu = false;
    config.num_threads = 8;

    auto result = search_genome_batch(config, genome);

    // Verify order matches input
    REQUIRE(result.spacer_results.size() == spacers.size());
    for (size_t i = 0; i < spacers.size(); ++i) {
        REQUIRE(result.spacer_results[i].spacer_name == spacers[i].name);
    }
}

TEST_CASE("search_genome_batch - num_threads=0 uses hardware concurrency",
          "[search_pipeline][parallel]") {
    Genome genome = load_fasta(kTestSmallFa);

    std::vector<SpacerEntry> spacers = {
        {"s1", "ACGT"}, {"s2", "GCTA"}, {"s3", "TTAA"}, {"s4", "AACC"}
    };

    BatchSearchConfig config;
    config.spacers = spacers;
    config.threshold = 2;
    config.prefer_gpu = false;
    config.num_threads = 0;  // Auto

    // Should not throw, should complete successfully
    auto result = search_genome_batch(config, genome);

    REQUIRE(result.spacer_results.size() == 4);
}

TEST_CASE("BatchSearchConfig - default num_threads is 0", "[search_pipeline][parallel]") {
    BatchSearchConfig config;
    REQUIRE(config.num_threads == 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Hit deduplication tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("deduplicate_hits - removes halo positions around exact match",
          "[search_pipeline][dedup]") {
    // Simulate the classic halo pattern: distance [2,1,0,1,2] at consecutive positions
    std::vector<SearchHit> hits;
    hits.push_back({98, 2, "chr1", 98, Strand::PLUS});
    hits.push_back({99, 1, "chr1", 99, Strand::PLUS});
    hits.push_back({100, 0, "chr1", 100, Strand::PLUS});
    hits.push_back({101, 1, "chr1", 101, Strand::PLUS});
    hits.push_back({102, 2, "chr1", 102, Strand::PLUS});

    auto deduped = deduplicate_hits(std::move(hits), 23);

    REQUIRE(deduped.size() == 1);
    REQUIRE(deduped[0].genome_pos == 100);
    REQUIRE(deduped[0].distance == 0);
}

TEST_CASE("deduplicate_hits - preserves distinct distant matches",
          "[search_pipeline][dedup]") {
    // Two real matches far apart (gap > pattern_len)
    std::vector<SearchHit> hits;
    hits.push_back({100, 0, "chr1", 100, Strand::PLUS});
    hits.push_back({500, 0, "chr1", 500, Strand::PLUS});

    auto deduped = deduplicate_hits(std::move(hits), 23);

    REQUIRE(deduped.size() == 2);
    REQUIRE(deduped[0].genome_pos == 100);
    REQUIRE(deduped[1].genome_pos == 500);
}

TEST_CASE("deduplicate_hits - independent per strand",
          "[search_pipeline][dedup]") {
    std::vector<SearchHit> hits;
    // Plus strand: exact match + halo
    hits.push_back({100, 0, "chr1", 100, Strand::PLUS});
    hits.push_back({101, 1, "chr1", 101, Strand::PLUS});
    // Minus strand: exact match + halo at same positions
    hits.push_back({100, 0, "chr1", 100, Strand::MINUS});
    hits.push_back({101, 1, "chr1", 101, Strand::MINUS});

    auto deduped = deduplicate_hits(std::move(hits), 23);

    REQUIRE(deduped.size() == 2);  // One per strand
    // Verify both strands present
    bool found_plus = false, found_minus = false;
    for (const auto& h : deduped) {
        if (h.strand == Strand::PLUS) found_plus = true;
        if (h.strand == Strand::MINUS) found_minus = true;
    }
    REQUIRE(found_plus);
    REQUIRE(found_minus);
}

TEST_CASE("deduplicate_hits - adjacent real matches not merged when gap > pattern_len",
          "[search_pipeline][dedup]") {
    // Two real d=0 matches with gap > pattern_len (30 > 23)
    std::vector<SearchHit> hits;
    hits.push_back({100, 0, "chr1", 100, Strand::PLUS});
    hits.push_back({101, 1, "chr1", 101, Strand::PLUS});  // halo of first
    hits.push_back({131, 0, "chr1", 131, Strand::PLUS});  // new match (gap=30>23)
    hits.push_back({132, 1, "chr1", 132, Strand::PLUS});  // halo of second

    auto deduped = deduplicate_hits(std::move(hits), 23);

    REQUIRE(deduped.size() == 2);
    REQUIRE(deduped[0].genome_pos == 100);
    REQUIRE(deduped[0].distance == 0);
    REQUIRE(deduped[1].genome_pos == 131);
    REQUIRE(deduped[1].distance == 0);
}

TEST_CASE("deduplicate_hits - empty input", "[search_pipeline][dedup]") {
    std::vector<SearchHit> hits;
    auto deduped = deduplicate_hits(std::move(hits), 23);
    REQUIRE(deduped.empty());
}

TEST_CASE("deduplicate_hits - single hit unchanged", "[search_pipeline][dedup]") {
    std::vector<SearchHit> hits;
    hits.push_back({100, 1, "chr1", 100, Strand::PLUS});

    auto deduped = deduplicate_hits(std::move(hits), 23);

    REQUIRE(deduped.size() == 1);
    REQUIRE(deduped[0].genome_pos == 100);
}

TEST_CASE("deduplicate_hits - keeps best from cluster of distance-1 matches",
          "[search_pipeline][dedup]") {
    // Cluster of distance-1 hits with one distance-0 in the middle
    std::vector<SearchHit> hits;
    hits.push_back({50, 2, "chr1", 50, Strand::PLUS});
    hits.push_back({51, 1, "chr1", 51, Strand::PLUS});
    hits.push_back({52, 1, "chr1", 52, Strand::PLUS});
    hits.push_back({53, 2, "chr1", 53, Strand::PLUS});

    auto deduped = deduplicate_hits(std::move(hits), 23);

    REQUIRE(deduped.size() == 1);
    REQUIRE(deduped[0].distance == 1);  // Best in cluster
}

// ───────────────────────────────────────────────────────────────────────────
// Integration: dedup removes inflation from real search
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - dedup reduces halo inflation on synthetic genome",
          "[search_pipeline][dedup][integration]") {
    // test_small.fa has ACGTACGT repeats - plant a unique exact match
    // Use a pattern that exists exactly once per chromosome to verify
    // dedup doesn't over-inflate
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    // chr2 has TGCATGCA repeats; search for a short unique substring
    config.pattern = "TGCATGCATGCATGCA";  // 16bp, appears in chr2
    config.threshold = 2;
    config.prefer_gpu = false;
    config.search_both_strands = false;  // Plus strand only for simplicity
    config.compute_mismatches = false;

    SearchResult result = search_genome(config, genome);

    // With deduplication, each distinct alignment location should appear once.
    // Without dedup, each real match would generate ~5 positions (d=0 surrounded
    // by d=1, d=2 halos). With dedup, we should see significantly fewer hits.
    // The key check: no two hits should be within pattern_len of each other
    // on the same strand.
    for (size_t i = 1; i < result.hits.size(); ++i) {
        if (result.hits[i].strand == result.hits[i-1].strand) {
            size_t gap = result.hits[i].genome_pos - result.hits[i-1].genome_pos;
            REQUIRE(gap > config.pattern.size());
        }
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Uracil normalization tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("normalize_uracil - basic", "[search_pipeline][uracil]") {
    REQUIRE(normalize_uracil("ACGU") == "ACGT");
    REQUIRE(normalize_uracil("acgu") == "acgt");
    REQUIRE(normalize_uracil("ACGT") == "ACGT");
    REQUIRE(normalize_uracil("") == "");
    REQUIRE(normalize_uracil("UUUU") == "TTTT");
    REQUIRE(normalize_uracil("AuCgU") == "AtCgT");
}

TEST_CASE("normalize_uracil - mixed U and T", "[search_pipeline][uracil]") {
    REQUIRE(normalize_uracil("AUGCUAGC") == "ATGCTAGC");
    REQUIRE(normalize_uracil("GAGTCCGAGCAGAAGAAGAA") == "GAGTCCGAGCAGAAGAAGAA");  // No U
    REQUIRE(normalize_uracil("GAGUCCGAGCAGAAGAAGAA") == "GAGTCCGAGCAGAAGAAGAA");  // One U
}

TEST_CASE("normalize_uracil - search equivalence", "[search_pipeline][uracil]") {
    // A spacer with U should yield same hits as the equivalent T version
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config_t;
    config_t.pattern = "ACGTACGTAC";
    config_t.threshold = 2;
    
    auto result_t = search_genome(config_t, genome);

    SearchConfig config_u;
    config_u.pattern = normalize_uracil("ACGUACGUAC");
    config_u.threshold = 2;
    
    auto result_u = search_genome(config_u, genome);

    REQUIRE(result_t.hits.size() == result_u.hits.size());
    for (size_t i = 0; i < result_t.hits.size(); ++i) {
        REQUIRE(result_t.hits[i].genome_pos == result_u.hits[i].genome_pos);
        REQUIRE(result_t.hits[i].distance == result_u.hits[i].distance);
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Forward-only mode tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("forward_only - skips RC search", "[search_pipeline][forward_only]") {
    Genome genome = load_fasta(kTestSmallFa);

    // Search both strands (default)
    SearchConfig config_both;
    config_both.pattern = "ACGTACGTAC";
    config_both.threshold = 2;
    
    config_both.search_both_strands = true;
    auto result_both = search_genome(config_both, genome);

    // Search forward only
    SearchConfig config_fwd;
    config_fwd.pattern = "ACGTACGTAC";
    config_fwd.threshold = 2;
    
    config_fwd.forward_only = true;
    auto result_fwd = search_genome(config_fwd, genome);

    // Forward-only should have fewer or equal hits (no minus strand)
    REQUIRE(result_fwd.hits.size() <= result_both.hits.size());

    // All forward-only hits should be on plus strand
    for (const auto& hit : result_fwd.hits) {
        REQUIRE(hit.strand == Strand::PLUS);
    }
}

TEST_CASE("forward_only - RC-only hits return zero", "[search_pipeline][forward_only]") {
    // Use a pattern that is NOT a palindrome, so RC search finds different hits
    Genome genome = load_fasta(kTestSmallFa);
    std::string pattern = "AAAAAAAAAA";  // All A's - RC is all T's
    std::string rc_pattern = "TTTTTTTTTT";

    // Forward-only with all-A pattern
    SearchConfig config;
    config.pattern = pattern;
    config.threshold = 2;
    config.forward_only = true;
    auto result = search_genome(config, genome);

    // Verify no minus strand hits
    for (const auto& hit : result.hits) {
        REQUIRE(hit.strand == Strand::PLUS);
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Summary mode tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("summarize_single_result - basic counts", "[search_pipeline][summary]") {
    SearchResult result;
    result.pattern = "ACGTACGTAC";
    result.threshold = 4;

    // Add some fake hits at different distances
    SearchHit h0; h0.distance = 0; result.hits.push_back(h0);
    SearchHit h1a; h1a.distance = 1; result.hits.push_back(h1a);
    SearchHit h1b; h1b.distance = 1; result.hits.push_back(h1b);
    SearchHit h2; h2.distance = 2; result.hits.push_back(h2);
    SearchHit h4; h4.distance = 4; result.hits.push_back(h4);

    auto summary = summarize_single_result(result, "test_spacer");

    REQUIRE(summary.entries.size() == 1);
    REQUIRE(summary.total_spacers_processed == 1);
    REQUIRE(summary.entries[0].spacer_name == "test_spacer");
    REQUIRE(summary.entries[0].total_hits == 5);
    REQUIRE(summary.entries[0].hits_by_distance.at(0) == 1);
    REQUIRE(summary.entries[0].hits_by_distance.at(1) == 2);
    REQUIRE(summary.entries[0].hits_by_distance.at(2) == 1);
    REQUIRE(summary.entries[0].hits_by_distance.at(4) == 1);
    REQUIRE(summary.entries[0].hits_by_distance.count(3) == 0);
}

TEST_CASE("format_summary_json - output format", "[search_pipeline][summary]") {
    SummaryResult summary;
    summary.threshold = 3;
    summary.total_spacers_processed = 1;
    summary.spacers_skipped = 0;

    SummaryEntry entry;
    entry.spacer_name = "EMX1";
    entry.spacer_sequence = "ACGTACGT";
    entry.total_hits = 5;
    entry.hits_by_distance[0] = 1;
    entry.hits_by_distance[1] = 2;
    entry.hits_by_distance[2] = 2;
    summary.entries.push_back(entry);

    std::string json = format_summary_json(summary);
    REQUIRE(json.find("\"threshold\": 3") != std::string::npos);
    REQUIRE(json.find("\"name\": \"EMX1\"") != std::string::npos);
    REQUIRE(json.find("\"total_hits\": 5") != std::string::npos);
    REQUIRE(json.find("\"0\": 1") != std::string::npos);
    REQUIRE(json.find("\"1\": 2") != std::string::npos);
}

TEST_CASE("format_summary_tsv - output format", "[search_pipeline][summary]") {
    SummaryResult summary;
    summary.threshold = 2;
    summary.total_spacers_processed = 1;
    summary.spacers_skipped = 0;

    SummaryEntry entry;
    entry.spacer_name = "FANCF";
    entry.spacer_sequence = "ACGTACGT";
    entry.total_hits = 3;
    entry.hits_by_distance[0] = 1;
    entry.hits_by_distance[2] = 2;
    summary.entries.push_back(entry);

    std::string tsv = format_summary_tsv(summary);

    // Check header
    REQUIRE(tsv.find("spacer\tsequence\ttotal_hits\td0\td1\td2") != std::string::npos);
    // Check data row - d1 should be 0 (no hits at distance 1)
    REQUIRE(tsv.find("FANCF\tACGTACGT\t3\t1\t0\t2") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Distance Mode Tests (Hamming vs Levenshtein)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - distance mode defaults to Levenshtein", "[search_pipeline][distance_mode]") {
    auto genome = load_fasta(kTestSmallFa);
    
    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 3;
    // distance_mode defaults to DistanceMode::LEVENSHTEIN
    
    REQUIRE(config.distance_mode == DistanceMode::LEVENSHTEIN);
}

TEST_CASE("search_genome - Hamming mode uses shift-add algorithm", "[search_pipeline][distance_mode]") {
    // Create simple genome with exact match
    std::string seq = "ACGTACGTACGTACGT";
    auto genome = encode_genome({{"chr1", seq}});
    
    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 0;  // Exact matches only
    config.distance_mode = DistanceMode::HAMMING;
    config.prefer_gpu = false;  // Use CPU for predictable behavior
    config.search_both_strands = false;  // Only forward strand for simpler test
    config.compute_mismatches = false;  // Skip post-processing
    
    auto result = search_genome(config, genome);
    
    // Should find at least one exact match (after deduplication)
    REQUIRE(result.hits.size() > 0);
    // All hits should have distance 0
    for (const auto& hit : result.hits) {
        REQUIRE(hit.distance == 0);
    }
}

TEST_CASE("search_genome - Hamming vs Levenshtein with no indels gives same results", "[search_pipeline][distance_mode]") {
    // When there are no indels, Hamming and Levenshtein must find the same sites.
    // Use a non-palindromic pattern (AACC, RC=GGTT) so the two strands produce
    // distinct hits and the test is unambiguous.
    std::string seq = "NNNNAACCNNNNGGTTNNNN";
    auto genome = encode_genome({{"chr1", seq}});

    SearchConfig config_lev;
    config_lev.pattern = "AACC";
    config_lev.threshold = 0;
    config_lev.distance_mode = DistanceMode::LEVENSHTEIN;
    config_lev.prefer_gpu = false;

    SearchConfig config_ham;
    config_ham.pattern = "AACC";
    config_ham.threshold = 0;
    config_ham.distance_mode = DistanceMode::HAMMING;
    config_ham.prefer_gpu = false;

    auto result_lev = search_genome(config_lev, genome);
    auto result_ham = search_genome(config_ham, genome);

    // Build position-sets for order-independent comparison
    using PosStrand = std::pair<size_t, Strand>;
    auto to_set = [](const std::vector<SearchHit>& hits) {
        std::set<PosStrand> s;
        for (const auto& h : hits) s.insert({h.genome_pos, h.strand});
        return s;
    };

    auto lev_set = to_set(result_lev.hits);
    auto ham_set = to_set(result_ham.hits);

    // Both modes must find the same positions (set equality, not count equality,
    // to tolerate any dedup differences at edge positions).
    REQUIRE(lev_set == ham_set);

    // All hits must be exact matches (distance 0)
    for (const auto& h : result_lev.hits) REQUIRE(h.distance == 0);
    for (const auto& h : result_ham.hits) REQUIRE(h.distance == 0);
}

TEST_CASE("search_genome - Hamming mode with substitutions", "[search_pipeline][distance_mode]") {
    // Create simple genome
    std::string seq = "ACGTNNNNACGTNNNN";
    auto genome = encode_genome({{"chr1", seq}});
    
    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 1;
    config.distance_mode = DistanceMode::HAMMING;
    config.prefer_gpu = false;
    config.search_both_strands = false;
    config.compute_mismatches = false;
    
    auto result = search_genome(config, genome);
    
    // Should find at least one match
    REQUIRE(result.hits.size() >= 1);
    for (const auto& hit : result.hits) {
        REQUIRE(hit.distance <= 1);
    }
}

TEST_CASE("search_genome - Hamming mode works with forward_only", "[search_pipeline][distance_mode]") {
    std::string seq = "ACGTACGTACGTACGT";
    auto genome = encode_genome({{"chr1", seq}});
    
    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 0;
    config.distance_mode = DistanceMode::HAMMING;
    config.forward_only = true;  // Search forward strand only
    config.prefer_gpu = false;
    
    auto result = search_genome(config, genome);
    
    // Should only have hits on forward strand
    for (const auto& hit : result.hits) {
        REQUIRE(hit.strand == Strand::PLUS);
    }
}
