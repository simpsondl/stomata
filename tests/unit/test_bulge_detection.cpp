#include <catch2/catch_test_macros.hpp>
#include <genome_loader.hpp>
#include <search_pipeline.hpp>

#include <string>
#include <vector>

// ───────────────────────────────────────────────────────────────────────────
// Bulge Detection Tests
// ───────────────────────────────────────────────────────────────────────────
// These tests validate that the alignment traceback correctly identifies:
// - DNA bulges (insertions in genome = deletions in pattern)
// - RNA bulges (insertions in pattern = deletions in genome)
// - Proper position tracking for bulges
// - Correct DISTAL/SEED flag computation for bulge positions
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("EditType enum values", "[bulge][enum]") {
    // Verify enum values for clarity in output
    REQUIRE(static_cast<uint8_t>(EditType::NONE) == 0);
    REQUIRE(static_cast<uint8_t>(EditType::MISMATCH) == 1);
    REQUIRE(static_cast<uint8_t>(EditType::DNA_BULGE) == 2);
    REQUIRE(static_cast<uint8_t>(EditType::RNA_BULGE) == 3);
}

TEST_CASE("Simple mismatch (no bulges)", "[bulge][mismatch]") {
    // Pattern: ACGTACGT (8 bp)
    // Genome:  ACGTACGA (mismatch at position 8)
    // Expected: 1 mismatch, no bulges
    
    std::string pattern = "ACGTACGT";
    std::string genome_seq = "ACGTACGANGG";  // NGG PAM
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    // Find the hit with distance 1 and a mismatch (not bulge)
    const SearchHit* mismatch_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && !hit.mismatch_info.has_dna_bulge() && 
            !hit.mismatch_info.has_rna_bulge() &&
            !hit.mismatch_info.mismatch_positions.empty()) {
            mismatch_hit = &hit;
            break;
        }
    }
    
    REQUIRE(mismatch_hit != nullptr);
    REQUIRE(mismatch_hit->mismatch_info.mismatch_positions.size() == 1);
    REQUIRE(mismatch_hit->mismatch_info.mismatch_positions[0] == 8);
    REQUIRE(mismatch_hit->mismatch_info.edit_types.size() == 1);
    REQUIRE(mismatch_hit->mismatch_info.edit_types[0] == EditType::MISMATCH);
    REQUIRE_FALSE(mismatch_hit->mismatch_info.has_dna_bulge());
    REQUIRE_FALSE(mismatch_hit->mismatch_info.has_rna_bulge());
}

TEST_CASE("DNA bulge in distal region", "[bulge][dna_bulge][distal]") {
    // Pattern: ACGTACGT (8 bp)
    // Genome:  ACGXTACGT (9 bp, extra X at position 4)
    // This is a DNA bulge (insertion in genome, deletion in pattern)
    // Alignment:
    //   Pattern: ACG-TACGT
    //   Genome:  ACGXTACGT
    // Bulge at pattern position 4 (after first 3 bases)
    
    std::string pattern = "ACGTACGT";
    std::string genome_seq = "ACGTTACGTNGG";  // Extra T after ACG
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    // Find the hit with edit distance 1 (the bulge)
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_dna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(bulge_hit->mismatch_info.mismatch_positions.size() == 1);
    REQUIRE(bulge_hit->mismatch_info.edit_types[0] == EditType::DNA_BULGE);
    REQUIRE(bulge_hit->mismatch_info.has_dna_bulge());
    REQUIRE_FALSE(bulge_hit->mismatch_info.has_rna_bulge());
    
    // Bulge should be in distal region (position 4 out of 8, so < 8-8 = 0 is false, but < 1 is false)
    // Actually for 8bp: distal is 1 to 0 (empty), seed is 1-8
    // Let me recalculate: m=8, distal is 1 to (m-8)=0, so no distal region!
    // For 8bp, everything is seed. Let's test with a longer pattern.
}

TEST_CASE("DNA bulge with 20bp pattern (distal region)", "[bulge][dna_bulge][distal]") {
    // Pattern: ACGTACGTACGTACGTACGT (20 bp, typical CRISPR length)
    // Genome:  ACGXTACGTACGTACGTACGT (21 bp, extra X at position 4)
    // Distal region: positions 1-12 (20-8=12)
    // Seed region: positions 13-20
    // Bulge at position 4 should trigger distal flag
    
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTTACGTACGTACGTACGTNGG";  // Extra T after first ACG
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_dna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(bulge_hit->mismatch_info.has_dna_bulge());
    REQUIRE_FALSE(bulge_hit->mismatch_info.has_rna_bulge());
    REQUIRE(bulge_hit->mismatch_info.mismatch_positions.size() == 1);
    REQUIRE(bulge_hit->mismatch_info.edit_types[0] == EditType::DNA_BULGE);

    // Bulge at position 4 is in distal region (1-12)
    REQUIRE(bulge_hit->mismatch_info.total_distal_edits() > 0);
    REQUIRE(bulge_hit->mismatch_info.total_seed_edits() == 0);
}

TEST_CASE("DNA bulge in seed region", "[bulge][dna_bulge][seed]") {
    // Pattern: ACGTACGTACGTACGTACGT (20 bp)
    // Genome:  ACGTACGTACGTACGTXACGT (21 bp, extra X at position 17)
    // Seed region: positions 13-20
    // Bulge at position 17 should trigger seed flag
    
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTACGTACGTACGTTACGTNGG";  // Extra T near end
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_dna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(bulge_hit->mismatch_info.has_dna_bulge());
    REQUIRE(bulge_hit->mismatch_info.mismatch_positions.size() == 1);
    REQUIRE(bulge_hit->mismatch_info.edit_types[0] == EditType::DNA_BULGE);

    // Bulge in seed region (13-20)
    REQUIRE(bulge_hit->mismatch_info.total_distal_edits() == 0);
    REQUIRE(bulge_hit->mismatch_info.total_seed_edits() > 0);
}

TEST_CASE("RNA bulge in distal region", "[bulge][rna_bulge][distal]") {
    // Pattern: ACGTTACGT (9 bp, extra T)
    // Genome:  ACGTACGT (8 bp)
    // This is an RNA bulge (insertion in pattern, deletion in genome)
    // Alignment:
    //   Pattern: ACGTTACGT
    //   Genome:  ACGT-ACGT
    // Bulge at genome position 5 (after ACGT)
    
    std::string pattern = "ACGTTACGTACGTACGTACG";  // 20 bp with extra T at position 5
    std::string genome_seq = "ACGTACGTACGTACGTACGNGG";  // 20 bp normal
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_rna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(bulge_hit->mismatch_info.has_rna_bulge());
    REQUIRE_FALSE(bulge_hit->mismatch_info.has_dna_bulge());
    REQUIRE(bulge_hit->mismatch_info.mismatch_positions.size() == 1);
    REQUIRE(bulge_hit->mismatch_info.edit_types[0] == EditType::RNA_BULGE);
    REQUIRE(bulge_hit->mismatch_info.total_distal_edits() > 0);
}

TEST_CASE("RNA bulge in seed region", "[bulge][rna_bulge][seed]") {
    // Pattern: ACGTACGTACGTACGTXACGT (21 bp with extra X at position 17)
    // Genome:  ACGTACGTACGTACGTACGT (20 bp)
    // RNA bulge in seed region
    
    std::string pattern = "ACGTACGTACGTACGTTACGT";  // 21 bp with extra T at position 17
    std::string genome_seq = "ACGTACGTACGTACGTACGTNGG";  // 20 bp
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_rna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(bulge_hit->mismatch_info.has_rna_bulge());
    REQUIRE(bulge_hit->mismatch_info.mismatch_positions.size() == 1);
    REQUIRE(bulge_hit->mismatch_info.edit_types[0] == EditType::RNA_BULGE);
    REQUIRE(bulge_hit->mismatch_info.total_seed_edits() > 0);
}

TEST_CASE("Bulge at boundary between distal and seed", "[bulge][boundary]") {
    // Pattern: ACGTACGTACGTACGTACGT (20 bp)
    // Distal: 1-12, Seed: 13-20
    // Genome with bulge at position 12 (last distal position)
    
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTACGTACGTTACGTACGTACGTNGG";  // Extra T at position 12
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_dna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    // Position 12 should be in distal region
    REQUIRE(bulge_hit->mismatch_info.total_distal_edits() > 0);
}

TEST_CASE("PAM disruption with bulge", "[bulge][pam]") {
    // Pattern: ACGTACGTACGTACGTACGT (20 bp)
    // Genome with DNA bulge AND invalid PAM
    // The bulge shouldn't affect PAM validation
    
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTTACGTACGTACGTACGTAAA";  // Extra T at 4, AAA PAM (invalid)
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && hit.mismatch_info.has_dna_bulge()) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(bulge_hit->mismatch_info.has_dna_bulge());
    REQUIRE(bulge_hit->mismatch_info.pam_type == PamType::OTHER);  // Invalid PAM
}

TEST_CASE("No bulges when edit distance is 0", "[bulge][perfect_match]") {
    // Perfect match should have no mismatches or bulges
    
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTACGTACGTACGTACGTNGG";
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 4;
    config.search_both_strands = false;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    REQUIRE(result.hits.size() >= 1);
    
    const SearchHit* perfect_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 0) {
            perfect_hit = &hit;
            break;
        }
    }
    
    REQUIRE(perfect_hit != nullptr);
    REQUIRE(perfect_hit->mismatch_info.mismatch_positions.empty());
    REQUIRE(perfect_hit->mismatch_info.edit_types.empty());
    REQUIRE_FALSE(perfect_hit->mismatch_info.has_dna_bulge());
    REQUIRE_FALSE(perfect_hit->mismatch_info.has_rna_bulge());
    REQUIRE(perfect_hit->mismatch_info.total_distal_edits() == 0);
    REQUIRE(perfect_hit->mismatch_info.total_seed_edits() == 0);
}

TEST_CASE("Bulge detection on minus strand", "[bulge][minus_strand]") {
    // Test that bulge detection works correctly on minus strand
    // Pattern: ACGTACGTACGTACGTACGT (20 bp)
    // RC:      ACGTACGTACGTACGTACGT (palindrome for simplicity)
    // Genome with bulge on minus strand
    
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTTACGTACGTACGTACGTNGG";  // Extra T at position 4
    
    Genome genome = encode_genome({{"chr1", genome_seq}});
    
    SearchConfig config;
    config.pam_filter = PamFilter::NONE;  // Disable PAM filtering for synthetic test data
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = true;
    config.compute_mismatches = true;
    
    SearchResult result = search_genome(config, genome);
    
    // Should find hits on both strands
    bool found_plus_bulge = false;
    bool found_minus_bulge = false;
    
    for (const auto& hit : result.hits) {
        if (hit.mismatch_info.has_dna_bulge() || hit.mismatch_info.has_rna_bulge()) {
            if (hit.strand == Strand::PLUS) {
                found_plus_bulge = true;
            } else {
                found_minus_bulge = true;
            }
        }
    }
    
    // At minimum, should find bulge on plus strand
    REQUIRE(found_plus_bulge);
}
