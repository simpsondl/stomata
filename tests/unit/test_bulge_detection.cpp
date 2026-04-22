#include <catch2/catch_test_macros.hpp>
#include <genome_loader.hpp>
#include <search_pipeline.hpp>

#include <string>
#include <vector>

// ───────────────────────────────────────────────────────────────────────────
// Bulge Detection Tests
// ───────────────────────────────────────────────────────────────────────────
// These tests validate that the alignment traceback correctly identifies:
// - DNA bulges (insertions in genome = 'D' in CIGAR)
// - RNA bulges (insertions in pattern = 'I' in CIGAR)
// - Perfect matches (no indels)
// ───────────────────────────────────────────────────────────────────────────

namespace {

bool cigar_has_dna_bulge(const std::string& cigar) {
    return cigar.find('D') != std::string::npos;
}

bool cigar_has_rna_bulge(const std::string& cigar) {
    return cigar.find('I') != std::string::npos;
}

}  // namespace

TEST_CASE("Simple mismatch produces all-M CIGAR", "[bulge][mismatch]") {
    // Pattern: ACGTACGT (8 bp), genome has a single mismatch at the end.
    std::string pattern    = "ACGTACGT";
    std::string genome_seq = "ACGTACGANGG";

    Genome genome = encode_genome({{"chr1", genome_seq}});

    SearchConfig config;
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);
    REQUIRE(result.hits.size() >= 1);

    const SearchHit* mismatch_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && !cigar_has_indel(hit.mismatch_info.cigar)) {
            mismatch_hit = &hit;
            break;
        }
    }

    REQUIRE(mismatch_hit != nullptr);
    REQUIRE_FALSE(cigar_has_dna_bulge(mismatch_hit->mismatch_info.cigar));
    REQUIRE_FALSE(cigar_has_rna_bulge(mismatch_hit->mismatch_info.cigar));
}

TEST_CASE("DNA bulge appears as 'D' in CIGAR", "[bulge][dna_bulge]") {
    // Pattern: ACGTACGT (8 bp)
    // Genome has an extra T inserted after the first ACG -> DNA bulge (gap in pattern).
    std::string pattern    = "ACGTACGT";
    std::string genome_seq = "ACGTTACGTNGG";

    Genome genome = encode_genome({{"chr1", genome_seq}});

    SearchConfig config;
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);
    REQUIRE(result.hits.size() >= 1);

    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && cigar_has_dna_bulge(hit.mismatch_info.cigar)) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(cigar_has_dna_bulge(bulge_hit->mismatch_info.cigar));
    REQUIRE_FALSE(cigar_has_rna_bulge(bulge_hit->mismatch_info.cigar));
}

TEST_CASE("DNA bulge with 20bp pattern", "[bulge][dna_bulge][spacer]") {
    std::string pattern    = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTTACGTACGTACGTACGTNGG";

    Genome genome = encode_genome({{"chr1", genome_seq}});

    SearchConfig config;
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);
    REQUIRE(result.hits.size() >= 1);

    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && cigar_has_dna_bulge(hit.mismatch_info.cigar)) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(cigar_has_indel(bulge_hit->mismatch_info.cigar));
}

TEST_CASE("RNA bulge appears as 'I' in CIGAR", "[bulge][rna_bulge]") {
    // Pattern has an extra T not present in the genome -> RNA bulge (gap in genome).
    std::string pattern    = "ACGTTACGTACGTACGTACG";  // 20 bp
    std::string genome_seq = "ACGTACGTACGTACGTACGNGG";

    Genome genome = encode_genome({{"chr1", genome_seq}});

    SearchConfig config;
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);
    REQUIRE(result.hits.size() >= 1);

    const SearchHit* bulge_hit = nullptr;
    for (const auto& hit : result.hits) {
        if (hit.distance == 1 && cigar_has_rna_bulge(hit.mismatch_info.cigar)) {
            bulge_hit = &hit;
            break;
        }
    }

    REQUIRE(bulge_hit != nullptr);
    REQUIRE(cigar_has_rna_bulge(bulge_hit->mismatch_info.cigar));
    REQUIRE_FALSE(cigar_has_dna_bulge(bulge_hit->mismatch_info.cigar));
}

TEST_CASE("Perfect match has no indels", "[bulge][perfect_match]") {
    std::string pattern    = "ACGTACGTACGTACGTACGT";
    std::string genome_seq = "ACGTACGTACGTACGTACGTNGG";

    Genome genome = encode_genome({{"chr1", genome_seq}});

    SearchConfig config;
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
    REQUIRE_FALSE(cigar_has_indel(perfect_hit->mismatch_info.cigar));
    REQUIRE(perfect_hit->mismatch_info.aligned_sequence == pattern);
}

TEST_CASE("Bulge detection runs on minus strand", "[bulge][minus_strand]") {
    std::string pattern    = "ACGTACGTACGTACGTACGT";  // palindrome
    std::string genome_seq = "ACGTTACGTACGTACGTACGTNGG";

    Genome genome = encode_genome({{"chr1", genome_seq}});

    SearchConfig config;
    config.pattern = pattern;
    config.threshold = 1;
    config.search_both_strands = true;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);

    bool found_plus_bulge = false;
    for (const auto& hit : result.hits) {
        if (cigar_has_indel(hit.mismatch_info.cigar) && hit.strand == Strand::PLUS) {
            found_plus_bulge = true;
            break;
        }
    }
    REQUIRE(found_plus_bulge);
}
