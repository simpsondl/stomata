#include <catch2/catch_test_macros.hpp>
#include <genome_loader.hpp>
#include <search_pipeline.hpp>

#include <stdexcept>
#include <string>
#include <vector>

#ifndef TEST_DATA_DIR
#error "TEST_DATA_DIR not defined. Set via CMake target_compile_definitions."
#endif

static const std::string kTestSmallFa = TEST_DATA_DIR "/synthetic/test_small.fa";

// ───────────────────────────────────────────────────────────────────────────
// reverse_complement tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("reverse_complement - basic sequences", "[genome_loader][rc]") {
    SECTION("Palindromic sequence") {
        REQUIRE(reverse_complement("ACGT") == "ACGT");
    }

    SECTION("All A's become all T's") {
        REQUIRE(reverse_complement("AAAA") == "TTTT");
    }

    SECTION("All T's become all A's") {
        REQUIRE(reverse_complement("TTTT") == "AAAA");
    }

    SECTION("All G's become all C's") {
        REQUIRE(reverse_complement("GGGG") == "CCCC");
    }

    SECTION("All C's become all G's") {
        REQUIRE(reverse_complement("CCCC") == "GGGG");
    }

    SECTION("Mixed sequence") {
        // ATCG -> complement = TAGC, reverse = CGAT
        REQUIRE(reverse_complement("ATCG") == "CGAT");
    }

    SECTION("Another mixed sequence") {
        // AACG -> complement = TTGC, reverse = CGTT
        REQUIRE(reverse_complement("AACG") == "CGTT");
    }

    SECTION("20-mer CRISPR-like pattern") {
        // This is a realistic spacer sequence
        std::string spacer = "ACGTACGTACGTACGTACGT";
        std::string rc = reverse_complement(spacer);
        // ACGTACGTACGTACGTACGT -> RC is same (symmetric)
        REQUIRE(rc == "ACGTACGTACGTACGTACGT");
    }

    SECTION("Asymmetric 20-mer") {
        std::string spacer = "AAAACCCCGGGGTTTTAAAA";
        std::string rc = reverse_complement(spacer);
        // Forward:  AAAACCCCGGGGTTTTAAAA
        // Compl:    TTTTGGGGCCCCAAAATTTT
        // Reverse:  TTTTAAAACCCCGGGGTTTT
        REQUIRE(rc == "TTTTAAAACCCCGGGGTTTT");
    }
}

TEST_CASE("reverse_complement - with N wildcards", "[genome_loader][rc]") {
    SECTION("N is preserved") {
        REQUIRE(reverse_complement("ACNGT") == "ACNGT");
    }

    SECTION("All N's") {
        REQUIRE(reverse_complement("NNN") == "NNN");
    }

    SECTION("N at ends") {
        REQUIRE(reverse_complement("NACGN") == "NCGTN");
    }
}

TEST_CASE("reverse_complement - case preservation", "[genome_loader][rc]") {
    SECTION("Lowercase preserved") {
        REQUIRE(reverse_complement("acgt") == "acgt");
    }

    SECTION("Mixed case") {
        REQUIRE(reverse_complement("AcGt") == "aCgT");
    }
}

TEST_CASE("reverse_complement - empty string", "[genome_loader][rc]") {
    REQUIRE(reverse_complement("") == "");
}

TEST_CASE("reverse_complement - invalid character throws", "[genome_loader][rc]") {
    REQUIRE_THROWS_AS(reverse_complement("ACXGT"), std::invalid_argument);
    REQUIRE_THROWS_AS(reverse_complement("123"), std::invalid_argument);
}

// ───────────────────────────────────────────────────────────────────────────
// extract_genome_slice tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("extract_genome_slice - basic extraction", "[genome_loader][slice]") {
    Genome genome = load_fasta(kTestSmallFa);

    SECTION("Extract first 10 bases") {
        std::string slice = extract_genome_slice(genome, 0, 10);
        REQUIRE(slice.size() == 10);
        // test_small.fa chr1 starts with (ACGT)×14
        REQUIRE(slice == "ACGTACGTAC");
    }

    SECTION("Extract from middle of chr1") {
        std::string slice = extract_genome_slice(genome, 56, 4);
        // Position 56-59 should be in the NNNN region (positions 56-59 of chr1)
        REQUIRE(slice == "NNNN");
    }

    SECTION("Extract spanning chr1 N region") {
        std::string slice = extract_genome_slice(genome, 54, 8);
        // Positions 54-55 are GT, 56-59 are NNNN, 60-61 would be AC
        REQUIRE(slice == "GTNNNNAC");
    }
}

TEST_CASE("extract_genome_slice - out of range throws", "[genome_loader][slice]") {
    Genome genome = load_fasta(kTestSmallFa);

    REQUIRE_THROWS_AS(extract_genome_slice(genome, 0, genome.total_bases + 1),
                      std::out_of_range);
    REQUIRE_THROWS_AS(extract_genome_slice(genome, genome.total_bases, 1),
                      std::out_of_range);
}

// ───────────────────────────────────────────────────────────────────────────
// Strand enum tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Strand enum - character values", "[search_pipeline][strand]") {
    REQUIRE(static_cast<char>(Strand::PLUS) == '+');
    REQUIRE(static_cast<char>(Strand::MINUS) == '-');
}

// ───────────────────────────────────────────────────────────────────────────
// SearchConfig new fields tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("SearchConfig - strand and mismatch defaults", "[search_pipeline][config]") {
    SearchConfig config;
    REQUIRE(config.search_both_strands == true);
    REQUIRE(config.compute_mismatches == true);
}

// ───────────────────────────────────────────────────────────────────────────
// SearchHit strand comparison tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("SearchHit - sorts by position then strand", "[search_pipeline][strand]") {
    SearchHit hit_plus{100, 2, "chr1", 100, Strand::PLUS, {}};
    SearchHit hit_minus{100, 2, "chr1", 100, Strand::MINUS, {}};
    SearchHit hit_later{200, 1, "chr1", 200, Strand::PLUS, {}};

    // Different positions: position wins
    REQUIRE(hit_plus < hit_later);
    REQUIRE(hit_minus < hit_later);

    // Same position: strand comparison ('+' < '-')
    REQUIRE(hit_plus < hit_minus);
}

// ───────────────────────────────────────────────────────────────────────────
// Dual-strand search tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - finds hits on forward strand", "[search_pipeline][strand]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.search_both_strands = false;  // Forward only

    SearchResult result = search_genome(config, genome);

    // Should find exact matches of ACGT on forward strand
    REQUIRE(result.hits.size() > 0);

    // All hits should be on + strand
    for (const auto& hit : result.hits) {
        REQUIRE(hit.strand == Strand::PLUS);
    }
}

TEST_CASE("search_genome - finds hits on both strands", "[search_pipeline][strand]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGT";  // RC(ACGT) = ACGT (palindrome)
    config.threshold = 0;
    config.prefer_gpu = false;
    config.search_both_strands = true;

    SearchResult result = search_genome(config, genome);

    // For a palindromic pattern like ACGT, both strands will find the same positions
    // So we should have hits on both + and - strands
    bool has_plus = false, has_minus = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS) has_plus = true;
        if (hit.strand == Strand::MINUS) has_minus = true;
    }

    REQUIRE(has_plus);
    REQUIRE(has_minus);
}

TEST_CASE("search_genome - asymmetric pattern finds different positions per strand", "[search_pipeline][strand]") {
    // Create a simple genome with known sequences
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "AAACCCGGGAAACCCGGG"  // 18 bp
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "AAAC";  // RC = GTTT
    config.threshold = 0;
    config.prefer_gpu = false;
    config.search_both_strands = true;

    SearchResult result = search_genome(config, genome);

    // AAAC appears at positions ending at 3 and 12 on forward strand
    // GTTT (RC) doesn't appear in this genome on forward strand
    // So we should only have + strand hits

    size_t plus_count = 0, minus_count = 0;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS) plus_count++;
        if (hit.strand == Strand::MINUS) minus_count++;
    }

    REQUIRE(plus_count == 2);  // AAAC at two positions
    REQUIRE(minus_count == 0); // GTTT not in genome
}

// ───────────────────────────────────────────────────────────────────────────
// Mismatch info tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - mismatch info populated for hits", "[search_pipeline][mismatch]") {
    Genome genome = load_fasta(kTestSmallFa);

    SearchConfig config;
    config.pattern = "ACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);

    REQUIRE(result.hits.size() > 0);

    // For exact matches, CIGAR should be all-M and contain no indels.
    for (const auto& hit : result.hits) {
        if (hit.distance == 0) {
            REQUIRE(!cigar_has_indel(hit.mismatch_info.cigar));
            REQUIRE(hit.mismatch_info.aligned_sequence == "ACGT");
        }
    }
}

TEST_CASE("search_genome - mismatch positions computed correctly", "[search_pipeline][mismatch]") {
    // Create a genome with a known single mismatch
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACGTACGTNGG"  // 20bp spacer + NGG PAM (N at position 20 for exact)
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";  // 20bp pattern
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;
    config.pam.extract_length = 3;  // annotate 3'-PAM (no filter)

    SearchResult result = search_genome(config, genome);

    // Should find exact match ending at position 19
    REQUIRE(result.hits.size() >= 1);

    // Find the hit at position 19 (exact match)
    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.genome_pos == 19 && hit.strand == Strand::PLUS) {
            found = true;
            REQUIRE(hit.distance == 0);
            REQUIRE(!cigar_has_indel(hit.mismatch_info.cigar));
            REQUIRE(hit.mismatch_info.aligned_sequence == "ACGTACGTACGTACGTACGT");
            REQUIRE(hit.mismatch_info.pam_sequence == "NGG");
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - detects mismatch in pattern region", "[search_pipeline][mismatch]") {
    // 20bp spacer with mismatch at position 15 (in seed region: positions 13-20)
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACXTACGTNGG"  // X at position 14 (0-based), position 15 (1-based) in spacer
    }};

    // But wait, X is not a valid nucleotide. Let's use a real mismatch.
    // Pattern: ACGTACGTACGTACGTACGT
    // Genome:  ACGTACGTACGTACATACGT + NGG (mismatch at position 15: G->A)
    entries[0].sequence = "ACGTACGTACGTACATACGTNGG";
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";  // 20bp pattern
    config.threshold = 1;  // Allow 1 mismatch
    config.prefer_gpu = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);

    // Should find a hit with distance 1
    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.distance == 1) {
            found = true;
            // Aligned sequence should reflect the genome bases (with the mismatch).
            REQUIRE(hit.mismatch_info.aligned_sequence == "ACGTACGTACGTACATACGT");
            REQUIRE(!cigar_has_indel(hit.mismatch_info.cigar));
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - detects mismatch near pattern start", "[search_pipeline][mismatch]") {
    // 20bp spacer with mismatch at position 3 (in distal region: positions 1-12)
    // Pattern: ACGTACGTACGTACGTACGT
    // Genome:  ACATACGTACGTACGTACGT + NGG (mismatch at position 3: G->A)
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACATACGTACGTACGTACGTNGG"
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";  // 20bp pattern
    config.threshold = 1;  // Allow 1 mismatch
    config.prefer_gpu = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);

    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.distance == 1) {
            found = true;
            REQUIRE(hit.mismatch_info.aligned_sequence == "ACATACGTACGTACGTACGT");
            REQUIRE(!cigar_has_indel(hit.mismatch_info.cigar));
        }
    }
    REQUIRE(found);
}

// ───────────────────────────────────────────────────────────────────────────
// PAM validation tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - valid NGG PAM", "[search_pipeline][pam]") {
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACGTACGTAGG"  // 20bp + AGG (valid NGG)
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;
    config.pam.extract_length = 3;

    SearchResult result = search_genome(config, genome);

    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
            found = true;
            REQUIRE(hit.mismatch_info.pam_sequence == "AGG");
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - valid NAG PAM", "[search_pipeline][pam]") {
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACGTACGTCAG"  // 20bp + CAG (valid NAG)
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;
    config.pam.extract_length = 3;

    SearchResult result = search_genome(config, genome);

    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
            found = true;
            REQUIRE(hit.mismatch_info.pam_sequence == "CAG");
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - non-NGG PAM extracted", "[search_pipeline][pam]") {
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACGTACGTATG"  // 20bp + ATG (invalid PAM)
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;
    config.pam.extract_length = 3;

    SearchResult result = search_genome(config, genome);

    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
            found = true;
            REQUIRE(hit.mismatch_info.pam_sequence == "ATG");
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - various PAM sequences extracted verbatim", "[search_pipeline][pam]") {
    std::vector<std::string> pams = {"AAA", "TTT", "CCC", "ACG", "TGA", "GAC"};

    for (const auto& pam : pams) {
        std::vector<FastaEntry> entries = {{
            "test_chr",
            "ACGTACGTACGTACGTACGT" + pam
        }};
        Genome genome = encode_genome(entries);

        SearchConfig config;
        config.pattern = "ACGTACGTACGTACGTACGT";
        config.threshold = 0;
        config.prefer_gpu = false;
        config.compute_mismatches = true;
        config.pam.extract_length = 3;

        SearchResult result = search_genome(config, genome);

        for (const auto& hit : result.hits) {
            if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
                REQUIRE(hit.mismatch_info.pam_sequence == pam);
            }
        }
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Output format tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("format_hits_tsv - includes expected columns", "[search_pipeline][format]") {
    std::vector<SearchHit> hits;

    SearchHit hit;
    hit.genome_pos = 99;
    hit.distance = 1;
    hit.chrom_name = "chr1";
    hit.chrom_offset = 99;
    hit.strand = Strand::PLUS;
    hit.mismatch_info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    hit.mismatch_info.pam_sequence = "AGG";
    hit.mismatch_info.cigar = "20M";

    hits.push_back(hit);

    std::string tsv = format_hits_tsv(hits, "ACGTACGTACGTACGTACGT");

    // Header columns (new layout).
    REQUIRE(tsv.find("chrom\tstart\tend\tpattern\tdistance\tstrand\t") != std::string::npos);
    REQUIRE(tsv.find("aligned_seq\tcigar\tpam_seq\t") != std::string::npos);
    REQUIRE(tsv.find("alignment_ambiguous\tn_ambiguous_cells\t") != std::string::npos);
    REQUIRE(tsv.find("cfd_score") != std::string::npos);

    // Data row contains strand, CIGAR, and PAM.
    REQUIRE(tsv.find("\t+\t") != std::string::npos);
    REQUIRE(tsv.find("20M") != std::string::npos);
    REQUIRE(tsv.find("AGG") != std::string::npos);
}

TEST_CASE("format_hits_bed - uses actual strand", "[search_pipeline][format]") {
    std::vector<SearchHit> hits;

    SearchHit hit_plus;
    hit_plus.genome_pos = 99;
    hit_plus.distance = 0;
    hit_plus.chrom_name = "chr1";
    hit_plus.chrom_offset = 99;
    hit_plus.strand = Strand::PLUS;

    SearchHit hit_minus;
    hit_minus.genome_pos = 199;
    hit_minus.distance = 0;
    hit_minus.chrom_name = "chr1";
    hit_minus.chrom_offset = 199;
    hit_minus.strand = Strand::MINUS;

    hits.push_back(hit_plus);
    hits.push_back(hit_minus);

    std::string bed = format_hits_bed(hits, "ACGTACGTACGTACGTACGT");

    // Should contain both + and - strands
    size_t plus_count = 0, minus_count = 0;
    for (size_t i = 0; i < bed.size(); ++i) {
        if (bed[i] == '+') plus_count++;
        if (bed[i] == '-') minus_count++;
    }

    REQUIRE(plus_count >= 1);
    REQUIRE(minus_count >= 1);
}

// ───────────────────────────────────────────────────────────────────────────
// Additional edge case tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("search_genome - multiple mismatches in both regions", "[search_pipeline][mismatch]") {
    // Pattern: ACGTACGTACGTACGTACGT (20bp)
    // Genome:  ACATACATACGTACATACGT + NGG
    // Mismatches at positions 3 (distal), 7 (distal), 15 (seed)
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACATACATACGTACATACGTNGG"
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";
    config.threshold = 3;  // Allow 3 mismatches
    config.prefer_gpu = false;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);

    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.distance == 3) {
            found = true;
            // Aligned sequence should reflect all three mismatches.
            REQUIRE(hit.mismatch_info.aligned_sequence == "ACATACATACGTACATACGT");
            REQUIRE(!cigar_has_indel(hit.mismatch_info.cigar));
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - PAM at genome boundary", "[search_pipeline][pam]") {
    // Test PAM extraction when spacer is near end of genome
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACGTACGTNGG"  // Exactly 20bp + 3bp PAM
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;
    config.pam.extract_length = 3;

    SearchResult result = search_genome(config, genome);

    // Should successfully extract PAM even at the end
    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
            found = true;
            REQUIRE(hit.mismatch_info.pam_sequence == "NGG");
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - PAM beyond genome returns empty", "[search_pipeline][pam]") {
    // Test PAM extraction when spacer ends exactly at genome end
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "ACGTACGTACGTACGTACGT"  // Exactly 20bp, no PAM available
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "ACGTACGTACGTACGTACGT";
    config.threshold = 0;
    config.prefer_gpu = false;
    config.compute_mismatches = true;
    config.pam.extract_length = 3;

    SearchResult result = search_genome(config, genome);

    // PAM extraction should handle gracefully (empty or truncated PAM)
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
            // PAM should be empty or contain fewer than 3 bases
            REQUIRE(hit.mismatch_info.pam_sequence.size() < 3);
        }
    }
}

TEST_CASE("search_genome - minus strand mismatch info", "[search_pipeline][mismatch][strand]") {
    // Create genome where RC of pattern appears with a mismatch
    // Pattern: ACGTACGTACGTACGTACGT (20bp)
    // RC:      ACGTACGTACGTACGTACGT (palindromic)
    // Let's use a non-palindromic pattern
    // Pattern: AAAACCCCGGGGTTTTAAAA
    // RC:      TTTTAAAACCCCGGGGTTTT
    std::vector<FastaEntry> entries = {{
        "test_chr",
        "TTTTAAAACCCCGAGGTTTTCAG"  // RC with mismatch at position 14 (G->A in RC coords)
    }};
    Genome genome = encode_genome(entries);

    SearchConfig config;
    config.pattern = "AAAACCCCGGGGTTTTAAAA";
    config.threshold = 1;
    config.prefer_gpu = false;
    config.search_both_strands = true;
    config.compute_mismatches = true;

    SearchResult result = search_genome(config, genome);

    // Should find hit on minus strand with mismatch info
    bool found = false;
    for (const auto& hit : result.hits) {
        if (hit.strand == Strand::MINUS && hit.distance == 1) {
            found = true;
            // Verify mismatch info is populated for minus strand
            REQUIRE(!hit.mismatch_info.aligned_sequence.empty());
            REQUIRE(!hit.mismatch_info.cigar.empty());
        }
    }
    REQUIRE(found);
}

TEST_CASE("search_genome - assorted PAMs are extracted verbatim", "[search_pipeline][pam]") {
    std::vector<std::string> pams = {
        "AGG", "CGG", "GGG", "TGG",
        "AAG", "CAG", "GAG", "TAG",
    };

    for (const auto& pam : pams) {
        std::vector<FastaEntry> entries = {{
            "test_chr",
            "ACGTACGTACGTACGTACGT" + pam
        }};
        Genome genome = encode_genome(entries);

        SearchConfig config;
        config.pattern = "ACGTACGTACGTACGTACGT";
        config.threshold = 0;
        config.prefer_gpu = false;
        config.compute_mismatches = true;
        config.pam.extract_length = 3;

        SearchResult result = search_genome(config, genome);

        for (const auto& hit : result.hits) {
            if (hit.strand == Strand::PLUS && hit.genome_pos == 19) {
                REQUIRE(hit.mismatch_info.pam_sequence == pam);
            }
        }
    }
}

TEST_CASE("format_hits_tsv - exact match shows full-M CIGAR", "[search_pipeline][format]") {
    std::vector<SearchHit> hits;

    SearchHit hit;
    hit.genome_pos = 99;
    hit.distance = 0;
    hit.chrom_name = "chr1";
    hit.chrom_offset = 99;
    hit.strand = Strand::PLUS;
    hit.mismatch_info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    hit.mismatch_info.pam_sequence = "AGG";
    hit.mismatch_info.cigar = "20M";

    hits.push_back(hit);

    std::string tsv = format_hits_tsv(hits, "ACGTACGTACGTACGTACGT");

    // For an exact match the CIGAR column should hold "20M".
    REQUIRE(tsv.find("\t20M\t") != std::string::npos);
    REQUIRE(tsv.find("AGG") != std::string::npos);
}
