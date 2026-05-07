#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cfd_scorer.hpp>
#include <search_pipeline.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ───────────────────────────────────────────────────────────────────────────
// Mismatch penalty tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CFD mismatch penalty - matching bases return 1.0", "[cfd][mismatch]") {
    // Matching bases should have no penalty (score 1.0)
    REQUIRE(cfd::get_mismatch_penalty(1, 'A', 'A') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(10, 'C', 'C') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(20, 'G', 'G') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(15, 'T', 'T') == 1.0);
}

TEST_CASE("CFD mismatch penalty - out of range positions return 1.0", "[cfd][mismatch]") {
    // Position 0 and positions > 20 should return 1.0 (no penalty)
    REQUIRE(cfd::get_mismatch_penalty(0, 'A', 'T') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(21, 'A', 'T') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(255, 'G', 'C') == 1.0);
}

TEST_CASE("CFD mismatch penalty - invalid bases return 1.0", "[cfd][mismatch]") {
    // Invalid bases (N, etc.) should return 1.0
    REQUIRE(cfd::get_mismatch_penalty(10, 'N', 'A') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(10, 'A', 'N') == 1.0);
    REQUIRE(cfd::get_mismatch_penalty(10, 'X', 'A') == 1.0);
}

TEST_CASE("CFD mismatch penalty - case insensitive", "[cfd][mismatch]") {
    // Should handle lowercase bases
    double upper = cfd::get_mismatch_penalty(10, 'A', 'T');
    double lower = cfd::get_mismatch_penalty(10, 'a', 't');
    REQUIRE(upper == lower);
}

TEST_CASE("CFD mismatch penalty - known values from Doench 2016", "[cfd][mismatch]") {
    // Verify some known values from the Doench 2016 paper

    // Position 1 (PAM-distal): Most mismatches are tolerated (close to 1.0)
    REQUIRE_THAT(cfd::get_mismatch_penalty(1, 'T', 'A'), WithinAbs(1.0, 0.01));
    REQUIRE_THAT(cfd::get_mismatch_penalty(1, 'G', 'A'), WithinAbs(0.9, 0.01));

    // Position 16: Position with some zeros (sensitive position)
    REQUIRE_THAT(cfd::get_mismatch_penalty(16, 'A', 'T'), WithinAbs(0.0, 0.01));
    REQUIRE_THAT(cfd::get_mismatch_penalty(16, 'G', 'C'), WithinAbs(0.0, 0.01));

    // Position 20 (PAM-proximal): Variable sensitivity
    REQUIRE_THAT(cfd::get_mismatch_penalty(20, 'G', 'A'), WithinAbs(0.938, 0.01));
    REQUIRE_THAT(cfd::get_mismatch_penalty(20, 'T', 'G'), WithinAbs(0.176, 0.01));
}

TEST_CASE("CFD matrix regression — cells that exposed the v0.7.0 transposition bug",
          "[cfd][mismatch][regression]") {
    // These cells were flipped in Stomata v0.6.5/v0.7.0 vs. the crisprScore
    // reference. Ground-truth values taken from crisprScore::getCFDScores()
    // (Bioconductor, commit state as of 2026-04-23). If any of these fail
    // again the matrix has been corrupted — recheck against crisprScore.
    REQUIRE_THAT(cfd::get_mismatch_penalty(3,  'C', 'A'), WithinAbs(0.866667, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(3,  'C', 'T'), WithinAbs(0.687500, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(4,  'A', 'G'), WithinAbs(0.625000, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(4,  'G', 'A'), WithinAbs(0.900000, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(5,  'C', 'T'), WithinAbs(0.636364, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(5,  'T', 'C'), WithinAbs(1.000000, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(6,  'A', 'G'), WithinAbs(0.714286, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(6,  'G', 'A'), WithinAbs(1.000000, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(13, 'C', 'G'), WithinAbs(0.136364, 0.001));
    REQUIRE_THAT(cfd::get_mismatch_penalty(17, 'C', 'A'), WithinAbs(0.466667, 0.001));
}

TEST_CASE("CFD compound score — crisprScore reference pairs", "[cfd][score][regression]") {
    // End-to-end CFD scores on five (spacer, protospacer, PAM) pairs.
    // Ground-truth from crisprScore::getCFDScores(). Asymmetric purine- and
    // pyrimidine-transition mismatches are included to catch the v0.6.5
    // transposition bug in compound form.
    struct Case {
        std::string spacer;
        std::string proto;
        std::string pam;
        double expected;
    };
    std::vector<Case> cases = {
        {"GTCACCAATCCTGTCCCTAG", "GTCACCAATCCTGTCCCTAG", "AGG", 1.0000000},
        {"GTCACCAATCCTGTCCCTAG", "GTAACCAATCCTGTCCATAG", "TGG", 0.4044444},
        {"GTCACCAATCCTGTCCCTAG", "GCCACCAGTCCGGTCCCAAG", "CGG", 0.2095238},
        {"GAGTCCGAGCAGAAGAAGAA", "GAGTCCGAGCAGAAGAAGAA", "AGG", 1.0000000},
        {"GAGTCCGAGCAGAAGAAGAA", "GAGTCAGAGCAGAAGAAGAA", "AGG", 0.9285714},
    };
    for (const auto& c : cases) {
        double actual = cfd::compute_cfd_score(c.spacer, c.proto, c.pam);
        REQUIRE_THAT(actual, WithinAbs(c.expected, 0.001));
    }
}

TEST_CASE("CFD mismatch penalty - seed vs distal region differences", "[cfd][mismatch]") {
    // Mismatches in seed region (positions 13-20) generally have lower penalties
    // than distal region (positions 1-12) for some mismatch types

    // rG:dC mismatch at different positions
    double distal = cfd::get_mismatch_penalty(5, 'G', 'C');  // Position 5
    double seed = cfd::get_mismatch_penalty(17, 'G', 'C');   // Position 17

    // Both should be less than 1.0 (there is a penalty)
    REQUIRE(distal < 1.0);
    REQUIRE(seed < 1.0);
}

// ───────────────────────────────────────────────────────────────────────────
// PAM penalty tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CFD PAM penalty - NGG is canonical (1.0)", "[cfd][pam]") {
    REQUIRE(cfd::get_pam_penalty(PamType::NGG) == 1.0);
    REQUIRE(cfd::get_pam_penalty("GG") == 1.0);
}

TEST_CASE("CFD PAM penalty - NAG is alternative (~0.26)", "[cfd][pam]") {
    REQUIRE_THAT(cfd::get_pam_penalty(PamType::NAG), WithinAbs(0.259, 0.01));
    REQUIRE_THAT(cfd::get_pam_penalty("AG"), WithinAbs(0.259, 0.01));
}

TEST_CASE("CFD PAM penalty - OTHER returns 0.0", "[cfd][pam]") {
    REQUIRE(cfd::get_pam_penalty(PamType::OTHER) == 0.0);
}

TEST_CASE("CFD PAM penalty - INCOMPLETE returns 0.0", "[cfd][pam]") {
    REQUIRE(cfd::get_pam_penalty(PamType::INCOMPLETE) == 0.0);
    REQUIRE(cfd::get_pam_penalty("G") == 0.0);  // Only one base
    REQUIRE(cfd::get_pam_penalty("") == 0.0);   // Empty
}

TEST_CASE("CFD PAM penalty - dinucleotide lookup", "[cfd][pam]") {
    // Check various dinucleotide PAMs
    REQUIRE_THAT(cfd::get_pam_penalty("CG"), WithinAbs(0.107, 0.01));
    REQUIRE_THAT(cfd::get_pam_penalty("GA"), WithinAbs(0.069, 0.01));
    REQUIRE_THAT(cfd::get_pam_penalty("GC"), WithinAbs(0.022, 0.01));
    REQUIRE_THAT(cfd::get_pam_penalty("GT"), WithinAbs(0.016, 0.01));
    REQUIRE_THAT(cfd::get_pam_penalty("TG"), WithinAbs(0.039, 0.01));

    // Non-functional PAMs (score 0)
    REQUIRE(cfd::get_pam_penalty("AA") == 0.0);
    REQUIRE(cfd::get_pam_penalty("TT") == 0.0);
    REQUIRE(cfd::get_pam_penalty("CC") == 0.0);
}

// ───────────────────────────────────────────────────────────────────────────
// Full CFD score computation tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CFD score - perfect match with NGG PAM", "[cfd][score]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";  // 20 bp
    std::string aligned = "ACGTACGTACGTACGTACGT";  // Perfect match
    std::string pam = "AGG";  // NGG PAM

    double score = cfd::compute_cfd_score(pattern, aligned, pam);

    // Perfect match with NGG should be 1.0
    REQUIRE_THAT(score, WithinAbs(1.0, 0.001));
}

TEST_CASE("CFD score - perfect match with NAG PAM", "[cfd][score]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string aligned = "ACGTACGTACGTACGTACGT";
    std::string pam = "AAG";  // NAG PAM

    double score = cfd::compute_cfd_score(pattern, aligned, pam);

    // Perfect match with NAG should be ~0.26
    REQUIRE_THAT(score, WithinAbs(0.259, 0.01));
}

TEST_CASE("CFD score - single mismatch", "[cfd][score]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string aligned = "TCGTACGTACGTACGTACGT";  // A->T at position 1
    std::string pam = "AGG";

    double score = cfd::compute_cfd_score(pattern, aligned, pam);

    // Should be mismatch_penalty(1, A, T) * 1.0 (NGG)
    double expected = cfd::get_mismatch_penalty(1, 'A', 'T') * 1.0;
    REQUIRE_THAT(score, WithinAbs(expected, 0.001));
}

TEST_CASE("CFD score - multiple mismatches multiply", "[cfd][score]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string aligned = "TCGTACGTACGTACGTACGA";  // Mismatches at pos 1 and 20
    std::string pam = "AGG";

    double score = cfd::compute_cfd_score(pattern, aligned, pam);

    // Should be product of individual penalties
    double penalty1 = cfd::get_mismatch_penalty(1, 'A', 'T');
    double penalty20 = cfd::get_mismatch_penalty(20, 'T', 'A');
    double expected = penalty1 * penalty20 * 1.0;  // * NGG PAM

    REQUIRE_THAT(score, WithinAbs(expected, 0.001));
}

TEST_CASE("CFD score - empty input returns 0.0", "[cfd][score]") {
    REQUIRE(cfd::compute_cfd_score("", "ACGT", "AGG") == 0.0);
    REQUIRE(cfd::compute_cfd_score("ACGT", "", "AGG") == 0.0);
}

TEST_CASE("CFD score - incomplete PAM returns 0.0", "[cfd][score]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string aligned = "ACGTACGTACGTACGTACGT";

    REQUIRE(cfd::compute_cfd_score(pattern, aligned, "") == 0.0);
    REQUIRE(cfd::compute_cfd_score(pattern, aligned, "A") == 0.0);
}

TEST_CASE("CFD score - N bases are treated as matches", "[cfd][score]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";
    std::string aligned = "NCGTACGTACGTACGTACGT";  // N at position 1
    std::string pam = "AGG";

    double score = cfd::compute_cfd_score(pattern, aligned, pam);

    // N should be treated as a match (no penalty)
    REQUIRE_THAT(score, WithinAbs(1.0, 0.001));
}

// ───────────────────────────────────────────────────────────────────────────
// MismatchInfo-based CFD score tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CFD score from MismatchInfo - perfect match", "[cfd][mismatch_info]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";

    MismatchInfo info;
    info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    info.pam_sequence = "AGG";
    info.pam_type = PamType::NGG;

    double score = cfd::compute_cfd_score(pattern, info);
    REQUIRE_THAT(score, WithinAbs(1.0, 0.001));
}

TEST_CASE("CFD score from MismatchInfo - DNA bulge (insertion) returns 0.0", "[cfd][mismatch_info]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";

    MismatchInfo info;
    info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    info.pam_sequence = "AGG";
    info.cigar = "10M1I9M";  // Insertion = DNA bulge in aligner terms

    double score = cfd::compute_cfd_score(pattern, info);

    // CFD is undefined for indels, should return 0
    REQUIRE(score == 0.0);
}

TEST_CASE("CFD score from MismatchInfo - RNA bulge (deletion) returns 0.0", "[cfd][mismatch_info]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";

    MismatchInfo info;
    info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    info.pam_sequence = "AGG";
    info.cigar = "5M1D14M";  // Deletion = RNA bulge

    double score = cfd::compute_cfd_score(pattern, info);
    REQUIRE(score == 0.0);
}

// ───────────────────────────────────────────────────────────────────────────
// ScoringInfo struct tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ScoringInfo - default values", "[cfd][struct]") {
    ScoringInfo info;
    REQUIRE(info.cfd_score == 0.0);
}

