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

TEST_CASE("CFD score from MismatchInfo - with bulge returns 0.0", "[cfd][mismatch_info]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";

    MismatchInfo info;
    info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    info.pam_sequence = "AGG";
    info.pam_type = PamType::NGG;
    info.seed_dna_bulges = 1;  // Has a DNA bulge

    double score = cfd::compute_cfd_score(pattern, info);

    // CFD is undefined for bulges, should return 0
    REQUIRE(score == 0.0);
}

TEST_CASE("CFD score from MismatchInfo - RNA bulge returns 0.0", "[cfd][mismatch_info]") {
    std::string pattern = "ACGTACGTACGTACGTACGT";

    MismatchInfo info;
    info.aligned_sequence = "ACGTACGTACGTACGTACGT";
    info.pam_sequence = "AGG";
    info.pam_type = PamType::NGG;
    info.distal_rna_bulges = 1;  // Has an RNA bulge

    double score = cfd::compute_cfd_score(pattern, info);
    REQUIRE(score == 0.0);
}

// ───────────────────────────────────────────────────────────────────────────
// Risk tier classification tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CFD risk tier - high risk (>= 0.1)", "[cfd][risk]") {
    REQUIRE(cfd::classify_risk_tier(0.1) == 2);
    REQUIRE(cfd::classify_risk_tier(0.5) == 2);
    REQUIRE(cfd::classify_risk_tier(1.0) == 2);
}

TEST_CASE("CFD risk tier - medium risk (0.01 - 0.1)", "[cfd][risk]") {
    REQUIRE(cfd::classify_risk_tier(0.01) == 1);
    REQUIRE(cfd::classify_risk_tier(0.05) == 1);
    REQUIRE(cfd::classify_risk_tier(0.099) == 1);
}

TEST_CASE("CFD risk tier - low risk (< 0.01)", "[cfd][risk]") {
    REQUIRE(cfd::classify_risk_tier(0.0) == 0);
    REQUIRE(cfd::classify_risk_tier(0.001) == 0);
    REQUIRE(cfd::classify_risk_tier(0.009) == 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Predicted active tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CFD is_predicted_active - default threshold 0.1", "[cfd][active]") {
    REQUIRE(cfd::is_predicted_active(0.1) == true);
    REQUIRE(cfd::is_predicted_active(0.5) == true);
    REQUIRE(cfd::is_predicted_active(0.09) == false);
    REQUIRE(cfd::is_predicted_active(0.0) == false);
}

TEST_CASE("CFD is_predicted_active - custom threshold", "[cfd][active]") {
    REQUIRE(cfd::is_predicted_active(0.05, 0.05) == true);
    REQUIRE(cfd::is_predicted_active(0.04, 0.05) == false);
    REQUIRE(cfd::is_predicted_active(0.01, 0.001) == true);
}

// ───────────────────────────────────────────────────────────────────────────
// ScoringInfo struct tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("ScoringInfo - default values", "[cfd][struct]") {
    ScoringInfo info;
    REQUIRE(info.cfd_score == 0.0);
    REQUIRE(info.risk_tier == 0);
}

TEST_CASE("ScoringInfo - is_predicted_active method", "[cfd][struct]") {
    ScoringInfo info;
    info.cfd_score = 0.15;

    REQUIRE(info.is_predicted_active() == true);
    REQUIRE(info.is_predicted_active(0.1) == true);
    REQUIRE(info.is_predicted_active(0.2) == false);
}

// ───────────────────────────────────────────────────────────────────────────
// MIT Specificity Score tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("MIT specificity - empty input returns 100.0", "[mit][specificity]") {
    std::vector<double> empty_scores;
    REQUIRE(cfd::compute_mit_specificity_score(empty_scores) == 100.0);
}

TEST_CASE("MIT specificity - single off-target with low CFD", "[mit][specificity]") {
    // Formula: 100 / (100 + sum(CFD))
    // With one CFD=0.1: 100 / (100 + 0.1) = 100 / 100.1 ≈ 99.9
    std::vector<double> scores = {0.1};
    double expected = 100.0 / (100.0 + 0.1);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - single off-target with high CFD", "[mit][specificity]") {
    // With one CFD=0.9: 100 / (100 + 0.9) = 100 / 100.9 ≈ 99.1
    std::vector<double> scores = {0.9};
    double expected = 100.0 / (100.0 + 0.9);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - multiple off-targets", "[mit][specificity]") {
    // With CFD scores: [0.5, 0.3, 0.2]
    // Sum = 1.0
    // Score = 100 / (100 + 1.0) = 100 / 101 ≈ 99.01
    std::vector<double> scores = {0.5, 0.3, 0.2};
    double expected = 100.0 / (100.0 + 1.0);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - many off-targets reduces score", "[mit][specificity]") {
    // 10 off-targets with CFD=0.5 each
    // Sum = 5.0
    // Score = 100 / (100 + 5.0) = 100 / 105 ≈ 95.24
    std::vector<double> scores(10, 0.5);
    double expected = 100.0 / (100.0 + 5.0);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - high off-target burden", "[mit][specificity]") {
    // 100 off-targets with CFD=0.2 each
    // Sum = 20.0
    // Score = 100 / (100 + 20.0) = 100 / 120 ≈ 83.33
    std::vector<double> scores(100, 0.2);
    double expected = 100.0 / (100.0 + 20.0);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - perfect matches excluded by default", "[mit][specificity]") {
    // By default, perfect matches (CFD=1.0) should NOT be excluded
    // This differs from some implementations that exclude the on-target
    std::vector<double> scores = {1.0, 0.5, 0.3};
    double sum_with_perfect = 1.0 + 0.5 + 0.3;  // = 1.8
    double expected = 100.0 / (100.0 + sum_with_perfect);
    
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores, false),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - include perfect matches option", "[mit][specificity]") {
    // When include_perfect_matches=true, all scores are included
    std::vector<double> scores = {1.0, 0.5, 0.3};
    double sum_all = 1.8;
    double expected = 100.0 / (100.0 + sum_all);
    
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores, true),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - zero CFD scores", "[mit][specificity]") {
    // Off-targets with CFD=0.0 don't contribute to sum
    std::vector<double> scores = {0.0, 0.0, 0.5, 0.0};
    double expected = 100.0 / (100.0 + 0.5);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

TEST_CASE("MIT specificity - mixed CFD scores realistic scenario", "[mit][specificity]") {
    // Realistic scenario: various off-targets with different CFD scores
    std::vector<double> scores = {
        0.8, 0.6, 0.4, 0.3, 0.2,  // 5 high-scoring off-targets
        0.1, 0.1, 0.05, 0.05, 0.01 // 5 low-scoring off-targets
    };
    double sum = 0.8 + 0.6 + 0.4 + 0.3 + 0.2 + 0.1 + 0.1 + 0.05 + 0.05 + 0.01;
    double expected = 100.0 / (100.0 + sum);
    REQUIRE_THAT(cfd::compute_mit_specificity_score(scores),
                 WithinAbs(expected, 0.01));
}

// ───────────────────────────────────────────────────────────────────────────
// MIT Specificity tier classification tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("MIT tier - poor specificity (< 20)", "[mit][tier]") {
    REQUIRE(cfd::classify_specificity_tier(0.0) == 0);
    REQUIRE(cfd::classify_specificity_tier(10.0) == 0);
    REQUIRE(cfd::classify_specificity_tier(19.9) == 0);
}

TEST_CASE("MIT tier - fair specificity (20-50)", "[mit][tier]") {
    REQUIRE(cfd::classify_specificity_tier(20.0) == 1);
    REQUIRE(cfd::classify_specificity_tier(35.0) == 1);
    REQUIRE(cfd::classify_specificity_tier(49.9) == 1);
}

TEST_CASE("MIT tier - good specificity (50-80)", "[mit][tier]") {
    REQUIRE(cfd::classify_specificity_tier(50.0) == 2);
    REQUIRE(cfd::classify_specificity_tier(65.0) == 2);
    REQUIRE(cfd::classify_specificity_tier(79.9) == 2);
}

TEST_CASE("MIT tier - excellent specificity (>= 80)", "[mit][tier]") {
    REQUIRE(cfd::classify_specificity_tier(80.0) == 3);
    REQUIRE(cfd::classify_specificity_tier(90.0) == 3);
    REQUIRE(cfd::classify_specificity_tier(99.9) == 3);
    REQUIRE(cfd::classify_specificity_tier(100.0) == 3);
}

TEST_CASE("MIT tier - boundary conditions", "[mit][tier]") {
    // Test exact boundaries
    REQUIRE(cfd::classify_specificity_tier(19.999) == 0);
    REQUIRE(cfd::classify_specificity_tier(20.0) == 1);
    REQUIRE(cfd::classify_specificity_tier(49.999) == 1);
    REQUIRE(cfd::classify_specificity_tier(50.0) == 2);
    REQUIRE(cfd::classify_specificity_tier(79.999) == 2);
    REQUIRE(cfd::classify_specificity_tier(80.0) == 3);
}

