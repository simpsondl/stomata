#pragma once

#include <string>
#include <vector>
#include <cstdint>

// Forward declarations
struct SearchHit;
struct MismatchInfo;
enum class PamType : uint8_t;

// CFD (Cutting Frequency Determination) Score
//
// The CFD score predicts off-target cleavage activity based on the Doench et al.
// 2016 paper "Optimized sgRNA design to maximize activity and minimize off-target
// effects of CRISPR-Cas9" (Nature Biotechnology).
//
// CFD scores range from 0.0 to 1.0:
//   - 1.0 = predicted to cleave as efficiently as the on-target
//   - 0.0 = predicted to have no cleavage activity
//
// The score is computed as the product of:
//   1. Individual mismatch penalties (position and nucleotide-specific)
//   2. PAM penalty (for non-NGG PAMs)
//
// References:
//   - Doench et al. (2016) Nature Biotechnology 34:184-191
//   - Supplementary Table 19 (mismatch scores)
//   - crisprScore R package (crisprVerse)
// ───────────────────────────────────────────────────────────────────────────

namespace cfd {

// CFD model was trained on 20 nt SpCas9 spacers (Doench 2016). Positions beyond
// 20 have no published penalty table, so we truncate silently.
inline constexpr size_t CFD_MAX_SPACER_LENGTH = 20;

// Risk tier cutoffs (fraction of on-target activity).
inline constexpr double CFD_HIGH_RISK_THRESHOLD   = 0.1;   // >=10% activity
inline constexpr double CFD_MEDIUM_RISK_THRESHOLD = 0.01;  // >=1% activity

// MIT specificity tier cutoffs (100 / (100 + ΣCFD)).
inline constexpr double MIT_EXCELLENT_THRESHOLD = 80.0;
inline constexpr double MIT_GOOD_THRESHOLD      = 50.0;
inline constexpr double MIT_FAIR_THRESHOLD      = 20.0;

// Core scoring functions

// Get the CFD penalty for a single mismatch at a specific position.
//
// Parameters:
//   position  - 1-based position in the spacer (1 = PAM-distal, 20 = PAM-proximal)
//   rna_base  - The RNA/gRNA base at this position (A, C, G, T)
//   dna_base  - The DNA/target base at this position (A, C, G, T)
//
// Returns:
//   Penalty score between 0.0 and 1.0
//   - 1.0 means fully tolerated (no effect on activity)
//   - 0.0 means completely abolishes activity
//   - Returns 1.0 for matching bases (no mismatch)
//   - Returns 1.0 for positions outside 1-20 range
//   - Returns 1.0 for invalid/unknown bases (N, etc.)
double get_mismatch_penalty(uint8_t position, char rna_base, char dna_base);

// Get the CFD penalty for a PAM sequence.
//
// Parameters:
//   pam_type - The classified PAM type (NGG, NAG, OTHER, INCOMPLETE)
//
// Returns:
//   Penalty score between 0.0 and 1.0
//   - 1.0 for NGG (canonical PAM)
//   - ~0.26 for NAG (alternative PAM)
//   - Variable for OTHER based on actual dinucleotide
//   - 0.0 for INCOMPLETE
double get_pam_penalty(PamType pam_type);

// Get the CFD penalty for a specific PAM dinucleotide sequence.
//
// Parameters:
//   pam_dinucleotide - The 2nd and 3rd bases of the PAM (e.g., "GG" from "NGG")
//
// Returns:
//   Penalty score between 0.0 and 1.0
double get_pam_penalty(const std::string& pam_dinucleotide);

// High-level scoring functions

// Compute the CFD score for an off-target hit.
//
// The CFD score is calculated as:
//   CFD = (product of mismatch penalties) * (PAM penalty)
//
// Parameters:
//   pattern       - The original spacer/gRNA sequence (used to get RNA bases)
//   aligned_seq   - The aligned genomic sequence at the hit (DNA bases)
//   pam_sequence  - The 3-base PAM sequence (e.g., "AGG", "AAG")
//
// Returns:
//   CFD score between 0.0 and 1.0
//
// Notes:
//   - Sequences should be the same length (typically 20bp for SpCas9)
//   - N bases in either sequence are treated as matches (no penalty)
//   - Only mismatches are penalized; bulges (indels) are currently scored as 0
double compute_cfd_score(const std::string& pattern,
                         const std::string& aligned_seq,
                         const std::string& pam_sequence);

// Compute the CFD score using a SearchHit's MismatchInfo.
//
// This is the primary interface for scoring hits from the search pipeline.
//
// Parameters:
//   pattern       - The original spacer/gRNA sequence
//   mismatch_info - The detailed mismatch information from a SearchHit
//
// Returns:
//   CFD score between 0.0 and 1.0
//
// Notes:
//   - If the hit contains DNA or RNA bulges, CFD score is set to 0.0
//     (CFD model only handles mismatches, not indels)
//   - Perfect matches (distance 0) return 1.0 * PAM_penalty
double compute_cfd_score(const std::string& pattern,
                         const MismatchInfo& mismatch_info);

// Utility functions

// Classify a CFD score into a risk tier.
//
// Returns:
//   0 = Low risk (CFD < 0.01) - very unlikely to cause off-target cleavage
//   1 = Medium risk (0.01 <= CFD < 0.1)
//   2 = High risk (CFD >= 0.1) - significant off-target potential
uint8_t classify_risk_tier(double cfd_score);

// Check if a CFD score exceeds the "active" threshold
// (default: CFD_HIGH_RISK_THRESHOLD = 10% of on-target activity).
bool is_predicted_active(double cfd_score, double threshold = CFD_HIGH_RISK_THRESHOLD);

// MIT Specificity Score
//
// The MIT Specificity Score is an aggregate metric that predicts guide RNA
// specificity by considering all potential off-targets in the genome.
// Developed by the Hsu/Zhang lab at MIT, it's used in their CRISPR design tool.
//
// The formula is:
//   Specificity = 100 / (100 + Σ(CFD_i))
//
// Where CFD_i is the CFD score of each off-target site (excluding on-target).
//
// Specificity scores range from 0 to 100:
//   - 100 = perfect specificity (no predicted off-targets)
//   - >50 = generally acceptable guide
//   - <50 = high off-target burden, consider alternative guide
//   - 0 = approaches as off-target burden increases
//
// References:
//   - Hsu et al. (2013) Nature Biotechnology 31:827-832
//   - MIT CRISPR Design Tool: http://crispr.mit.edu
// ───────────────────────────────────────────────────────────────────────────

// Compute MIT Specificity Score from a collection of CFD scores.
//
// Parameters:
//   cfd_scores - Vector of CFD scores for all off-target hits
//   include_perfect_matches - Whether to include CFD=1.0 scores in the sum
//                             (typically false, as perfect matches are
//                             excluded from specificity calculation)
//
// Returns:
//   Specificity score between 0.0 and 100.0
//   - 100.0 = no off-targets (or empty input)
//   - Lower scores indicate higher off-target burden
double compute_mit_specificity_score(
    const std::vector<double>& cfd_scores,
    bool include_perfect_matches = false);

// Classify a MIT specificity score into quality tiers.
//
// Returns:
//   0 = Poor specificity (score < 20) - high off-target burden
//   1 = Fair specificity (20 <= score < 50) - moderate off-target burden
//   2 = Good specificity (50 <= score < 80) - acceptable for most uses
//   3 = Excellent specificity (score >= 80) - very few predicted off-targets
uint8_t classify_specificity_tier(double mit_score);

} // namespace cfd
