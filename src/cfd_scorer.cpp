#include <cfd_scorer.hpp>
#include <search_pipeline.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <unordered_map>

namespace cfd {

// CFD Mismatch Matrix (Doench et al. 2016, Supplementary Table 19)
//
// The matrix is indexed by:
//   - Position: 1-20 (1 = PAM-distal, 20 = PAM-proximal)
//   - RNA base: A, C, G, T (the gRNA base)
//   - DNA base: A, C, G, T (the target base, must be different from RNA base)
//
// Format: mismatch_scores[position-1][rna_base_idx][dna_base_idx]
// Base indices: A=0, C=1, G=2, T=3
//
// Data source: crisprScore R package (crisprVerse)
// https://github.com/crisprVerse/crisprScore
// ───────────────────────────────────────────────────────────────────────────

// Helper to convert base to index (A=0, C=1, G=2, T=3)
static inline int base_to_idx(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;  // Invalid base
    }
}

// CFD mismatch penalties: [position 0-19][rna 0-3][dna 0-3]
// Value of 1.0 means no penalty (match or fully tolerated mismatch)
// Diagonal entries (matches) are 1.0
static const std::array<std::array<std::array<double, 4>, 4>, 20> MISMATCH_SCORES = {{
    // Position 1 (PAM-distal)
    {{
        //       dA      dC      dG      dT      (DNA base)
        /* rA */ {1.0,    0.857,  0.9,    1.0   },
        /* rC */ {1.0,    1.0,    0.913,  1.0   },
        /* rG */ {0.9,    0.714,  1.0,    1.0   },
        /* rT */ {1.0,    0.857,  0.957,  1.0   }
    }},
    // Position 2
    {{
        /* rA */ {1.0,    0.786,  0.846,  0.727 },
        /* rC */ {0.909,  1.0,    0.696,  0.857 },
        /* rG */ {0.8,    0.692,  1.0,    0.636 },
        /* rT */ {0.846,  0.84,   0.727,  1.0   }
    }},
    // Position 3
    {{
        /* rA */ {1.0,    0.429,  0.75,   0.706 },
        /* rC */ {0.688,  1.0,    0.5,    0.429 },
        /* rG */ {0.611,  0.385,  1.0,    0.5   },
        /* rT */ {0.714,  0.5,    0.867,  1.0   }
    }},
    // Position 4
    {{
        /* rA */ {1.0,    0.353,  0.9,    0.636 },
        /* rC */ {0.8,    1.0,    0.5,    0.647 },
        /* rG */ {0.625,  0.529,  1.0,    0.364 },
        /* rT */ {0.476,  0.625,  0.842,  1.0   }
    }},
    // Position 5
    {{
        /* rA */ {1.0,    0.5,    0.867,  0.364 },
        /* rC */ {0.636,  1.0,    0.6,    1.0   },
        /* rG */ {0.72,   0.786,  1.0,    0.3   },
        /* rT */ {0.5,    0.64,   0.571,  1.0   }
    }},
    // Position 6
    {{
        /* rA */ {1.0,    0.455,  1.0,    0.714 },
        /* rC */ {0.929,  1.0,    0.5,    0.909 },
        /* rG */ {0.714,  0.682,  1.0,    0.667 },
        /* rT */ {0.867,  0.571,  0.929,  1.0   }
    }},
    // Position 7
    {{
        /* rA */ {1.0,    0.438,  1.0,    0.438 },
        /* rC */ {0.813,  1.0,    0.471,  0.688 },
        /* rG */ {0.706,  0.688,  1.0,    0.571 },
        /* rT */ {0.875,  0.588,  0.75,   1.0   }
    }},
    // Position 8
    {{
        /* rA */ {1.0,    0.429,  1.0,    0.429 },
        /* rC */ {0.875,  1.0,    0.643,  1.0   },
        /* rG */ {0.733,  0.615,  1.0,    0.625 },
        /* rT */ {0.8,    0.733,  0.65,   1.0   }
    }},
    // Position 9
    {{
        /* rA */ {1.0,    0.571,  0.643,  0.6   },
        /* rC */ {0.875,  1.0,    0.619,  0.923 },
        /* rG */ {0.667,  0.538,  1.0,    0.533 },
        /* rT */ {0.929,  0.619,  0.857,  1.0   }
    }},
    // Position 10
    {{
        /* rA */ {1.0,    0.333,  0.933,  0.882 },
        /* rC */ {0.941,  1.0,    0.389,  0.533 },
        /* rG */ {0.556,  0.4,    1.0,    0.813 },
        /* rT */ {0.857,  0.5,    0.867,  1.0   }
    }},
    // Position 11
    {{
        /* rA */ {1.0,    0.4,    1.0,    0.308 },
        /* rC */ {0.308,  1.0,    0.25,   0.667 },
        /* rG */ {0.65,   0.429,  1.0,    0.385 },
        /* rT */ {0.75,   0.4,    0.75,   1.0   }
    }},
    // Position 12
    {{
        /* rA */ {1.0,    0.263,  0.933,  0.333 },
        /* rC */ {0.538,  1.0,    0.444,  0.947 },
        /* rG */ {0.722,  0.529,  1.0,    0.385 },
        /* rT */ {0.8,    0.5,    0.714,  1.0   }
    }},
    // Position 13
    {{
        /* rA */ {1.0,    0.211,  0.923,  0.3   },
        /* rC */ {0.7,    1.0,    0.136,  0.789 },
        /* rG */ {0.652,  0.421,  1.0,    0.3   },
        /* rT */ {0.692,  0.261,  0.385,  1.0   }
    }},
    // Position 14
    {{
        /* rA */ {1.0,    0.214,  0.75,   0.533 },
        /* rC */ {0.733,  1.0,    0.0,    0.286 },
        /* rG */ {0.467,  0.429,  1.0,    0.267 },
        /* rT */ {0.619,  0.0,    0.35,   1.0   }
    }},
    // Position 15
    {{
        /* rA */ {1.0,    0.273,  0.941,  0.2   },
        /* rC */ {0.067,  1.0,    0.05,   0.273 },
        /* rG */ {0.65,   0.273,  1.0,    0.143 },
        /* rT */ {0.579,  0.05,   0.222,  1.0   }
    }},
    // Position 16
    {{
        /* rA */ {1.0,    0.0,    1.0,    0.0   },
        /* rC */ {0.308,  1.0,    0.154,  0.667 },
        /* rG */ {0.192,  0.0,    1.0,    0.0   },
        /* rT */ {0.909,  0.346,  1.0,    1.0   }
    }},
    // Position 17
    {{
        /* rA */ {1.0,    0.176,  0.933,  0.133 },
        /* rC */ {0.467,  1.0,    0.059,  0.706 },
        /* rG */ {0.176,  0.235,  1.0,    0.25  },
        /* rT */ {0.533,  0.118,  0.467,  1.0   }
    }},
    // Position 18
    {{
        /* rA */ {1.0,    0.19,   0.692,  0.5   },
        /* rC */ {0.643,  1.0,    0.133,  0.429 },
        /* rG */ {0.4,    0.476,  1.0,    0.667 },
        /* rT */ {0.667,  0.333,  0.538,  1.0   }
    }},
    // Position 19
    {{
        /* rA */ {1.0,    0.207,  0.714,  0.538 },
        /* rC */ {0.462,  1.0,    0.125,  0.276 },
        /* rG */ {0.375,  0.448,  1.0,    0.667 },
        /* rT */ {0.286,  0.25,   0.429,  1.0   }
    }},
    // Position 20 (PAM-proximal)
    // Raw data: TA20=0.5625, CA20=0.5, GA20=0.9375, TG20=0.176, etc.
    {{
        /* rA */ {1.0,    0.227,  0.765,  0.6   },
        /* rC */ {0.3,    1.0,    0.059,  0.091 },
        /* rG */ {0.938,  0.429,  1.0,    0.7   },
        /* rT */ {0.563,  0.091,  0.176,  1.0   }
    }}
}};

// CFD PAM Penalties (Doench et al. 2016)
//
// PAM dinucleotide scores (the 2nd and 3rd bases after the spacer)
// For a canonical NGG PAM, the dinucleotide is "GG"
//
// Data source: crisprScore R package (crisprVerse)
// ───────────────────────────────────────────────────────────────────────────

// PAM penalties indexed by dinucleotide
// First base index * 4 + second base index (A=0, C=1, G=2, T=3)
static const std::array<double, 16> PAM_SCORES = {{
    // AA    AC    AG    AT
    0.0,    0.0,  0.259, 0.0,
    // CA    CC    CG    CT
    0.0,    0.0,  0.107, 0.0,
    // GA    GC    GG    GT
    0.069,  0.022, 1.0,  0.016,
    // TA    TC    TG    TT
    0.0,    0.0,  0.039, 0.0
}};

// Implementation

double get_mismatch_penalty(uint8_t position, char rna_base, char dna_base) {
    // Validate position (1-20)
    if (position < 1 || position > 20) {
        return 1.0;  // No penalty for out-of-range positions
    }

    // Convert bases to indices
    int rna_idx = base_to_idx(rna_base);
    int dna_idx = base_to_idx(dna_base);

    // Invalid bases get no penalty
    if (rna_idx < 0 || dna_idx < 0) {
        return 1.0;
    }

    // Matching bases have no penalty
    if (rna_idx == dna_idx) {
        return 1.0;
    }

    // Look up the mismatch penalty
    return MISMATCH_SCORES[position - 1][rna_idx][dna_idx];
}

double get_pam_penalty(const std::string& pam_dinucleotide) {
    if (pam_dinucleotide.size() < 2) {
        return 0.0;  // Incomplete PAM
    }

    int idx1 = base_to_idx(pam_dinucleotide[0]);
    int idx2 = base_to_idx(pam_dinucleotide[1]);

    if (idx1 < 0 || idx2 < 0) {
        return 0.0;  // Invalid bases
    }

    return PAM_SCORES[idx1 * 4 + idx2];
}

// PamType lacks the actual dinucleotide, so OTHER is scored as 0.0 (inactive)
// — callers that know the dinucleotide should call the string overload.
double get_pam_penalty(PamType pam_type) {
    switch (pam_type) {
        case PamType::NGG: return get_pam_penalty("GG");
        case PamType::NAG: return get_pam_penalty("AG");
        case PamType::OTHER:
        case PamType::INCOMPLETE:
        default:           return 0.0;
    }
}

double compute_cfd_score(const std::string& pattern,
                         const std::string& aligned_seq,
                         const std::string& pam_sequence) {
    // Validate inputs
    if (pattern.empty() || aligned_seq.empty()) {
        return 0.0;
    }

    // Use the shorter length to avoid out-of-bounds access, then cap at the
    // CFD model's trained spacer length.
    size_t len = std::min(pattern.size(), aligned_seq.size());
    if (len > CFD_MAX_SPACER_LENGTH) {
        len = CFD_MAX_SPACER_LENGTH;
    }

    // Compute product of mismatch penalties
    double score = 1.0;

    for (size_t i = 0; i < len; ++i) {
        char rna_base = pattern[i];
        char dna_base = aligned_seq[i];

        // N bases are treated as matches (no penalty)
        char rna_upper = static_cast<char>(std::toupper(static_cast<unsigned char>(rna_base)));
        char dna_upper = static_cast<char>(std::toupper(static_cast<unsigned char>(dna_base)));

        if (rna_upper == 'N' || dna_upper == 'N') {
            continue;  // No penalty for N
        }

        if (rna_upper != dna_upper) {
            // Mismatch - apply penalty
            // Position is 1-based, from PAM-distal to PAM-proximal
            uint8_t position = static_cast<uint8_t>(i + 1);
            double penalty = get_mismatch_penalty(position, rna_base, dna_base);
            score *= penalty;
        }
    }

    // Apply PAM penalty
    if (pam_sequence.size() >= 3) {
        // Extract dinucleotide (2nd and 3rd bases of PAM)
        std::string dinucleotide = pam_sequence.substr(1, 2);
        score *= get_pam_penalty(dinucleotide);
    } else if (pam_sequence.size() == 2) {
        // Assume we have just the dinucleotide
        score *= get_pam_penalty(pam_sequence);
    } else {
        // Incomplete PAM
        score *= 0.0;
    }

    return score;
}

double compute_cfd_score(const std::string& pattern,
                         const MismatchInfo& mismatch_info) {
    // CFD model only handles mismatches, not bulges (indels)
    // If there are any bulges, CFD is not applicable
    if (mismatch_info.has_dna_bulge() || mismatch_info.has_rna_bulge()) {
        return 0.0;  // CFD undefined for indels
    }

    // Use aligned sequence for comparison
    return compute_cfd_score(pattern,
                             mismatch_info.aligned_sequence,
                             mismatch_info.pam_sequence);
}

uint8_t classify_risk_tier(double cfd_score) {
    if (cfd_score >= CFD_HIGH_RISK_THRESHOLD) {
        return 2;  // High risk
    } else if (cfd_score >= CFD_MEDIUM_RISK_THRESHOLD) {
        return 1;  // Medium risk
    } else {
        return 0;  // Low risk
    }
}

bool is_predicted_active(double cfd_score, double threshold) {
    return cfd_score >= threshold;
}

// MIT Specificity Score implementation

double compute_mit_specificity_score(
    const std::vector<double>& cfd_scores,
    bool include_perfect_matches) {
    
    // Empty input = perfect specificity
    if (cfd_scores.empty()) {
        return 100.0;
    }

    // Sum up all CFD scores (optionally excluding perfect matches)
    double sum = 0.0;
    for (double cfd : cfd_scores) {
        // Skip perfect matches if requested
        if (!include_perfect_matches && cfd >= 0.9999) {
            continue;
        }
        sum += cfd;
    }

    // MIT formula: 100 / (100 + sum(CFD))
    return 100.0 / (100.0 + sum);
}

uint8_t classify_specificity_tier(double mit_score) {
    if (mit_score >= MIT_EXCELLENT_THRESHOLD) {
        return 3;  // Excellent specificity
    } else if (mit_score >= MIT_GOOD_THRESHOLD) {
        return 2;  // Good specificity
    } else if (mit_score >= MIT_FAIR_THRESHOLD) {
        return 1;  // Fair specificity
    } else {
        return 0;  // Poor specificity
    }
}

} // namespace cfd

