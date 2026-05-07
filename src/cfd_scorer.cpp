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
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.857143, 1.000000, 1.000000},
        /* rC */ {1.000000, 1.000000, 0.913043, 1.000000},
        /* rG */ {0.900000, 0.714286, 1.000000, 1.000000},
        /* rT */ {1.000000, 0.857143, 0.956522, 1.000000}
    }},
    // Position 2
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.785714, 0.800000, 0.727273},
        /* rC */ {0.727273, 1.000000, 0.695652, 0.909091},
        /* rG */ {0.846154, 0.692308, 1.000000, 0.636364},
        /* rT */ {0.846154, 0.857143, 0.840000, 1.000000}
    }},
    // Position 3
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.428571, 0.611111, 0.705882},
        /* rC */ {0.866667, 1.000000, 0.500000, 0.687500},
        /* rG */ {0.750000, 0.384615, 1.000000, 0.500000},
        /* rT */ {0.714286, 0.428571, 0.500000, 1.000000}
    }},
    // Position 4
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.352941, 0.625000, 0.636364},
        /* rC */ {0.842105, 1.000000, 0.500000, 0.800000},
        /* rG */ {0.900000, 0.529412, 1.000000, 0.363636},
        /* rT */ {0.476190, 0.647059, 0.625000, 1.000000}
    }},
    // Position 5
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.500000, 0.720000, 0.363636},
        /* rC */ {0.571429, 1.000000, 0.600000, 0.636364},
        /* rG */ {0.866667, 0.785714, 1.000000, 0.300000},
        /* rT */ {0.500000, 1.000000, 0.640000, 1.000000}
    }},
    // Position 6
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.454545, 0.714286, 0.714286},
        /* rC */ {0.928571, 1.000000, 0.500000, 0.928571},
        /* rG */ {1.000000, 0.681818, 1.000000, 0.666667},
        /* rT */ {0.866667, 0.909091, 0.571429, 1.000000}
    }},
    // Position 7
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.437500, 0.705882, 0.437500},
        /* rC */ {0.750000, 1.000000, 0.470588, 0.812500},
        /* rG */ {1.000000, 0.687500, 1.000000, 0.571429},
        /* rT */ {0.875000, 0.687500, 0.588235, 1.000000}
    }},
    // Position 8
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.428571, 0.733333, 0.428571},
        /* rC */ {0.650000, 1.000000, 0.642857, 0.875000},
        /* rG */ {1.000000, 0.615385, 1.000000, 0.625000},
        /* rT */ {0.800000, 1.000000, 0.733333, 1.000000}
    }},
    // Position 9
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.571429, 0.666667, 0.600000},
        /* rC */ {0.857143, 1.000000, 0.619048, 0.875000},
        /* rG */ {0.642857, 0.538462, 1.000000, 0.533333},
        /* rT */ {0.928571, 0.923077, 0.619048, 1.000000}
    }},
    // Position 10
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.333333, 0.555556, 0.882353},
        /* rC */ {0.866667, 1.000000, 0.388889, 0.941176},
        /* rG */ {0.933333, 0.400000, 1.000000, 0.812500},
        /* rT */ {0.857143, 0.533333, 0.500000, 1.000000}
    }},
    // Position 11
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.400000, 0.650000, 0.307692},
        /* rC */ {0.750000, 1.000000, 0.250000, 0.307692},
        /* rG */ {1.000000, 0.428571, 1.000000, 0.384615},
        /* rT */ {0.750000, 0.666667, 0.400000, 1.000000}
    }},
    // Position 12
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.263158, 0.722222, 0.333333},
        /* rC */ {0.714286, 1.000000, 0.444444, 0.538462},
        /* rG */ {0.933333, 0.529412, 1.000000, 0.384615},
        /* rT */ {0.800000, 0.947368, 0.500000, 1.000000}
    }},
    // Position 13
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.210526, 0.652174, 0.300000},
        /* rC */ {0.384615, 1.000000, 0.136364, 0.700000},
        /* rG */ {0.923077, 0.421053, 1.000000, 0.300000},
        /* rT */ {0.692308, 0.789474, 0.260870, 1.000000}
    }},
    // Position 14
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.214286, 0.466667, 0.533333},
        /* rC */ {0.350000, 1.000000, 0.000000, 0.733333},
        /* rG */ {0.750000, 0.428571, 1.000000, 0.266667},
        /* rT */ {0.619048, 0.285714, 0.000000, 1.000000}
    }},
    // Position 15
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.272727, 0.650000, 0.200000},
        /* rC */ {0.222222, 1.000000, 0.050000, 0.066667},
        /* rG */ {0.941176, 0.272727, 1.000000, 0.142857},
        /* rT */ {0.578947, 0.272727, 0.050000, 1.000000}
    }},
    // Position 16
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.000000, 0.192308, 0.000000},
        /* rC */ {1.000000, 1.000000, 0.153846, 0.307692},
        /* rG */ {1.000000, 0.000000, 1.000000, 0.000000},
        /* rT */ {0.909091, 0.666667, 0.346154, 1.000000}
    }},
    // Position 17
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.176471, 0.176471, 0.133333},
        /* rC */ {0.466667, 1.000000, 0.058824, 0.466667},
        /* rG */ {0.933333, 0.235294, 1.000000, 0.250000},
        /* rT */ {0.533333, 0.705882, 0.117647, 1.000000}
    }},
    // Position 18
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.190476, 0.400000, 0.500000},
        /* rC */ {0.538462, 1.000000, 0.133333, 0.642857},
        /* rG */ {0.692308, 0.476190, 1.000000, 0.666667},
        /* rT */ {0.666667, 0.428571, 0.333333, 1.000000}
    }},
    // Position 19
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.206897, 0.375000, 0.538462},
        /* rC */ {0.428571, 1.000000, 0.125000, 0.461538},
        /* rG */ {0.714286, 0.448276, 1.000000, 0.666667},
        /* rT */ {0.285714, 0.275862, 0.250000, 1.000000}
    }},
    // Position 20 (PAM-proximal)
    {{
        //       dA         dC         dG         dT        (DNA base)
        /* rA */ {1.000000, 0.227273, 0.764706, 0.600000},
        /* rC */ {0.500000, 1.000000, 0.058824, 0.300000},
        /* rG */ {0.937500, 0.428571, 1.000000, 0.700000},
        /* rT */ {0.562500, 0.090909, 0.176471, 1.000000}
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
    // CFD model only handles mismatches, not bulges (indels).
    if (cigar_has_indel(mismatch_info.cigar)) {
        return 0.0;
    }

    return compute_cfd_score(pattern,
                             mismatch_info.aligned_sequence,
                             mismatch_info.pam_sequence);
}

} // namespace cfd

