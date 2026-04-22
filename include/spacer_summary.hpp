#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "annotation_bed.hpp"
#include "search_pipeline.hpp"

// Per-spacer summary layer.
//
// Aggregate CFD (aCFD) follows CRISPick 2025 (Doench lab):
//   - SDR = guide positions 4-20 (1-indexed, PAM-proximal; PAM-distal positions 1-3 tolerated)
//   - aCFD = sum of per-hit CFD scores for hits with SDR mismatches <= 1
//   - Default promiscuity threshold: 4.8
//   - Applicable only to Cas9 with 20 bp spacers; marked NA otherwise
//
// Genomic feature annotation: pass zero or more named BedForest pointers via
// AnnotationBed. Results appear as n_<label>_overlaps columns in TSV output.

namespace spacer {

// Default aCFD threshold for promiscuity classification (CRISPick 2025).
inline constexpr double DEFAULT_ACFD_THRESHOLD = 4.8;

// SDR bounds (1-indexed, inclusive). SDR = positions 4-20.
inline constexpr int SDR_START = 4;
inline constexpr int SDR_END   = 20;

// Return value from count_sdr_mismatches when the alignment has indels.
inline constexpr int SDR_INDEL = -1;

// Count mismatches between pattern and aligned_seq restricted to SDR positions
// (guide positions 4-20, 1-indexed, where 1 = PAM-distal / 5' end).
//
// Returns SDR_INDEL (-1) if the CIGAR contains any I or D.
// Positions are counted on the spacer string: index 0 = position 1.
int count_sdr_mismatches(const std::string& pattern,
                         const std::string& aligned_seq,
                         const std::string& cigar);

// Summary for one spacer.
struct SpacerSummaryRow {
    std::string spacer_name;
    std::string spacer_sequence;

    // Off-target counts by edit distance (all hits, including any on-target d=0)
    size_t n_d0  = 0;   // distance == 0
    size_t n_d1  = 0;   // distance == 1
    size_t n_d2  = 0;   // distance == 2
    size_t n_d3p = 0;   // distance >= 3

    // aCFD fields (Cas9 / 20 bp spacers only)
    bool   acfd_available = false;  // false when spacer length != 20
    size_t n_sdr1_sites   = 0;      // hits with SDR mismatches <= 1
    double aggregate_cfd  = 0.0;    // sum of CFD scores for those hits
    bool   is_promiscuous = false;  // aggregate_cfd > acfd_threshold

    // Genomic feature overlap counts, parallel to the AnnotationBed vector
    // passed to compute_spacer_summary. Empty if no annotations were supplied.
    std::vector<int> annotation_counts;
};

// Named BED annotation for feature overlap counting.
// label appears as n_<label>_overlaps in TSV output.
struct AnnotationBed {
    std::string              label;
    const annot::BedForest*  forest;  // non-owning; caller keeps alive
};

// Configuration for summary computation.
struct SpacerSummaryConfig {
    double acfd_threshold = DEFAULT_ACFD_THRESHOLD;
};

// Compute the summary row for one spacer's hits.
// annotations may be empty to skip all feature overlap counting.
SpacerSummaryRow compute_spacer_summary(
    const std::string&              spacer_name,
    const std::string&              spacer_sequence,
    const std::vector<SearchHit>&   hits,
    const SpacerSummaryConfig&      config,
    const std::vector<AnnotationBed>& annotations = {});

// Format a vector of summary rows as a TSV string (with header).
// annotation_labels must match the order used in compute_spacer_summary.
std::string format_spacer_summary_tsv(const std::vector<SpacerSummaryRow>& rows,
                                      const std::vector<std::string>& annotation_labels = {});

} // namespace spacer
