#include "spacer_summary.hpp"
#include "annotation_bed.hpp"

#include <sstream>

namespace spacer {

int count_sdr_mismatches(const std::string& pattern,
                         const std::string& aligned_seq,
                         const std::string& cigar) {
    // Any indel disqualifies the hit from aCFD (paper: bulge sites inactive).
    for (char c : cigar) {
        if (c == 'I' || c == 'D') return SDR_INDEL;
    }

    // SDR = positions 4-20 (1-indexed) = string indices [3, 20).
    // Position 1 is the 5'/PAM-distal end (index 0).
    const size_t sdr_start_idx = static_cast<size_t>(SDR_START - 1);  // 3
    const size_t sdr_end_idx   = static_cast<size_t>(SDR_END);        // 20 (exclusive)

    const size_t len = std::min(std::min(pattern.size(), aligned_seq.size()), sdr_end_idx);
    if (len <= sdr_start_idx) return 0;

    int mismatches = 0;
    for (size_t i = sdr_start_idx; i < len; ++i) {
        char p = pattern[i];
        char a = aligned_seq[i];
        // N in either sequence = wildcard, no penalty (mirrors CFD scorer).
        if (p == 'N' || a == 'N') continue;
        // Case-insensitive compare (aligned_seq may be lowercase for visual diff).
        if ((p | 0x20) != (a | 0x20)) ++mismatches;
    }
    return mismatches;
}

SpacerSummaryRow compute_spacer_summary(
    const std::string&              spacer_name,
    const std::string&              spacer_sequence,
    const std::vector<SearchHit>&   hits,
    const SpacerSummaryConfig&      config,
    const std::vector<AnnotationBed>& annotations)
{
    SpacerSummaryRow row;
    row.spacer_name     = spacer_name;
    row.spacer_sequence = spacer_sequence;

    // aCFD is only meaningful for 20 bp spacers (Cas9 SDR model).
    row.acfd_available = (spacer_sequence.size() == 20);

    const size_t pattern_len = spacer_sequence.size();

    row.annotation_counts.assign(annotations.size(), 0);

    for (const auto& hit : hits) {
        // Distance buckets
        switch (hit.distance) {
            case 0:  ++row.n_d0;  break;
            case 1:  ++row.n_d1;  break;
            case 2:  ++row.n_d2;  break;
            default: ++row.n_d3p; break;
        }

        // aCFD: only for Cas9-length spacers, only when CFD was computed
        if (row.acfd_available && hit.scoring_info.cfd_score >= 0.0) {
            int sdr_mm = count_sdr_mismatches(spacer_sequence,
                                               hit.mismatch_info.aligned_sequence,
                                               hit.mismatch_info.cigar);
            if (sdr_mm != SDR_INDEL && sdr_mm <= 1) {
                ++row.n_sdr1_sites;
                row.aggregate_cfd += hit.scoring_info.cfd_score;
            }
        }

        // Genomic feature overlaps
        if (!annotations.empty()) {
            size_t start = (hit.chrom_offset >= pattern_len - 1)
                           ? hit.chrom_offset - pattern_len + 1
                           : 0;
            size_t end = hit.chrom_offset + 1;

            for (size_t i = 0; i < annotations.size(); ++i) {
                if (annotations[i].forest &&
                    annot::overlaps_any(hit.chrom_name, start, end, *annotations[i].forest))
                    ++row.annotation_counts[i];
            }
        }
    }

    row.is_promiscuous = row.acfd_available && (row.aggregate_cfd > config.acfd_threshold);

    return row;
}

std::string format_spacer_summary_tsv(const std::vector<SpacerSummaryRow>& rows,
                                      const std::vector<std::string>& annotation_labels) {
    std::ostringstream oss;

    oss << "spacer\tsequence\t"
        << "n_d0\tn_d1\tn_d2\tn_d3p\t"
        << "n_sdr1_sites\taggregate_cfd\tis_promiscuous";
    for (const auto& label : annotation_labels)
        oss << "\tn_" << label << "_overlaps";
    oss << '\n';

    for (const auto& row : rows) {
        oss << row.spacer_name << '\t'
            << row.spacer_sequence << '\t'
            << row.n_d0  << '\t'
            << row.n_d1  << '\t'
            << row.n_d2  << '\t'
            << row.n_d3p << '\t';

        if (row.acfd_available) {
            oss << row.n_sdr1_sites << '\t'
                << row.aggregate_cfd  << '\t'
                << (row.is_promiscuous ? "true" : "false");
        } else {
            oss << ".\t.\t.";
        }

        for (size_t i = 0; i < annotation_labels.size(); ++i) {
            oss << '\t';
            if (i < row.annotation_counts.size())
                oss << row.annotation_counts[i];
            else
                oss << '.';
        }
        oss << '\n';
    }

    return oss.str();
}

} // namespace spacer
