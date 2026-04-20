#include <search_pipeline.hpp>
#include <cfd_scorer.hpp>
#include <shift_add_cpu.hpp>
#include <myers_cpu.hpp>
#include <gpu_engine.hpp>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace {

// Reserve hits.size()/HIT_RESERVE_DIVISOR_GENERAL slots when filtering a
// distances vector — typical hit density is well under 10%.
constexpr size_t HIT_RESERVE_DIVISOR_GENERAL = 10;

// Halo clustering (deduplicate_hits) typically collapses ~3x, so reserve
// hits.size()/DEDUP_EXPECTED_REDUCTION for the output.
constexpr size_t DEDUP_EXPECTED_REDUCTION = 3;

// Below PARALLEL_HITS_MIN_COUNT hits, the thread-launch overhead outweighs the
// per-hit annotation work — stay sequential.
constexpr size_t PARALLEL_HITS_MIN_COUNT = 100;

// Need at least this many hits per worker thread for parallelism to pay off.
constexpr size_t PARALLEL_MIN_CHUNK_SIZE = 50;

// CRISPR seed region: the PAM-proximal 8 nt (positions m-8..m-1 in pattern
// coordinates). Edits here disrupt Cas9 activity more than distal edits.
constexpr size_t SEED_REGION_LENGTH = 8;

// SpCas9 PAM length (NGG / NAG).
constexpr size_t SPCAS9_PAM_LENGTH = 3;

}  // namespace

// Filtering

std::vector<size_t> filter_by_threshold(const std::vector<uint8_t>& distances,
                                        uint8_t threshold) {
    std::vector<size_t> positions;
    positions.reserve(distances.size() / HIT_RESERVE_DIVISOR_GENERAL);

    for (size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] <= threshold) {
            positions.push_back(i);
        }
    }

    return positions;
}

// Hit deduplication

std::vector<SearchHit> deduplicate_hits(std::vector<SearchHit> hits,
                                         size_t pattern_len) {
    if (hits.size() <= 1) return hits;

    // Sort by strand first, then genome position
    std::sort(hits.begin(), hits.end(),
              [](const SearchHit& a, const SearchHit& b) {
                  if (a.strand != b.strand) return a.strand < b.strand;
                  return a.genome_pos < b.genome_pos;
              });

    std::vector<SearchHit> deduped;
    deduped.reserve(hits.size() / DEDUP_EXPECTED_REDUCTION);

    size_t i = 0;
    while (i < hits.size()) {
        // Start a new cluster with hits[i]
        size_t best_idx = i;
        size_t cluster_end = i + 1;

        // Extend cluster: consecutive hits on same strand AND same chromosome
        // that are within (pattern_len - 1) positions of each other.
        // This gap ensures overlapping alignments (halos) are merged, but
        // non-overlapping matches (even if pattern_len apart) are kept separate.
        while (cluster_end < hits.size() &&
               hits[cluster_end].strand == hits[i].strand &&
               hits[cluster_end].chrom_name == hits[cluster_end - 1].chrom_name &&
               hits[cluster_end].genome_pos - hits[cluster_end - 1].genome_pos < pattern_len) {

            // Track the best hit in cluster using PAM-aware selection:
            // 1. Prefer lower distance
            // 2. If tied on distance, prefer better PAM (NGG > NAG > OTHER)
            // 3. If still tied, prefer rightmost position
            bool is_better = false;
            
            if (hits[cluster_end].distance < hits[best_idx].distance) {
                // Lower distance always wins
                is_better = true;
            } else if (hits[cluster_end].distance == hits[best_idx].distance) {
                // Same distance: prefer better PAM
                auto current_pam = hits[cluster_end].mismatch_info.pam_type;
                auto best_pam = hits[best_idx].mismatch_info.pam_type;
                
                // PAM priority: NGG > NAG > OTHER > INCOMPLETE
                auto pam_priority = [](PamType p) {
                    if (p == PamType::NGG) return 3;
                    if (p == PamType::NAG) return 2;
                    if (p == PamType::OTHER) return 1;
                    return 0; // INCOMPLETE
                };
                
                int current_priority = pam_priority(current_pam);
                int best_priority = pam_priority(best_pam);
                
                if (current_priority > best_priority) {
                    // Better PAM
                    is_better = true;
                } else if (current_priority == best_priority) {
                    // Same PAM: prefer rightmost position
                    is_better = true;
                }
            }
            
            if (is_better) {
                best_idx = cluster_end;
            }
            ++cluster_end;
        }

        // Emit only the best hit from this cluster
        deduped.push_back(std::move(hits[best_idx]));
        i = cluster_end;
    }

    return deduped;
}

// Chromosome resolution

std::tuple<std::string, size_t> resolve_chromosome(
    const std::vector<Chromosome>& chromosomes,
    size_t pos) {

    for (const auto& chrom : chromosomes) {
        if (pos >= chrom.start && pos < chrom.start + chrom.length) {
            return {chrom.name, pos - chrom.start};
        }
    }

    throw std::out_of_range("Position " + std::to_string(pos) +
                            " is beyond the genome bounds");
}

// Input normalization

std::string normalize_uracil(const std::string& sequence) {
    std::string result = sequence;
    for (char& c : result) {
        if (c == 'U') c = 'T';
        else if (c == 'u') c = 't';
    }
    return result;
}

// Validation helpers

static void validate_pattern(const std::string& pattern) {
    if (pattern.empty()) {
        throw std::invalid_argument("Pattern cannot be empty");
    }

    for (char c : pattern) {
        switch (c) {
            case 'A': case 'a':
            case 'C': case 'c':
            case 'G': case 'g':
            case 'T': case 't':
            case 'N': case 'n':
                break;
            default:
                throw std::invalid_argument(
                    "Invalid character in pattern: '" + std::string(1, c) + "'");
        }
    }
}

// PAM validation

// Check if a 3-base PAM sequence is valid (NGG or NAG).
// Returns true for valid PAMs, false otherwise.
// Classify a PAM sequence into its type (NGG, NAG, OTHER, or INCOMPLETE).
static PamType classify_pam(const std::string& pam) {
    if (pam.size() < 3) {
        return PamType::INCOMPLETE;
    }

    // N at first position is always valid
    // Second position must be G (for NGG) or A (for NAG)
    // Third position must be G
    char p2 = static_cast<char>(std::toupper(static_cast<unsigned char>(pam[1])));
    char p3 = static_cast<char>(std::toupper(static_cast<unsigned char>(pam[2])));

    // NGG: second and third are G
    if (p2 == 'G' && p3 == 'G') {
        return PamType::NGG;
    }

    // NAG: second is A, third is G
    if (p2 == 'A' && p3 == 'G') {
        return PamType::NAG;
    }

    return PamType::OTHER;
}

// Check if a PAM type matches the specified filter.
// Returns true if the PAM should be included based on the filter.
static bool matches_pam_filter(PamType pam_type, PamFilter filter) {
    if (filter == PamFilter::NONE || filter == PamFilter::ANY) {
        return true;  // No filtering
    }

    switch (filter) {
        case PamFilter::NGG_ONLY:
            return pam_type == PamType::NGG;
        case PamFilter::NAG_ONLY:
            return pam_type == PamType::NAG;
        case PamFilter::NGG_OR_NAG:
            return pam_type == PamType::NGG || pam_type == PamType::NAG;
        default:
            return true;
    }
}

// Parallel annotation helper: work-stealing over hits for post-search enrichment.

// Work-stealing parallel map over a vector of hits. Stays sequential unless
// hit_count and per-thread workload both clear the overhead thresholds.
template<typename Func>
static void parallel_process_hits(std::vector<SearchHit>& hits,
                                  Func func,
                                  size_t num_threads) {
    const size_t hit_count = hits.size();

    if (hit_count < PARALLEL_HITS_MIN_COUNT || num_threads <= 1) {
        for (auto& hit : hits) {
            func(hit);
        }
        return;
    }

    size_t effective_threads = std::min(num_threads, hit_count / PARALLEL_MIN_CHUNK_SIZE);
    if (effective_threads <= 1) {
        // Not enough work to justify threading
        for (auto& hit : hits) {
            func(hit);
        }
        return;
    }

    // Parallel path: Work-stealing pattern
    std::atomic<size_t> next_index{0};

    auto worker = [&]() {
        while (true) {
            size_t i = next_index.fetch_add(1);
            if (i >= hit_count) break;
            func(hits[i]);
        }
    };

    // Launch threads
    std::vector<std::thread> threads;
    threads.reserve(effective_threads);
    for (size_t t = 0; t < effective_threads; ++t) {
        threads.emplace_back(worker);
    }

    // Wait for completion
    for (auto& t : threads) {
        t.join();
    }
}

// Lightweight PAM-only extraction (runs before full annotation so PAM filters
// can reject hits before the expensive alignment work).

// Extract only the PAM sequence and classify it, without doing full alignment.
// This is the cheap part of extract_mismatch_info() - just genome slicing + classification.
// Used for early PAM filtering to avoid expensive annotation for hits that will be filtered out.
static void extract_pam_only(const GenomeView& view, SearchHit& hit, size_t pattern_len) {
    // For SpCas9, PAM (NGG) is 3' of the protospacer
    if (hit.strand == Strand::PLUS) {
        size_t pam_start = hit.genome_pos + 1;
        if (pam_start + SPCAS9_PAM_LENGTH <= view.total_bases) {
            hit.mismatch_info.pam_sequence = extract_genome_slice(view, pam_start, SPCAS9_PAM_LENGTH);
        } else if (pam_start < view.total_bases) {
            hit.mismatch_info.pam_sequence = extract_genome_slice(view, pam_start,
                                                                   view.total_bases - pam_start);
        } else {
            hit.mismatch_info.pam_sequence = "";
        }
    } else {
        // Minus strand: PAM is 5' of the RC alignment on the forward strand
        size_t align_begin = (hit.genome_pos >= pattern_len - 1) ? hit.genome_pos - pattern_len + 1 : 0;
        if (align_begin >= SPCAS9_PAM_LENGTH) {
            size_t pam_start = align_begin - SPCAS9_PAM_LENGTH;
            hit.mismatch_info.pam_sequence = reverse_complement(
                extract_genome_slice(view, pam_start, SPCAS9_PAM_LENGTH));
        } else if (align_begin > 0) {
            hit.mismatch_info.pam_sequence = reverse_complement(
                extract_genome_slice(view, 0, align_begin));
        } else {
            hit.mismatch_info.pam_sequence = "";
        }
    }

    hit.mismatch_info.pam_type = classify_pam(hit.mismatch_info.pam_sequence);
}

// Alignment computation

// Compute semi-global alignment between pattern and text using standard DP.
// Returns a tuple: alignment operations, ambiguity flag, and the text end-index
// the traceback was anchored at (so the caller can derive align_text_start).
// When multiple operations lead to the same edit distance at a cell,
// we use tie-breaking priority: MATCH > DNA_BULGE (INS_TEXT) > RNA_BULGE (INS_PATTERN)
// The ambiguity flag is set to true if alternative equally-optimal alignments exist.
//
// `anchor_right` selects where the alignment must end:
//   true  — alignment ends at text position n (plus-strand windows where
//           the caller placed the Myers end-position at the last text char)
//   false — alignment ends at the leftmost j minimizing dp[m][j]
//           (minus-strand windows where the caller RC'd the window, putting
//           the true alignment at the START of text with free end gaps)
struct AlignmentResult {
    std::vector<AlignOp> alignment;
    bool   has_ambiguity = false;
    size_t n_ambiguous_cells = 0;
    size_t best_j = 0;
};

static AlignmentResult compute_alignment(
        const std::string& pattern,
        const std::string& text,
        bool anchor_right) {
    size_t m = pattern.size();
    size_t n = text.size();
    
    // DP matrix: dp[i][j] = edit distance for pattern[0..i) vs text[0..j).
    // Border convention depends on strand semantics (see anchor_right):
    //   anchor_right == true  (plus strand):  free start gaps in text
    //                                          → dp[0][j] = 0
    //                                          → alignment lives at a suffix of text
    //   anchor_right == false (minus strand): no free start gaps in text
    //                                          → dp[0][j] = j
    //                                          → alignment must begin at text[0]
    //                                            (which is the Myers end-position after RC)
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));

    for (size_t i = 0; i <= m; ++i) {
        dp[i][0] = static_cast<int>(i);
    }

    for (size_t j = 0; j <= n; ++j) {
        dp[0][j] = anchor_right ? 0 : static_cast<int>(j);
    }
    
    // Fill DP table
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pattern[i-1])));
            char t = static_cast<char>(std::toupper(static_cast<unsigned char>(text[j-1])));
            
            // Match/mismatch cost
            // Pattern N = wildcard (matches any text base, cost 0)
            // Genome N = masked/unknown (matches nothing, cost 1)
            // This matches Myers' algorithm semantics where genome N produces Eq=0
            int match_cost;
            if (p == 'N') {
                match_cost = 0;  // Pattern N is wildcard
            } else if (t == 'N') {
                match_cost = 1;  // Genome N is masked (unknown base)
            } else {
                match_cost = (p == t) ? 0 : 1;
            }

            // Three options: match/mismatch, insert in text, insert in pattern
            int match = dp[i-1][j-1] + match_cost;
            int ins_text = dp[i][j-1] + 1;      // horizontal = gap in pattern
            int ins_pattern = dp[i-1][j] + 1;   // vertical = gap in text
            
            dp[i][j] = std::min({match, ins_text, ins_pattern});
        }
    }
    
    // Choose the traceback end-column based on strand semantics (see header).
    size_t best_j;
    if (anchor_right) {
        best_j = n;
    } else {
        int best_score = dp[m][n];
        best_j = n;
        for (size_t jj = 1; jj <= n; ++jj) {
            if (dp[m][jj] < best_score) {
                best_score = dp[m][jj];
                best_j = jj;
            }
        }
    }

    std::vector<AlignOp> alignment;
    alignment.reserve(m + best_j);

    bool has_ambiguity = false;  // Track if alternative alignments exist
    size_t n_ambiguous_cells = 0;

    size_t i = m;
    size_t j = best_j;
    
    while (i > 0 && j > 0) {
        char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pattern[i-1])));
        char t = static_cast<char>(std::toupper(static_cast<unsigned char>(text[j-1])));
        
        int match_cost;
        if (p == 'N') {
            match_cost = 0;
        } else if (t == 'N') {
            match_cost = 1;
        } else {
            match_cost = (p == t) ? 0 : 1;
        }
        int current = dp[i][j];
        
        // Count how many operations lead to the current cell with optimal cost
        int num_optimal = 0;
        bool can_match = false, can_ins_text = false, can_ins_pattern = false;
        
        if (i > 0 && j > 0 && current == dp[i-1][j-1] + match_cost) {
            can_match = true;
            ++num_optimal;
        }
        if (j > 0 && current == dp[i][j-1] + 1) {
            can_ins_text = true;
            ++num_optimal;
        }
        if (i > 0 && current == dp[i-1][j] + 1) {
            can_ins_pattern = true;
            ++num_optimal;
        }
        
        // If multiple operations are optimal, flag ambiguity
        if (num_optimal > 1) {
            has_ambiguity = true;
            ++n_ambiguous_cells;
        }
        
        // Apply tie-breaking priority: MATCH > INS_TEXT (DNA bulge) > INS_PATTERN (RNA bulge)
        if (can_match) {
            alignment.push_back(AlignOp::MATCH);
            --i;
            --j;
        } else if (can_ins_text) {
            alignment.push_back(AlignOp::INS_TEXT);
            --j;
        } else if (can_ins_pattern) {
            alignment.push_back(AlignOp::INS_PATTERN);
            --i;
        } else {
            // Should not happen with correct DP
            break;
        }
    }
    
    // Handle remaining pattern positions (if any)
    while (i > 0) {
        alignment.push_back(AlignOp::INS_PATTERN);
        --i;
    }
    
    // Alignment was built backwards, reverse it
    std::reverse(alignment.begin(), alignment.end());

    return {std::move(alignment), has_ambiguity, n_ambiguous_cells, best_j};
}

// Build a run-length CIGAR string from an alignment op vector.
// Uses M for match-or-mismatch (like SAM CIGAR without =/X), I for RNA bulge
// (pattern consumes a base that genome doesn't), D for DNA bulge (genome
// consumes a base the pattern doesn't). Empty alignment -> empty string.
static std::string format_cigar(const std::vector<AlignOp>& alignment) {
    std::string out;
    if (alignment.empty()) return out;

    auto op_char = [](AlignOp op) -> char {
        switch (op) {
            case AlignOp::MATCH:       return 'M';
            case AlignOp::INS_PATTERN: return 'I';  // gap in text (RNA bulge)
            case AlignOp::INS_TEXT:    return 'D';  // gap in pattern (DNA bulge)
        }
        return '?';
    };

    char current = op_char(alignment[0]);
    size_t run = 1;
    for (size_t k = 1; k < alignment.size(); ++k) {
        char c = op_char(alignment[k]);
        if (c == current) {
            ++run;
        } else {
            out += std::to_string(run);
            out += current;
            current = c;
            run = 1;
        }
    }
    out += std::to_string(run);
    out += current;
    return out;
}

// Mismatch extraction

// Extract detailed mismatch information for a hit.
// For minus strand hits, the aligned_sequence is stored in forward orientation
// but comparison and mismatch positions are computed in the pattern's orientation.
// Extract the 3' SpCas9 PAM (3 nt, NGG) adjacent to the protospacer.
//   Plus strand:  5'-[spacer]-PAM-3'   → PAM at genome_pos + 1
//   Minus strand: PAM is 5' of the RC alignment on the forward strand; extract
//                 from the forward strand and RC to get the PAM on the minus strand.
// Writes hit.mismatch_info.pam_sequence and pam_type.
static void extract_pam_sequence(GenomeView view, SearchHit& hit, size_t m) {
    if (hit.strand == Strand::PLUS) {
        size_t pam_start = hit.genome_pos + 1;
        if (pam_start + SPCAS9_PAM_LENGTH <= view.total_bases) {
            hit.mismatch_info.pam_sequence = extract_genome_slice(view, pam_start, SPCAS9_PAM_LENGTH);
        } else if (pam_start < view.total_bases) {
            hit.mismatch_info.pam_sequence = extract_genome_slice(view, pam_start,
                                                                   view.total_bases - pam_start);
        } else {
            hit.mismatch_info.pam_sequence = "";
        }
    } else {
        size_t align_begin = (hit.genome_pos >= m - 1) ? hit.genome_pos - m + 1 : 0;
        if (align_begin >= SPCAS9_PAM_LENGTH) {
            hit.mismatch_info.pam_sequence = reverse_complement(
                extract_genome_slice(view, align_begin - SPCAS9_PAM_LENGTH, SPCAS9_PAM_LENGTH));
        } else if (align_begin > 0) {
            hit.mismatch_info.pam_sequence = reverse_complement(
                extract_genome_slice(view, 0, align_begin));
        } else {
            hit.mismatch_info.pam_sequence = "";
        }
    }
    hit.mismatch_info.pam_type = classify_pam(hit.mismatch_info.pam_sequence);
}

// Distance-zero case: the aligned m-base window has no gaps, so we take the
// last (PLUS) or first (MINUS, RC'd) m bases of the extracted compare_sequence.
static std::string extract_aligned_region_perfect_match(const std::string& compare_sequence,
                                                        Strand strand, size_t m) {
    if (compare_sequence.size() < m) return compare_sequence;
    return (strand == Strand::PLUS)
           ? compare_sequence.substr(compare_sequence.size() - m)
           : compare_sequence.substr(0, m);
}

// Walk alignment ops, emitting the text base for each (MATCH, INS_TEXT) op and
// '-' for RNA bulges (INS_PATTERN). Produces the genome-side aligned string.
static std::string build_aligned_sequence_from_alignment(const std::vector<AlignOp>& alignment,
                                                         const std::string& compare_sequence,
                                                         size_t align_text_start) {
    std::string out;
    out.reserve(alignment.size());
    size_t seq_idx = align_text_start;
    for (const AlignOp& op : alignment) {
        if (op == AlignOp::MATCH || op == AlignOp::INS_TEXT) {
            if (seq_idx < compare_sequence.size()) out += compare_sequence[seq_idx];
            ++seq_idx;
        } else if (op == AlignOp::INS_PATTERN) {
            out += '-';
        }
    }
    return out;
}

// Split edits into PAM-proximal seed vs distal (positions 1..m-SEED_REGION_LENGTH).
// Pattern 'N' is a wildcard (no mismatch recorded); genome 'N' is masked (always
// a mismatch). Positions reported are 1-based in pattern coordinates.
static void classify_edits_by_region(const std::vector<AlignOp>& alignment,
                                     const std::string& pattern,
                                     const std::string& text,
                                     size_t align_text_start,
                                     size_t m,
                                     MismatchInfo& out) {
    const size_t distal_end = (m > SEED_REGION_LENGTH) ? m - SEED_REGION_LENGTH : 0;
    size_t pattern_pos = 0;
    size_t text_pos = align_text_start;

    auto record_edit = [&](EditType type) {
        out.mismatch_positions.push_back(static_cast<uint8_t>(pattern_pos + 1));
        out.edit_types.push_back(type);
        const bool distal = pattern_pos < distal_end;
        const bool in_seed = !distal && pattern_pos < m;
        switch (type) {
            case EditType::MISMATCH:
                if (distal)       ++out.distal_mismatches;
                else if (in_seed) ++out.seed_mismatches;
                break;
            case EditType::DNA_BULGE:
                if (distal)       ++out.distal_dna_bulges;
                else if (in_seed) ++out.seed_dna_bulges;
                break;
            case EditType::RNA_BULGE:
                if (distal)       ++out.distal_rna_bulges;
                else if (in_seed) ++out.seed_rna_bulges;
                break;
            case EditType::NONE:
                break;
        }
    };

    for (const AlignOp& op : alignment) {
        switch (op) {
            case AlignOp::MATCH: {
                if (pattern_pos < pattern.size() && text_pos < text.size()) {
                    char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pattern[pattern_pos])));
                    char t = static_cast<char>(std::toupper(static_cast<unsigned char>(text[text_pos])));
                    if (p != 'N' && (t == 'N' || p != t)) {
                        record_edit(EditType::MISMATCH);
                    }
                }
                ++pattern_pos;
                ++text_pos;
                break;
            }
            case AlignOp::INS_TEXT:
                record_edit(EditType::DNA_BULGE);
                ++text_pos;
                break;
            case AlignOp::INS_PATTERN:
                record_edit(EditType::RNA_BULGE);
                ++pattern_pos;
                break;
        }
    }
}

static void extract_mismatch_info(const std::string& original_pattern,
                                  GenomeView view,
                                  SearchHit& hit) {
    const size_t m = original_pattern.size();

    // Reset per-hit mutable fields (pam_sequence may already be set by the
    // early lazy-PAM path — preserved across the clear).
    hit.mismatch_info.mismatch_positions.clear();
    hit.mismatch_info.edit_types.clear();
    hit.mismatch_info.seed_mismatches = 0;
    hit.mismatch_info.seed_dna_bulges = 0;
    hit.mismatch_info.seed_rna_bulges = 0;
    hit.mismatch_info.distal_mismatches = 0;
    hit.mismatch_info.distal_dna_bulges = 0;
    hit.mismatch_info.distal_rna_bulges = 0;
    hit.mismatch_info.alignment_is_ambiguous = false;

    // Extract a window big enough to hold the alignment (up to m + distance bases
    // before the Myers end-position, to account for pattern-side insertions).
    const size_t max_align_len = m + hit.distance + 1;
    const size_t align_start = (hit.genome_pos >= max_align_len - 1)
                               ? hit.genome_pos - max_align_len + 1 : 0;
    const size_t align_len   = hit.genome_pos - align_start + 1;
    const std::string genome_seq = extract_genome_slice(view, align_start, align_len);

    if (hit.mismatch_info.pam_sequence.empty()) {
        extract_pam_sequence(view, hit, m);
    }

    // For the MINUS strand we searched with the RC pattern, so RC the window
    // before comparison. The PAM was already RC'd by extract_pam_sequence.
    const std::string& compare_pattern = original_pattern;
    std::string compare_sequence = (hit.strand == Strand::MINUS)
        ? reverse_complement(genome_seq) : genome_seq;

    if (hit.distance == 0) {
        hit.mismatch_info.aligned_sequence =
            extract_aligned_region_perfect_match(compare_sequence, hit.strand, m);
        hit.mismatch_info.cigar = std::to_string(m) + "M";
        return;
    }

    // distance > 0: run banded alignment. Plus strand anchors at j=n (Myers end),
    // minus strand (RC'd window) anchors at leftmost min of last row.
    const bool anchor_right = (hit.strand == Strand::PLUS);
    AlignmentResult ar = compute_alignment(compare_pattern, compare_sequence, anchor_right);
    hit.mismatch_info.alignment_is_ambiguous = ar.has_ambiguity;
    hit.mismatch_info.n_ambiguous_cells = static_cast<uint16_t>(
        std::min<size_t>(ar.n_ambiguous_cells, std::numeric_limits<uint16_t>::max()));
    hit.mismatch_info.cigar = format_cigar(ar.alignment);

    size_t text_consumed = 0;
    for (const AlignOp& op : ar.alignment) {
        if (op == AlignOp::MATCH || op == AlignOp::INS_TEXT) ++text_consumed;
    }
    const size_t align_text_start = ar.best_j - text_consumed;

    hit.mismatch_info.aligned_sequence = build_aligned_sequence_from_alignment(
        ar.alignment, compare_sequence, align_text_start);

    classify_edits_by_region(ar.alignment, compare_pattern, compare_sequence,
                             align_text_start, m, hit.mismatch_info);
}

// Search pipeline

// Assemble a SearchHit from absolute genome coords + distance + strand.
// Populates chrom_name / chrom_offset by looking up the chromosome that
// contains abs_pos. Downstream annotation (PAM, mismatch info, scores) is
// layered on top later.
static SearchHit make_search_hit(const GenomeView& view, size_t abs_pos,
                                 uint8_t distance, Strand strand) {
    auto [chrom_name, chrom_offset] = resolve_chromosome(*view.chromosomes, abs_pos);
    SearchHit hit;
    hit.genome_pos   = abs_pos;
    hit.distance     = distance;
    hit.chrom_name   = std::move(chrom_name);
    hit.chrom_offset = chrom_offset;
    hit.strand       = strand;
    return hit;
}

SearchResult search_genome(const SearchConfig& config, GenomeView view) {
    // Validate input
    validate_pattern(config.pattern);

    SearchResult result;
    result.pattern = config.pattern;
    result.threshold = config.threshold;
    result.total_positions = view.total_bases;
    result.used_gpu = false;

    // Prepare patterns for both strands
    std::string fwd_pattern = config.pattern;
    std::string rc_pattern = reverse_complement(config.pattern);

    // Decide whether to use GPU
    bool gpu_size_ok = config.pattern.size() <= 64;
    bool use_gpu = config.prefer_gpu &&
                   gpu_available() &&
                   gpu_size_ok;  // GPU kernel limit

    // Inform the user when GPU was requested but the pattern is too long
    if (config.prefer_gpu && gpu_available() && !gpu_size_ok) {
        std::cerr << "Note: pattern length " << config.pattern.size()
                  << " bp exceeds 64-bp GPU kernel limit; falling back to CPU multi-word Myers.\n";
    }

    // Resolve search window. search_end == 0 means "whole genome".
    size_t search_start = config.search_start;
    size_t search_end = (config.search_end == 0) ? view.total_bases : config.search_end;
    if (search_end > view.total_bases) search_end = view.total_bases;
    if (search_start > search_end) search_start = search_end;
    size_t search_length = search_end - search_start;

    // Lambda to run search on one strand
    auto search_strand = [&](const std::string& pattern, Strand strand)
        -> std::vector<SearchHit>
    {
        std::vector<uint8_t> distances;

        // Dispatch based on distance mode. Kernels take (start, length) and
        // return distances indexed relative to start.
        if (config.distance_mode == DistanceMode::HAMMING) {
            // Hamming distance (substitution-only)
            if (use_gpu) {
                // GPU shift-add implementation
                auto gpu_result = shift_add_gpu(pattern, view, search_start, search_length);
                distances = std::move(gpu_result.distances);
                result.used_gpu = true;
            } else {
                auto cpu_result = shift_add_hamming_distance_bv(pattern, view, search_start, search_length);
                distances = std::move(cpu_result.distances);
            }
        } else {
            // Levenshtein distance (default: substitutions + insertions + deletions)
            if (use_gpu) {
                auto gpu_result = myers_gpu(pattern, view, search_start, search_length);
                distances = std::move(gpu_result.distances);
                result.used_gpu = true;
            } else {
                auto cpu_result = myers_edit_distances_bv(pattern, view, search_start, search_length);
                distances = std::move(cpu_result.distances);
            }
        }

        // Filter by threshold (returned indices are relative to search_start)
        auto hit_positions = filter_by_threshold(distances, config.threshold);

        // Build SearchHit objects with chromosome info
        std::vector<SearchHit> hits;
        hits.reserve(hit_positions.size());

        for (size_t rel_pos : hit_positions) {
            size_t abs_pos = search_start + rel_pos;
            hits.push_back(make_search_hit(view, abs_pos, distances[rel_pos], strand));
        }

        return hits;
    };

    // Strand-pass gating:
    //   forward_only → run only forward pass (+ hits)
    //   reverse_only → run only RC pass (- hits)
    //   neither     → both passes (default; gated by search_both_strands)
    const bool do_fwd = !config.reverse_only;
    const bool do_rc  = config.reverse_only ||
                        (config.search_both_strands && !config.forward_only);

    if (do_fwd) {
        result.hits = search_strand(fwd_pattern, Strand::PLUS);
    }
    if (do_rc) {
        auto rc_hits = search_strand(rc_pattern, Strand::MINUS);
        result.hits.insert(result.hits.end(),
                           std::make_move_iterator(rc_hits.begin()),
                           std::make_move_iterator(rc_hits.end()));
    }

    // Sort by genome position and strand
    std::sort(result.hits.begin(), result.hits.end());

    // Thread count for parallel annotation.
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 1;

    // Lazy PAM extraction: extract PAM first (cheap), filter, then run the
    // full alignment annotation only for hits that pass the PAM filter.
    bool enable_lazy_pam = (config.pam_filter != PamFilter::NONE &&
                           config.pam_filter != PamFilter::ANY &&
                           config.compute_mismatches);

    if (enable_lazy_pam) {
        // Step 1: Extract PAM only for all hits (cheap - just genome slicing) - PARALLEL
        parallel_process_hits(result.hits,
            [&](SearchHit& hit) {
                extract_pam_only(view, hit, config.pattern.size());
            },
            num_threads);

        // Step 2: Filter by PAM early (discards ~93% of hits for NGG filter)
        std::vector<SearchHit> filtered_hits;
        filtered_hits.reserve(result.hits.size() / 10);  // Expect ~10% to pass

        for (auto& hit : result.hits) {
            if (matches_pam_filter(hit.mismatch_info.pam_type, config.pam_filter)) {
                filtered_hits.push_back(std::move(hit));
            }
        }

        result.hits = std::move(filtered_hits);

        // Step 3: Extract full mismatch info only for passing hits (expensive but only 6-7% of original set) - PARALLEL
        parallel_process_hits(result.hits,
            [&](SearchHit& hit) {
                extract_mismatch_info(config.pattern, view, hit);
            },
            num_threads);
    } else {
        // Original path: Extract full mismatch info for all hits - PARALLEL
        if (config.compute_mismatches) {
            parallel_process_hits(result.hits,
                [&](SearchHit& hit) {
                    extract_mismatch_info(config.pattern, view, hit);
                },
                num_threads);
        }
    }

    // Deduplicate overlapping semi-global alignment hits.
    // Myers' algorithm reports edit distance at every position, creating
    // "halos" of nearby low-distance positions around each real match.
    // Deduplication is PAM-aware: when clustering hits with equal distance,
    // it prefers hits with valid PAMs (NGG > NAG > OTHER).
    if (!config.disable_deduplication) {
        result.hits = deduplicate_hits(std::move(result.hits), config.pattern.size());
    }

    // Compute CFD scores if requested (requires mismatch info) - PARALLEL
    if (config.compute_scores && config.compute_mismatches) {
        parallel_process_hits(result.hits,
            [&](SearchHit& hit) {
                hit.scoring_info.cfd_score = cfd::compute_cfd_score(
                    config.pattern, hit.mismatch_info);
                hit.scoring_info.risk_tier = cfd::classify_risk_tier(
                    hit.scoring_info.cfd_score);
            },
            num_threads);
    }

    // PAM filtering already done above in lazy path

    // Apply max_hits limit if specified
    if (config.max_hits > 0 && result.hits.size() > config.max_hits) {
        result.hits.resize(config.max_hits);
    }

    return result;
}

// Output formatting

static const char* pam_type_str(PamType t) {
    switch (t) {
        case PamType::NGG:        return "NGG";
        case PamType::NAG:        return "NAG";
        case PamType::OTHER:      return "OTHER";
        case PamType::INCOMPLETE: return "INCOMPLETE";
    }
    return "OTHER";
}

static const char* edit_type_str(EditType t) {
    switch (t) {
        case EditType::NONE:      return "NONE";
        case EditType::MISMATCH:  return "MISMATCH";
        case EditType::DNA_BULGE: return "DNA_BULGE";
        case EditType::RNA_BULGE: return "RNA_BULGE";
    }
    return "NONE";
}

static const char* risk_tier_str(uint8_t tier) {
    switch (tier) {
        case 0: return "LOW";
        case 1: return "MEDIUM";
        case 2: return "HIGH";
    }
    return "UNKNOWN";
}

static const char* specificity_tier_str(uint8_t tier) {
    switch (tier) {
        case 0: return "POOR";
        case 1: return "FAIR";
        case 2: return "GOOD";
        case 3: return "EXCELLENT";
    }
    return "UNKNOWN";
}

static std::string format_mismatch_positions(const std::vector<uint8_t>& positions) {
    if (positions.empty()) return ".";
    std::string out;
    for (size_t i = 0; i < positions.size(); ++i) {
        if (i > 0) out += ",";
        out += std::to_string(positions[i]);
    }
    return out;
}

static std::string format_edit_types(const std::vector<EditType>& edits) {
    if (edits.empty()) return ".";
    std::string out;
    for (size_t i = 0; i < edits.size(); ++i) {
        if (i > 0) out += ",";
        out += edit_type_str(edits[i]);
    }
    return out;
}

std::string format_hits_tsv(const std::vector<SearchHit>& hits,
                            const std::string& pattern) {
    std::ostringstream oss;

    // Header with new biologically-meaningful columns
    oss << "chrom\tstart\tend\tpattern\tdistance\tstrand\t"
        << "aligned_seq\tmismatch_pos\tedit_types\tcigar\t"
        << "pam_seq\tpam_type\t"
        << "seed_edits\tdistal_edits\t"
        << "seed_mismatches\tseed_dna_bulges\tseed_rna_bulges\t"
        << "distal_mismatches\tdistal_dna_bulges\tdistal_rna_bulges\t"
        << "alignment_ambiguous\tn_ambiguous_cells\t"
        << "cfd_score\trisk_tier\n";

    // Data rows
    size_t pattern_len = pattern.size();
    for (const auto& hit : hits) {
        // Start position is where the alignment ends minus pattern length + 1
        // For semi-global alignment ending at hit.chrom_offset:
        // start = chrom_offset - pattern_len + 1 (but can't be negative)
        size_t start = (hit.chrom_offset >= pattern_len - 1)
                       ? hit.chrom_offset - pattern_len + 1
                       : 0;
        size_t end = hit.chrom_offset + 1;  // exclusive end

        oss << hit.chrom_name << '\t'
            << start << '\t'
            << end << '\t'
            << pattern << '\t'
            << static_cast<int>(hit.distance) << '\t'
            << static_cast<char>(hit.strand) << '\t'
            << (hit.mismatch_info.aligned_sequence.empty() ? "." : hit.mismatch_info.aligned_sequence) << '\t'
            << format_mismatch_positions(hit.mismatch_info.mismatch_positions) << '\t'
            << format_edit_types(hit.mismatch_info.edit_types) << '\t'
            << (hit.mismatch_info.cigar.empty() ? "." : hit.mismatch_info.cigar) << '\t'
            << (hit.mismatch_info.pam_sequence.empty() ? "." : hit.mismatch_info.pam_sequence) << '\t'
            << pam_type_str(hit.mismatch_info.pam_type) << '\t'
            << static_cast<int>(hit.mismatch_info.total_seed_edits()) << '\t'
            << static_cast<int>(hit.mismatch_info.total_distal_edits()) << '\t'
            << static_cast<int>(hit.mismatch_info.seed_mismatches) << '\t'
            << static_cast<int>(hit.mismatch_info.seed_dna_bulges) << '\t'
            << static_cast<int>(hit.mismatch_info.seed_rna_bulges) << '\t'
            << static_cast<int>(hit.mismatch_info.distal_mismatches) << '\t'
            << static_cast<int>(hit.mismatch_info.distal_dna_bulges) << '\t'
            << static_cast<int>(hit.mismatch_info.distal_rna_bulges) << '\t'
            << (hit.mismatch_info.alignment_is_ambiguous ? "true" : "false") << '\t'
            << static_cast<int>(hit.mismatch_info.n_ambiguous_cells) << '\t'
            << hit.scoring_info.cfd_score << '\t'
            << risk_tier_str(hit.scoring_info.risk_tier) << '\n';
    }

    return oss.str();
}

std::string format_hits_bed(const std::vector<SearchHit>& hits,
                            const std::string& pattern) {
    std::ostringstream oss;

    size_t pattern_len = pattern.size();
    for (const auto& hit : hits) {
        size_t start = (hit.chrom_offset >= pattern_len - 1)
                       ? hit.chrom_offset - pattern_len + 1
                       : 0;
        size_t end = hit.chrom_offset + 1;

        // BED score: higher is better, so use (pattern_len - distance)
        int score = static_cast<int>(pattern_len) - static_cast<int>(hit.distance);

        oss << hit.chrom_name << '\t'
            << start << '\t'
            << end << '\t'
            << pattern << '\t'
            << score << '\t'
            << static_cast<char>(hit.strand) << '\n';
    }

    return oss.str();
}

// JSON output helpers

namespace {

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b";  break;
            case '\f': out += "\\f";  break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x",
                                  static_cast<unsigned char>(c));
                    out += buf;
                } else {
                    out += c;
                }
        }
    }
    return out;
}

// Emit the body of one hit as JSON key/value pairs (no enclosing braces).
// Caller is responsible for "{" ... "}" and any leading fields like spacer name.
void write_hit_body_json(std::ostringstream& oss,
                          const SearchHit& hit,
                          const std::string& pattern) {
    const size_t pattern_len = pattern.size();
    const size_t start = (hit.chrom_offset >= pattern_len - 1)
                        ? hit.chrom_offset - pattern_len + 1
                        : 0;
    const size_t end = hit.chrom_offset + 1;

    oss << "\"chrom\":\""    << json_escape(hit.chrom_name) << "\","
        << "\"start\":"      << start << ","
        << "\"end\":"        << end << ","
        << "\"pattern\":\""  << json_escape(pattern) << "\","
        << "\"distance\":"   << static_cast<int>(hit.distance) << ","
        << "\"strand\":\""   << static_cast<char>(hit.strand) << "\",";

    const auto& mi = hit.mismatch_info;

    oss << "\"aligned_seq\":";
    if (mi.aligned_sequence.empty()) oss << "null";
    else oss << "\"" << json_escape(mi.aligned_sequence) << "\"";
    oss << ",";

    oss << "\"mismatch_pos\":[";
    for (size_t i = 0; i < mi.mismatch_positions.size(); ++i) {
        if (i) oss << ",";
        oss << static_cast<int>(mi.mismatch_positions[i]);
    }
    oss << "],";

    oss << "\"edit_types\":[";
    for (size_t i = 0; i < mi.edit_types.size(); ++i) {
        if (i) oss << ",";
        oss << "\"" << edit_type_str(mi.edit_types[i]) << "\"";
    }
    oss << "],";

    oss << "\"cigar\":";
    if (mi.cigar.empty()) oss << "null";
    else oss << "\"" << json_escape(mi.cigar) << "\"";
    oss << ",";

    oss << "\"pam_seq\":";
    if (mi.pam_sequence.empty()) oss << "null";
    else oss << "\"" << json_escape(mi.pam_sequence) << "\"";
    oss << ",";

    oss << "\"pam_type\":\""           << pam_type_str(mi.pam_type) << "\","
        << "\"seed_edits\":"           << static_cast<int>(mi.total_seed_edits()) << ","
        << "\"distal_edits\":"         << static_cast<int>(mi.total_distal_edits()) << ","
        << "\"seed_mismatches\":"      << static_cast<int>(mi.seed_mismatches) << ","
        << "\"seed_dna_bulges\":"      << static_cast<int>(mi.seed_dna_bulges) << ","
        << "\"seed_rna_bulges\":"      << static_cast<int>(mi.seed_rna_bulges) << ","
        << "\"distal_mismatches\":"    << static_cast<int>(mi.distal_mismatches) << ","
        << "\"distal_dna_bulges\":"    << static_cast<int>(mi.distal_dna_bulges) << ","
        << "\"distal_rna_bulges\":"    << static_cast<int>(mi.distal_rna_bulges) << ","
        << "\"alignment_ambiguous\":"  << (mi.alignment_is_ambiguous ? "true" : "false") << ","
        << "\"n_ambiguous_cells\":"    << static_cast<int>(mi.n_ambiguous_cells) << ","
        << "\"cfd_score\":"            << hit.scoring_info.cfd_score << ","
        << "\"risk_tier\":\""          << risk_tier_str(hit.scoring_info.risk_tier) << "\"";
}

}  // anonymous namespace

std::string format_hits_json(const std::vector<SearchHit>& hits,
                              const std::string& pattern) {
    std::ostringstream oss;
    oss.precision(6);
    oss << std::fixed;
    oss << "[";
    for (size_t i = 0; i < hits.size(); ++i) {
        if (i) oss << ",";
        oss << "\n  {";
        write_hit_body_json(oss, hits[i], pattern);
        oss << "}";
    }
    if (!hits.empty()) oss << "\n";
    oss << "]\n";
    return oss.str();
}

// Batch processing

std::vector<SpacerEntry> parse_spacer_file(const std::string& filepath) {
    // "-" means read from stdin (CLI convention).
    bool use_stdin = (filepath == "-");
    std::ifstream file;
    if (!use_stdin) {
        file.open(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open spacer file: " + filepath);
        }
    }
    std::istream& in = use_stdin ? std::cin : static_cast<std::istream&>(file);
    const std::string display_path = use_stdin ? "<stdin>" : filepath;

    std::vector<SpacerEntry> spacers;
    std::string line;
    size_t line_num = 0;
    size_t unnamed_count = 0;

    while (std::getline(in, line)) {
        ++line_num;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // Trim leading whitespace
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;  // Whitespace-only line

        // Parse tab/whitespace-delimited fields
        std::istringstream iss(line.substr(start));
        std::string field1, field2;

        if (!(iss >> field1)) continue;  // Empty after trim

        SpacerEntry entry;
        if (iss >> field2) {
            // Two or more fields: name<TAB>sequence[<TAB>...]
            entry.name = field1;
            entry.sequence = field2;
        } else {
            // Single field: just sequence
            entry.name = "spacer_" + std::to_string(++unnamed_count);
            entry.sequence = field1;
        }
        entry.source_file = display_path;
        entry.source_line = line_num;

        // Note: Validation deferred to search time to allow graceful skipping
        // of individual invalid spacers without aborting the entire batch

        spacers.push_back(std::move(entry));
    }

    if (spacers.empty()) {
        throw std::runtime_error(
            "No valid spacers found in " + display_path +
            " (no non-comment, non-blank lines seen)");
    }

    return spacers;
}

// Compute MIT specificity scores for all spacers in a batch result.
// This aggregates CFD scores across all off-targets for each spacer.
static void compute_batch_mit_scores(BatchSearchResult& batch_result) {
    for (auto& spacer_result : batch_result.spacer_results) {
        // Collect CFD scores from all hits for this spacer
        std::vector<double> cfd_scores;
        cfd_scores.reserve(spacer_result.result.hits.size());
        
        for (const auto& hit : spacer_result.result.hits) {
            cfd_scores.push_back(hit.scoring_info.cfd_score);
        }
        
        // Compute MIT specificity score
        // Note: We don't exclude perfect matches here because in batch mode,
        // we're typically looking at all off-targets, not including the on-target
        spacer_result.mit_specificity_score = 
            cfd::compute_mit_specificity_score(cfd_scores, false);
        
        // Classify specificity tier
        spacer_result.specificity_tier = 
            cfd::classify_specificity_tier(spacer_result.mit_specificity_score);
    }
}

// Sequential batch search implementation (used for single spacer or num_threads=1)
static BatchSearchResult search_genome_batch_sequential(
    const BatchSearchConfig& config, GenomeView view) {

    BatchSearchResult batch_result;
    batch_result.total_hits = 0;
    batch_result.used_gpu = false;
    batch_result.spacers_skipped = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < config.spacers.size(); ++i) {
        const auto& spacer = config.spacers[i];

        if (config.verbose) {
            std::cerr << "[" << (i + 1) << "/" << config.spacers.size()
                      << "] Searching: " << spacer.name << "..." << std::flush;
        }

        // Validate spacer sequence before searching
        try {
            validate_pattern(spacer.sequence);
        } catch (const std::invalid_argument& e) {
            if (config.verbose) {
                std::cerr << " SKIPPED (invalid sequence)\n";
            }
            std::cerr << "Warning: Skipping spacer '" << spacer.name << "'";
            if (!spacer.source_file.empty() && spacer.source_line > 0) {
                std::cerr << " (" << spacer.source_file << ":" << spacer.source_line << ")";
            }
            std::cerr << ": " << e.what() << "\n";
            batch_result.spacers_skipped++;
            continue;  // Skip this spacer, continue with next
        }

        // Build SearchConfig for this spacer
        SearchConfig single_config;
        single_config.pattern = spacer.sequence;
        single_config.threshold = config.threshold;
        single_config.prefer_gpu = config.prefer_gpu;
        single_config.search_both_strands = config.search_both_strands;
        single_config.forward_only = config.forward_only;
        single_config.reverse_only = config.reverse_only;
        single_config.compute_mismatches = config.compute_mismatches;
        single_config.compute_scores = config.compute_scores;
        single_config.pam_filter = config.pam_filter;
        single_config.max_hits = config.max_hits_per_spacer;
        single_config.distance_mode = config.distance_mode;
        single_config.search_start = config.search_start;
        single_config.search_end = config.search_end;

        // Reuse existing search_genome function
        SearchResult result = search_genome(single_config, view);

        if (config.verbose) {
            std::cerr << " " << result.hits.size() << " hits\n";
        }

        BatchSpacerResult spacer_result;
        spacer_result.spacer_name = spacer.name;
        spacer_result.spacer_sequence = spacer.sequence;
        spacer_result.result = std::move(result);

        batch_result.total_hits += spacer_result.result.hits.size();
        batch_result.used_gpu = batch_result.used_gpu || spacer_result.result.used_gpu;

        batch_result.spacer_results.push_back(std::move(spacer_result));
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    batch_result.total_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();

    // Compute MIT specificity scores if requested
    if (config.compute_mit_score && config.compute_scores) {
        compute_batch_mit_scores(batch_result);
    }

    return batch_result;
}

// GPU-optimized batch search implementation using TRUE SPARSE output
// Processes all patterns (+ reverse complements) in batched GPU calls
// This avoids the 3GB-per-pattern dense arrays entirely
static BatchSearchResult search_genome_batch_gpu(
    const BatchSearchConfig& config, GenomeView view) {

    BatchSearchResult batch_result;
    batch_result.total_hits = 0;
    batch_result.used_gpu = true;
    batch_result.spacers_skipped = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    // First pass: validate all spacers and collect valid patterns
    std::vector<std::string> all_patterns;
    std::vector<size_t> valid_spacer_indices;  // Track which spacers are valid
    const bool do_fwd = !config.reverse_only;
    const bool do_rc  = config.reverse_only ||
                        (config.search_both_strands && !config.forward_only);
    const size_t patterns_per_spacer = (do_fwd ? 1 : 0) + (do_rc ? 1 : 0);
    all_patterns.reserve(config.spacers.size() * patterns_per_spacer);
    valid_spacer_indices.reserve(config.spacers.size());

    for (size_t i = 0; i < config.spacers.size(); ++i) {
        const auto& spacer = config.spacers[i];

        // Validate spacer sequence
        try {
            validate_pattern(spacer.sequence);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Warning: Skipping spacer '" << spacer.name << "'";
            if (!spacer.source_file.empty() && spacer.source_line > 0) {
                std::cerr << " (" << spacer.source_file << ":" << spacer.source_line << ")";
            }
            std::cerr << ": " << e.what() << "\n";
            batch_result.spacers_skipped++;
            continue;  // Skip this spacer
        }

        // This spacer is valid
        valid_spacer_indices.push_back(i);
        if (do_fwd) all_patterns.push_back(spacer.sequence);
        if (do_rc)  all_patterns.push_back(reverse_complement(spacer.sequence));
    }
    
    // If no valid spacers, return empty result
    if (all_patterns.empty()) {
        batch_result.total_time_ms = 0.0;
        return batch_result;
    }

    // Resolve search window (0 => whole genome). Kernels take (start, length).
    size_t batch_start = config.search_start;
    size_t batch_end = (config.search_end == 0) ? view.total_bases : config.search_end;
    if (batch_end > view.total_bases) batch_end = view.total_bases;
    if (batch_start > batch_end) batch_start = batch_end;
    size_t batch_length = batch_end - batch_start;

    // Use TRUE SPARSE GPU batch search - no 3GB arrays allocated!
    // This is the key performance optimization: sparse output from GPU kernel
    // is processed directly without expansion to dense arrays.
    // Sparse hits carry ABSOLUTE positions, so no relative-to-absolute fix-up
    // is needed downstream.
    std::vector<GpuSparseResult> gpu_results =
        (config.distance_mode == DistanceMode::HAMMING)
            ? shift_add_gpu_batch_sparse(all_patterns, view, batch_start, batch_length, config.threshold)
            : myers_gpu_batch_sparse    (all_patterns, view, batch_start, batch_length, config.threshold);

    // Process each spacer's results
    batch_result.spacer_results.reserve(config.spacers.size());

    for (size_t valid_idx = 0; valid_idx < valid_spacer_indices.size(); ++valid_idx) {
        size_t spacer_idx = valid_spacer_indices[valid_idx];
        const auto& spacer = config.spacers[spacer_idx];
        
        // Index forward / RC results in the flat patterns list. Each spacer
        // contributed `patterns_per_spacer` patterns; forward (if present) is
        // first, then RC.
        const size_t base_idx = valid_idx * patterns_per_spacer;
        const size_t fwd_idx = base_idx;
        const size_t rc_idx  = base_idx + (do_fwd ? 1 : 0);

        SearchResult result;
        result.pattern = spacer.sequence;
        result.threshold = config.threshold;
        result.total_positions = view.total_bases;
        result.used_gpu = true;

        // Estimate hit count for reservation
        size_t estimated_hits = 0;
        if (do_fwd) estimated_hits += gpu_results[fwd_idx].hits.size();
        if (do_rc)  estimated_hits += gpu_results[rc_idx].hits.size();
        result.hits.reserve(estimated_hits);

        if (do_fwd) {
            for (const auto& sparse_hit : gpu_results[fwd_idx].hits) {
                result.hits.push_back(make_search_hit(
                    view, sparse_hit.position, sparse_hit.distance, Strand::PLUS));
            }
        }
        if (do_rc) {
            for (const auto& sparse_hit : gpu_results[rc_idx].hits) {
                result.hits.push_back(make_search_hit(
                    view, sparse_hit.position, sparse_hit.distance, Strand::MINUS));
            }
        }

        // Sort hits by position
        std::sort(result.hits.begin(), result.hits.end());

        // Thread count for parallel annotation.
        size_t num_threads = config.num_threads;
        if (num_threads == 0) {
            num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 1;
        }

        // Lazy PAM extraction with early filtering (see single-search path).
        bool enable_lazy_pam = (config.pam_filter != PamFilter::NONE &&
                               config.pam_filter != PamFilter::ANY &&
                               config.compute_mismatches);

        if (enable_lazy_pam) {
            // Step 1: Extract PAM only for all hits (cheap - just genome slicing) - PARALLEL
            parallel_process_hits(result.hits,
                [&](SearchHit& hit) {
                    extract_pam_only(view, hit, spacer.sequence.size());
                },
                num_threads);

            // Step 2: Filter by PAM early (discards ~93% of hits for NGG filter)
            std::vector<SearchHit> filtered_hits;
            filtered_hits.reserve(result.hits.size() / 10);  // Expect ~10% to pass

            for (auto& hit : result.hits) {
                if (matches_pam_filter(hit.mismatch_info.pam_type, config.pam_filter)) {
                    filtered_hits.push_back(std::move(hit));
                }
            }

            result.hits = std::move(filtered_hits);

            // Step 3: Extract full mismatch info only for passing hits - PARALLEL
            parallel_process_hits(result.hits,
                [&](SearchHit& hit) {
                    extract_mismatch_info(spacer.sequence, view, hit);
                },
                num_threads);
        } else {
            // Original path: Extract full mismatch info for all hits - PARALLEL
            if (config.compute_mismatches) {
                parallel_process_hits(result.hits,
                    [&](SearchHit& hit) {
                        extract_mismatch_info(spacer.sequence, view, hit);
                    },
                    num_threads);
            }
        }

        // Deduplicate overlapping semi-global alignment hits
        // Deduplication is PAM-aware: prefers hits with valid PAMs when distance is tied
        if (!config.disable_deduplication) {
            result.hits = deduplicate_hits(std::move(result.hits), spacer.sequence.size());
        }

        // Compute CFD scores if requested - PARALLEL
        if (config.compute_scores && config.compute_mismatches) {
            parallel_process_hits(result.hits,
                [&](SearchHit& hit) {
                    hit.scoring_info.cfd_score = cfd::compute_cfd_score(
                        spacer.sequence, hit.mismatch_info);
                    hit.scoring_info.risk_tier = cfd::classify_risk_tier(
                        hit.scoring_info.cfd_score);
                },
                num_threads);
        }

        // PAM filtering already done above in lazy path

        // Apply max_hits limit
        if (config.max_hits_per_spacer > 0 && 
            result.hits.size() > config.max_hits_per_spacer) {
            result.hits.resize(config.max_hits_per_spacer);
        }

        BatchSpacerResult spacer_result;
        spacer_result.spacer_name = spacer.name;
        spacer_result.spacer_sequence = spacer.sequence;
        spacer_result.result = std::move(result);

        batch_result.total_hits += spacer_result.result.hits.size();
        batch_result.spacer_results.push_back(std::move(spacer_result));

        if (config.verbose) {
            std::cerr << "[" << (spacer_idx + 1) << "/" << config.spacers.size()
                      << "] Completed: " << spacer.name
                      << " (" << batch_result.spacer_results.back().result.hits.size()
                      << " hits)\n";
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    batch_result.total_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();

    // Compute MIT specificity scores if requested
    if (config.compute_mit_score && config.compute_scores) {
        compute_batch_mit_scores(batch_result);
    }

    return batch_result;
}

BatchSearchResult search_genome_batch(const BatchSearchConfig& config,
                                      GenomeView view) {
    BatchSearchResult batch_result;
    batch_result.total_hits = 0;
    batch_result.used_gpu = false;

    const size_t num_spacers = config.spacers.size();

    // Handle empty input
    if (num_spacers == 0) {
        batch_result.total_time_ms = 0.0;
        return batch_result;
    }

    // Check if all spacers are GPU-eligible and GPU is preferred
    bool all_gpu_eligible = config.prefer_gpu && gpu_available();
    if (all_gpu_eligible) {
        for (const auto& spacer : config.spacers) {
            if (spacer.sequence.size() > 64) {
                all_gpu_eligible = false;
                break;
            }
        }
    }

    // Use optimized GPU batch path if all spacers are GPU-eligible
    if (all_gpu_eligible) {
        return search_genome_batch_gpu(config, view);
    }

    // Determine thread count
    size_t num_threads = config.num_threads;
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 1;  // Fallback if detection fails
    }

    // Cap threads to spacer count (no point having more threads than work)
    num_threads = std::min(num_threads, num_spacers);

    // Single spacer or single thread: use sequential path for simplicity
    if (num_spacers == 1 || num_threads == 1) {
        return search_genome_batch_sequential(config, view);
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Pre-allocate results vector with correct size (preserves ordering)
    batch_result.spacer_results.resize(num_spacers);

    // Synchronization primitives
    std::mutex gpu_mutex;           // Serialize GPU access
    std::mutex output_mutex;        // Protect verbose output
    std::atomic<size_t> next_index{0};      // Work queue index (atomic for work-stealing)
    std::atomic<size_t> completed{0};       // Progress counter
    std::atomic<size_t> total_hits{0};      // Accumulate total hits
    std::atomic<bool> any_used_gpu{false};  // Track if any spacer used GPU

    // Worker function - each thread runs this
    auto worker = [&]() {
        while (true) {
            // Fetch next work item (atomic increment = work stealing)
            size_t i = next_index.fetch_add(1);
            if (i >= num_spacers) break;  // No more work

            const auto& spacer = config.spacers[i];

            // Build SearchConfig for this spacer
            SearchConfig single_config;
            single_config.pattern = spacer.sequence;
            single_config.threshold = config.threshold;
            single_config.prefer_gpu = config.prefer_gpu;
            single_config.search_both_strands = config.search_both_strands;
            single_config.forward_only = config.forward_only;
            single_config.reverse_only = config.reverse_only;
            single_config.compute_mismatches = config.compute_mismatches;
            single_config.compute_scores = config.compute_scores;
            single_config.pam_filter = config.pam_filter;
            single_config.max_hits = config.max_hits_per_spacer;

            SearchResult result;

            // Check if this spacer will use GPU
            bool will_use_gpu = config.prefer_gpu &&
                               gpu_available() &&
                               spacer.sequence.size() <= 64;

            if (will_use_gpu) {
                // Serialize GPU access to avoid CUDA contention
                std::lock_guard<std::mutex> lock(gpu_mutex);
                result = search_genome(single_config, view);
            } else {
                // CPU path can run fully in parallel
                result = search_genome(single_config, view);
            }

            // Build spacer result
            BatchSpacerResult spacer_result;
            spacer_result.spacer_name = spacer.name;
            spacer_result.spacer_sequence = spacer.sequence;
            spacer_result.result = std::move(result);

            // Update atomic counters
            total_hits += spacer_result.result.hits.size();
            if (spacer_result.result.used_gpu) {
                any_used_gpu.store(true);
            }

            // Store result at correct index (no mutex needed - each index is unique)
            batch_result.spacer_results[i] = std::move(spacer_result);

            // Increment completed counter
            size_t done = completed.fetch_add(1) + 1;

            // Verbose output (thread-safe)
            if (config.verbose) {
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cerr << "[" << done << "/" << num_spacers
                          << "] Completed: " << spacer.name
                          << " (" << batch_result.spacer_results[i].result.hits.size()
                          << " hits)\n";
            }
        }
    };

    // Launch worker threads
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (size_t t = 0; t < num_threads; ++t) {
        threads.emplace_back(worker);
    }

    // Wait for all threads to complete
    for (auto& t : threads) {
        t.join();
    }

    // Finalize batch result from atomics
    batch_result.total_hits = total_hits.load();
    batch_result.used_gpu = any_used_gpu.load();

    auto end_time = std::chrono::high_resolution_clock::now();
    batch_result.total_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();

    // Compute MIT specificity scores if requested
    if (config.compute_mit_score && config.compute_scores) {
        compute_batch_mit_scores(batch_result);
    }

    return batch_result;
}

std::string format_batch_hits_tsv(const BatchSearchResult& result) {
    std::ostringstream oss;

    // Header with spacer column first and new biologically-meaningful columns
    oss << "spacer\tchrom\tstart\tend\tpattern\tdistance\tstrand\t"
        << "aligned_seq\tmismatch_pos\tedit_types\tcigar\t"
        << "pam_seq\tpam_type\t"
        << "seed_edits\tdistal_edits\t"
        << "seed_mismatches\tseed_dna_bulges\tseed_rna_bulges\t"
        << "distal_mismatches\tdistal_dna_bulges\tdistal_rna_bulges\t"
        << "alignment_ambiguous\tn_ambiguous_cells\t"
        << "cfd_score\trisk_tier\t"
        << "mit_specificity\tspecificity_tier\n";

    for (const auto& spacer_result : result.spacer_results) {
        const auto& pattern = spacer_result.spacer_sequence;
        size_t pattern_len = pattern.size();

        for (const auto& hit : spacer_result.result.hits) {
            size_t start = (hit.chrom_offset >= pattern_len - 1)
                           ? hit.chrom_offset - pattern_len + 1
                           : 0;
            size_t end = hit.chrom_offset + 1;

            oss << spacer_result.spacer_name << '\t'
                << hit.chrom_name << '\t'
                << start << '\t'
                << end << '\t'
                << pattern << '\t'
                << static_cast<int>(hit.distance) << '\t'
                << static_cast<char>(hit.strand) << '\t'
                << (hit.mismatch_info.aligned_sequence.empty() ? "." : hit.mismatch_info.aligned_sequence) << '\t'
                << format_mismatch_positions(hit.mismatch_info.mismatch_positions) << '\t'
                << format_edit_types(hit.mismatch_info.edit_types) << '\t'
                << (hit.mismatch_info.cigar.empty() ? "." : hit.mismatch_info.cigar) << '\t'
                << (hit.mismatch_info.pam_sequence.empty() ? "." : hit.mismatch_info.pam_sequence) << '\t'
                << pam_type_str(hit.mismatch_info.pam_type) << '\t'
                << static_cast<int>(hit.mismatch_info.total_seed_edits()) << '\t'
                << static_cast<int>(hit.mismatch_info.total_distal_edits()) << '\t'
                << static_cast<int>(hit.mismatch_info.seed_mismatches) << '\t'
                << static_cast<int>(hit.mismatch_info.seed_dna_bulges) << '\t'
                << static_cast<int>(hit.mismatch_info.seed_rna_bulges) << '\t'
                << static_cast<int>(hit.mismatch_info.distal_mismatches) << '\t'
                << static_cast<int>(hit.mismatch_info.distal_dna_bulges) << '\t'
                << static_cast<int>(hit.mismatch_info.distal_rna_bulges) << '\t'
                << (hit.mismatch_info.alignment_is_ambiguous ? "true" : "false") << '\t'
                << static_cast<int>(hit.mismatch_info.n_ambiguous_cells) << '\t'
                << hit.scoring_info.cfd_score << '\t'
                << risk_tier_str(hit.scoring_info.risk_tier) << '\t'
                << spacer_result.mit_specificity_score << '\t'
                << specificity_tier_str(spacer_result.specificity_tier) << '\n';
        }
    }

    return oss.str();
}

std::string format_batch_hits_bed(const BatchSearchResult& result) {
    std::ostringstream oss;

    for (const auto& spacer_result : result.spacer_results) {
        const auto& pattern = spacer_result.spacer_sequence;
        size_t pattern_len = pattern.size();

        for (const auto& hit : spacer_result.result.hits) {
            size_t start = (hit.chrom_offset >= pattern_len - 1)
                           ? hit.chrom_offset - pattern_len + 1
                           : 0;
            size_t end = hit.chrom_offset + 1;

            // BED score: higher is better, so use (pattern_len - distance)
            int score = static_cast<int>(pattern_len) - static_cast<int>(hit.distance);

            // BED6 + spacer name in column 7
            oss << hit.chrom_name << '\t'
                << start << '\t'
                << end << '\t'
                << pattern << '\t'
                << score << '\t'
                << static_cast<char>(hit.strand) << '\t'
                << spacer_result.spacer_name << '\n';
        }
    }

    return oss.str();
}

std::string format_batch_hits_json(const BatchSearchResult& result) {
    std::ostringstream oss;
    oss.precision(6);
    oss << std::fixed;
    oss << "[";

    bool first = true;
    for (const auto& spacer_result : result.spacer_results) {
        const auto& pattern = spacer_result.spacer_sequence;
        for (const auto& hit : spacer_result.result.hits) {
            if (!first) oss << ",";
            first = false;
            oss << "\n  {"
                << "\"spacer\":\""          << json_escape(spacer_result.spacer_name) << "\","
                << "\"spacer_sequence\":\"" << json_escape(pattern) << "\",";
            write_hit_body_json(oss, hit, pattern);
            oss << ","
                << "\"mit_specificity\":"   << spacer_result.mit_specificity_score << ","
                << "\"specificity_tier\":\"" << specificity_tier_str(spacer_result.specificity_tier) << "\""
                << "}";
        }
    }
    if (!first) oss << "\n";
    oss << "]\n";
    return oss.str();
}

// Summary mode

SummaryResult summarize_batch_results(const BatchSearchResult& batch_result,
                                       uint8_t threshold) {
    SummaryResult summary;
    summary.total_spacers_processed = batch_result.spacer_results.size();
    summary.spacers_skipped = batch_result.spacers_skipped;
    summary.threshold = threshold;

    for (const auto& spacer_result : batch_result.spacer_results) {
        SummaryEntry entry;
        entry.spacer_name = spacer_result.spacer_name;
        entry.spacer_sequence = spacer_result.spacer_sequence;
        entry.total_hits = spacer_result.result.hits.size();

        for (const auto& hit : spacer_result.result.hits) {
            entry.hits_by_distance[hit.distance]++;
        }

        summary.entries.push_back(std::move(entry));
    }

    return summary;
}

SummaryResult summarize_single_result(const SearchResult& result,
                                       const std::string& pattern_name) {
    SummaryResult summary;
    summary.total_spacers_processed = 1;
    summary.spacers_skipped = 0;
    summary.threshold = result.threshold;

    SummaryEntry entry;
    entry.spacer_name = pattern_name;
    entry.spacer_sequence = result.pattern;
    entry.total_hits = result.hits.size();

    for (const auto& hit : result.hits) {
        entry.hits_by_distance[hit.distance]++;
    }

    summary.entries.push_back(std::move(entry));
    return summary;
}

std::string format_summary_json(const SummaryResult& summary) {
    std::ostringstream oss;
    oss << "{\n";
    oss << "  \"threshold\": " << static_cast<int>(summary.threshold) << ",\n";
    oss << "  \"total_spacers_processed\": " << summary.total_spacers_processed << ",\n";
    oss << "  \"spacers_skipped\": " << summary.spacers_skipped << ",\n";
    oss << "  \"spacers\": [\n";

    for (size_t i = 0; i < summary.entries.size(); ++i) {
        const auto& entry = summary.entries[i];
        oss << "    {\n";
        oss << "      \"name\": \"" << entry.spacer_name << "\",\n";
        oss << "      \"sequence\": \"" << entry.spacer_sequence << "\",\n";
        oss << "      \"total_hits\": " << entry.total_hits << ",\n";
        oss << "      \"hits_by_distance\": {";

        bool first = true;
        for (const auto& [dist, count] : entry.hits_by_distance) {
            if (!first) oss << ", ";
            oss << "\"" << static_cast<int>(dist) << "\": " << count;
            first = false;
        }
        oss << "}\n";
        oss << "    }";
        if (i + 1 < summary.entries.size()) oss << ",";
        oss << "\n";
    }

    oss << "  ]\n";
    oss << "}\n";
    return oss.str();
}

std::string format_summary_tsv(const SummaryResult& summary) {
    std::ostringstream oss;

    // Header: spacer, sequence, total_hits, d0, d1, d2, ... d<threshold>
    oss << "spacer\tsequence\ttotal_hits";
    for (uint8_t d = 0; d <= summary.threshold; ++d) {
        oss << "\td" << static_cast<int>(d);
    }
    oss << "\n";

    for (const auto& entry : summary.entries) {
        oss << entry.spacer_name << "\t"
            << entry.spacer_sequence << "\t"
            << entry.total_hits;

        for (uint8_t d = 0; d <= summary.threshold; ++d) {
            auto it = entry.hits_by_distance.find(d);
            oss << "\t" << (it != entry.hits_by_distance.end() ? it->second : 0);
        }
        oss << "\n";
    }

    return oss.str();
}
