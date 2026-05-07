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

// Per-worker batch size for the work-stealing dispatch in
// parallel_process_hits. Each worker grabs PARALLEL_DISPATCH_CHUNK hits per
// atomic increment instead of one. At the 19M-hit / 64-thread regime the
// single-counter atomic was the contention bottleneck (5× observed
// parallelism on a 64-core machine; chunked dispatch lifts that toward the
// memory-bandwidth ceiling). A chunk of 64 cuts atomic traffic 64× while
// keeping the tail-imbalance bounded to <0.001% of the workload.
constexpr size_t PARALLEL_DISPATCH_CHUNK = 64;

}  // namespace

bool cigar_has_indel(const std::string& cigar) {
    return cigar.find('I') != std::string::npos ||
           cigar.find('D') != std::string::npos;
}

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
                                         size_t dedup_radius) {
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
        // within dedup_radius positions of each other.
        // For Levenshtein, dedup_radius = threshold+1 (exact Myers halo size).
        while (cluster_end < hits.size() &&
               hits[cluster_end].strand == hits[i].strand &&
               hits[cluster_end].chrom_name == hits[cluster_end - 1].chrom_name &&
               hits[cluster_end].genome_pos - hits[cluster_end - 1].genome_pos < dedup_radius) {

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

// PAM matching (IUPAC)

// Match one IUPAC pattern base against one observed base.
// Pattern '-' is rejected. Observed 'N' is treated as unknown and does not
// match any concrete pattern base (mirrors search_pipeline's masking semantics).
static bool iupac_match_one(char pat, char seq) {
    pat = static_cast<char>(std::toupper(static_cast<unsigned char>(pat)));
    seq = static_cast<char>(std::toupper(static_cast<unsigned char>(seq)));
    if (seq == 'N') return pat == 'N';
    switch (pat) {
        case 'A': return seq == 'A';
        case 'C': return seq == 'C';
        case 'G': return seq == 'G';
        case 'T': case 'U': return seq == 'T';
        case 'R': return seq == 'A' || seq == 'G';
        case 'Y': return seq == 'C' || seq == 'T';
        case 'S': return seq == 'C' || seq == 'G';
        case 'W': return seq == 'A' || seq == 'T';
        case 'K': return seq == 'G' || seq == 'T';
        case 'M': return seq == 'A' || seq == 'C';
        case 'B': return seq != 'A';
        case 'D': return seq != 'C';
        case 'H': return seq != 'G';
        case 'V': return seq != 'T';
        case 'N': return true;
        default:  return false;
    }
}

// Match an observed PAM (genomic bases) against an IUPAC pattern.
// Returns false when the observation is shorter than the pattern; any
// extraction overhang beyond the pattern length is ignored (so users can
// request wider PAM context than they filter on).
static bool pam_matches_pattern(const std::string& pam, const std::string& pattern) {
    if (pam.size() < pattern.size()) return false;
    for (size_t i = 0; i < pattern.size(); ++i) {
        if (!iupac_match_one(pattern[i], pam[i])) return false;
    }
    return true;
}

// Internal SpCas9 classifier retained for the halo dedup tie-breaker. Applied
// only for 3' PAMs; 5' PAM spacers receive OTHER so the dedup bias stays out of
// tech configurations it wasn't designed for. Slated for removal when the halo
// dedup heuristic is revisited.
static PamType classify_pam_for_dedup(const std::string& pam, PamPosition position) {
    if (position != PamPosition::THREE_PRIME) return PamType::OTHER;
    if (pam.size() < 3) return PamType::INCOMPLETE;
    char p2 = static_cast<char>(std::toupper(static_cast<unsigned char>(pam[1])));
    char p3 = static_cast<char>(std::toupper(static_cast<unsigned char>(pam[2])));
    if (p2 == 'G' && p3 == 'G') return PamType::NGG;
    if (p2 == 'A' && p3 == 'G') return PamType::NAG;
    return PamType::OTHER;
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

    // Parallel path: chunked work-stealing. Each worker grabs
    // PARALLEL_DISPATCH_CHUNK hits per atomic increment to amortize the
    // contention cost at high thread counts. Single-counter dispatch was
    // the bottleneck at 19M hits / 64 threads (only ~5× parallelism observed
    // on a 64-core machine despite no per-hit data dependency); chunked
    // dispatch reduces atomic traffic by PARALLEL_DISPATCH_CHUNK× while
    // preserving work-stealing's load-balancing property.
    std::atomic<size_t> next_chunk_start{0};

    auto worker = [&]() {
        while (true) {
            size_t chunk_start = next_chunk_start.fetch_add(PARALLEL_DISPATCH_CHUNK);
            if (chunk_start >= hit_count) break;
            size_t chunk_end = std::min(chunk_start + PARALLEL_DISPATCH_CHUNK, hit_count);
            for (size_t i = chunk_start; i < chunk_end; ++i) {
                func(hits[i]);
            }
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

// Extract the PAM bases adjacent to a hit, given a pluggable PamSpec.
// Works for both 3'-PAM technologies (Cas9: spacer-PAM) and 5'-PAM technologies
// (Cas12a: PAM-spacer), and for arbitrary PAM length. On the minus strand the
// returned sequence is already reverse-complemented to live in the hit's own
// orientation.
//
// The plus-strand alignment spans [align_begin, genome_pos]. On the minus
// strand the RC'd alignment spans that same window on the forward strand; the
// "3'-PAM" in the hit's orientation therefore sits *before* align_begin on the
// forward strand (RC'd to restore orientation), and the "5'-PAM" sits *after*
// the end of the window. Plus strand is the mirror of that.
// Free-function variant: read PAM bases given only (end_pos, strand). Used
// to pre-filter sparse GPU hits before materializing a full SearchHit.
static std::string read_pam_sequence(const GenomeView& view,
                                     size_t end_pos, Strand strand,
                                     size_t pattern_len,
                                     const PamSpec& pam_spec) {
    const size_t L = pam_spec.effective_extract_length();
    if (L == 0) return {};

    const size_t align_begin = (end_pos >= pattern_len - 1) ? end_pos - pattern_len + 1 : 0;
    const bool three_prime_in_hit_orientation = (pam_spec.position == PamPosition::THREE_PRIME);
    const bool plus = (strand == Strand::PLUS);

    // Forward-strand side where the PAM bytes live:
    //   plus + 3'-PAM  -> after the window
    //   plus + 5'-PAM  -> before the window
    //   minus + 3'-PAM -> before the window (then RC)
    //   minus + 5'-PAM -> after the window (then RC)
    const bool pam_after_window = (plus && three_prime_in_hit_orientation) ||
                                  (!plus && !three_prime_in_hit_orientation);

    std::string bases;
    if (pam_after_window) {
        size_t pam_start = end_pos + 1;
        if (pam_start >= view.total_bases) return {};
        size_t avail = view.total_bases - pam_start;
        size_t take  = std::min(avail, L);
        bases = extract_genome_slice(view, pam_start, take);
    } else {
        if (align_begin == 0) return {};
        size_t take = std::min(align_begin, L);
        bases = extract_genome_slice(view, align_begin - take, take);
    }

    if (!plus) bases = reverse_complement(bases);
    return bases;
}

static void extract_pam_bases(const GenomeView& view, SearchHit& hit,
                              size_t pattern_len, const PamSpec& pam_spec) {
    std::string bases = read_pam_sequence(view, hit.genome_pos, hit.strand,
                                           pattern_len, pam_spec);
    hit.mismatch_info.pam_sequence = std::move(bases);
    hit.mismatch_info.pam_type =
        classify_pam_for_dedup(hit.mismatch_info.pam_sequence, pam_spec.position);
}

// Verify that `pattern` admits an alignment of edit distance ≤ threshold to a
// target of *exactly* `tl` bases ending at `end_pos` in the hit's orientation.
// Standard global edit-distance DP with row-min early exit. Used by the Lev
// PAM post-filter: the widened PAM filter optimistically accepts hits whose
// PAM matches at any candidate target length in [pattern_len ± threshold],
// but does not check that the alignment at that length is actually valid.
// For substitution-only canonical hits whose first `|PAM|` aligned bases
// happen to look like a PAM (e.g. spacer 3' tail = NGG on minus strand), the
// optimistic accept becomes a false positive. This helper closes the gap by
// confirming that an alignment exists at the target length the PAM came from.
static bool alignment_exists_at_target_length(
        const std::string& pattern,
        const GenomeView& view,
        size_t end_pos,
        Strand strand,
        size_t tl,
        int threshold) {
    if (tl == 0) return false;
    if (end_pos + 1 < tl) return false;
    const size_t target_start = end_pos - tl + 1;
    if (target_start + tl > view.total_bases) return false;

    std::string target = extract_genome_slice(view, target_start, tl);
    if (strand == Strand::MINUS) target = reverse_complement(target);

    const size_t m = pattern.size();
    const size_t n = target.size();

    std::vector<int> prev(n + 1), curr(n + 1);
    for (size_t j = 0; j <= n; ++j) prev[j] = static_cast<int>(j);

    for (size_t i = 1; i <= m; ++i) {
        curr[0] = static_cast<int>(i);
        int row_min = curr[0];
        const char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pattern[i - 1])));
        for (size_t j = 1; j <= n; ++j) {
            const char t = static_cast<char>(std::toupper(static_cast<unsigned char>(target[j - 1])));
            int match_cost;
            if (p == 'N') match_cost = 0;
            else if (t == 'N') match_cost = 1;
            else match_cost = (p == t) ? 0 : 1;
            curr[j] = std::min({
                prev[j - 1] + match_cost,
                prev[j] + 1,
                curr[j - 1] + 1
            });
            if (curr[j] < row_min) row_min = curr[j];
        }
        if (row_min > threshold) return false;
        prev.swap(curr);
    }
    return prev[n] <= threshold;
}

// Forward declarations for use in repopulate_hit_at_tl.
static std::string format_cigar(const std::vector<AlignOp>& alignment);
static std::string build_aligned_sequence_from_alignment(const std::vector<AlignOp>& alignment,
                                                         const std::string& compare_sequence,
                                                         size_t align_text_start);

// Pinned-tl global alignment with traceback. Used by the v0.10.0 Lev
// post-filter: when a hit is accepted on the basis of an alt-tl PAM
// match, the canonical (best-d) alignment doesn't describe the tl that
// justifies keeping the hit, so we re-derive distance/CIGAR/aligned_seq
// from the alt-tl alignment instead.
//
// Both ends pinned (no free start/end gaps). Pattern fully consumed,
// target (length `tl`) fully consumed. Pattern N is wildcard (cost 0);
// genome N is masked (cost 1).
//
// Returns the alignment (in MATCH > INS_TEXT > INS_PATTERN tie-break
// order) and ambiguity flag. If d > threshold, alignment is empty and
// distance is set to threshold+1 as a sentinel.
struct PinnedAlignResult {
    int distance;
    std::vector<AlignOp> alignment;  // forward order; covers all pattern + all target
    bool has_ambiguity;
    size_t n_ambiguous_cells;
};

static PinnedAlignResult align_pinned_tl(const std::string& pattern,
                                          const std::string& target,
                                          int threshold) {
    PinnedAlignResult out{};
    out.distance = threshold + 1;

    const size_t m = pattern.size();
    const size_t n = target.size();

    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));
    for (size_t i = 0; i <= m; ++i) dp[i][0] = static_cast<int>(i);
    for (size_t j = 0; j <= n; ++j) dp[0][j] = static_cast<int>(j);

    for (size_t i = 1; i <= m; ++i) {
        const char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pattern[i - 1])));
        int row_min = std::numeric_limits<int>::max();
        for (size_t j = 1; j <= n; ++j) {
            const char t = static_cast<char>(std::toupper(static_cast<unsigned char>(target[j - 1])));
            int match_cost;
            if (p == 'N') match_cost = 0;
            else if (t == 'N') match_cost = 1;
            else match_cost = (p == t) ? 0 : 1;
            dp[i][j] = std::min({
                dp[i - 1][j - 1] + match_cost,
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1
            });
            if (dp[i][j] < row_min) row_min = dp[i][j];
        }
        // Early exit if every cell in this row already exceeds threshold.
        if (row_min > threshold) {
            return out;
        }
    }

    if (dp[m][n] > threshold) {
        return out;
    }
    out.distance = dp[m][n];

    // Traceback from (m, n) with tie-breaking: MATCH > INS_TEXT (DNA bulge)
    // > INS_PATTERN (RNA bulge). Mirrors compute_alignment().
    std::vector<AlignOp> alignment;
    alignment.reserve(m + n);
    bool has_ambig = false;
    size_t n_ambig = 0;
    size_t i = m, j = n;
    while (i > 0 && j > 0) {
        const char p = static_cast<char>(std::toupper(static_cast<unsigned char>(pattern[i - 1])));
        const char t = static_cast<char>(std::toupper(static_cast<unsigned char>(target[j - 1])));
        int match_cost;
        if (p == 'N') match_cost = 0;
        else if (t == 'N') match_cost = 1;
        else match_cost = (p == t) ? 0 : 1;
        const int current = dp[i][j];

        const bool can_match = (current == dp[i - 1][j - 1] + match_cost);
        const bool can_ins_text = (current == dp[i][j - 1] + 1);
        const bool can_ins_pattern = (current == dp[i - 1][j] + 1);
        const int n_optimal = (int)can_match + (int)can_ins_text + (int)can_ins_pattern;
        if (n_optimal > 1) { has_ambig = true; ++n_ambig; }

        if (can_match) {
            alignment.push_back(AlignOp::MATCH);
            --i; --j;
        } else if (can_ins_text) {
            alignment.push_back(AlignOp::INS_TEXT);
            --j;
        } else if (can_ins_pattern) {
            alignment.push_back(AlignOp::INS_PATTERN);
            --i;
        } else {
            break;  // shouldn't happen with correct DP
        }
    }
    while (i > 0) { alignment.push_back(AlignOp::INS_PATTERN); --i; }
    while (j > 0) { alignment.push_back(AlignOp::INS_TEXT);    --j; }
    std::reverse(alignment.begin(), alignment.end());

    out.alignment = std::move(alignment);
    out.has_ambiguity = has_ambig;
    out.n_ambiguous_cells = n_ambig;
    return out;
}

// Re-derive a hit's mismatch_info from a pinned target-length alignment.
// Replaces hit.distance, hit.mismatch_info.{cigar, aligned_sequence,
// alignment_is_ambiguous, n_ambiguous_cells, pam_sequence, pam_type}.
//
// Used by the v0.10.0 Lev post-filter (Option A): when canonical PAM
// (at text_consumed = best_d alignment's tl) doesn't match but an alt-tl
// alignment does, this re-populates the hit so that distance/CIGAR/PAM
// describe the SAME alignment topology — biologically the alignment
// that justifies keeping the hit. Without this, aCFD scoring uses the
// canonical (lower) distance while PAM display reflects the alt-tl
// alignment, inflating the score for hits that are actually d > 0.
//
// Returns true if a valid d ≤ threshold alignment exists at this tl;
// false (and leaves hit unchanged) otherwise.
static bool repopulate_hit_at_tl(
        const std::string& original_pattern,
        GenomeView view,
        const PamSpec& pam_spec,
        SearchHit& hit,
        size_t pinned_tl,
        int threshold) {
    if (pinned_tl == 0) return false;
    if (hit.genome_pos + 1 < pinned_tl) return false;
    const size_t target_start = hit.genome_pos - pinned_tl + 1;
    if (target_start + pinned_tl > view.total_bases) return false;

    std::string genome_seq = extract_genome_slice(view, target_start, pinned_tl);
    std::string compare_sequence = (hit.strand == Strand::MINUS)
        ? reverse_complement(genome_seq) : genome_seq;

    PinnedAlignResult ar = align_pinned_tl(original_pattern, compare_sequence, threshold);
    if (ar.distance > threshold) return false;

    hit.distance = static_cast<uint8_t>(ar.distance);
    hit.mismatch_info.cigar = format_cigar(ar.alignment);
    hit.mismatch_info.aligned_sequence =
        build_aligned_sequence_from_alignment(ar.alignment, compare_sequence, 0);
    hit.mismatch_info.alignment_is_ambiguous = ar.has_ambiguity;
    hit.mismatch_info.n_ambiguous_cells = static_cast<uint16_t>(
        std::min<size_t>(ar.n_ambiguous_cells, std::numeric_limits<uint16_t>::max()));
    // PAM at this pinned tl. extract_pam_bases reads pam_pat_len bases and
    // populates pam_sequence + pam_type.
    extract_pam_bases(view, hit, pinned_tl, pam_spec);
    return true;
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

static void extract_mismatch_info(const std::string& original_pattern,
                                  GenomeView view,
                                  const PamSpec& pam_spec,
                                  SearchHit& hit) {
    const size_t m = original_pattern.size();
    hit.mismatch_info.alignment_is_ambiguous = false;

    const size_t max_align_len = m + hit.distance + 1;
    const size_t align_start = (hit.genome_pos >= max_align_len - 1)
                               ? hit.genome_pos - max_align_len + 1 : 0;
    const size_t align_len   = hit.genome_pos - align_start + 1;
    const std::string genome_seq = extract_genome_slice(view, align_start, align_len);

    const std::string& compare_pattern = original_pattern;
    std::string compare_sequence = (hit.strand == Strand::MINUS)
        ? reverse_complement(genome_seq) : genome_seq;

    if (hit.distance == 0) {
        // No bulges possible → target length == pattern length.
        if (hit.mismatch_info.pam_sequence.empty()) {
            extract_pam_bases(view, hit, m, pam_spec);
        }
        hit.mismatch_info.aligned_sequence =
            extract_aligned_region_perfect_match(compare_sequence, hit.strand, m);
        hit.mismatch_info.cigar = std::to_string(m) + "M";
        return;
    }

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

    // v0.7.0: extract PAM using the actual aligned-target length (text_consumed),
    // which may differ from the nominal pattern length m when target-side bulges
    // shortened the target or query-side bulges lengthened it. Any pam_sequence
    // that might have been set earlier is overwritten here so it reflects the
    // true alignment geometry, not the filter's nominal-length guess.
    if (text_consumed > 0) {
        extract_pam_bases(view, hit, text_consumed, pam_spec);
    } else if (hit.mismatch_info.pam_sequence.empty()) {
        extract_pam_bases(view, hit, m, pam_spec);
    }

    hit.mismatch_info.aligned_sequence = build_aligned_sequence_from_alignment(
        ar.alignment, compare_sequence, align_text_start);
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

        // Build SearchHit objects with chromosome info.
        // For Hamming (bitap), positions rel_pos < pattern_len-1 are prefix-match
        // artifacts from the start of the text, not full-length alignments — skip them.
        const size_t min_valid_pos = (config.distance_mode == DistanceMode::HAMMING)
                                     ? fwd_pattern.size() - 1
                                     : 0;
        std::vector<SearchHit> hits;
        hits.reserve(hit_positions.size());

        for (size_t rel_pos : hit_positions) {
            if (rel_pos < min_valid_pos) continue;
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

    // PAM filter (cheap, per-hit) is independent of per-hit mismatch
    // extraction (CIGAR, aligned seq — expensive). Always filter when a
    // PAM is set; compute mismatch info separately only when requested.
    //
    // v0.7.1: for Levenshtein with target-side bulges, the nominal pattern
    // length is wrong for the PAM position on the bulge-sensitive side
    // (plus-strand 5' PAM or minus-strand 3' PAM). Check a window of
    // candidate target lengths and accept if any matches. See the batch
    // path for the same logic and correctness notes.
    if (config.pam.filtering_enabled()) {
        const bool is_lev = (config.distance_mode == DistanceMode::LEVENSHTEIN);
        const size_t pam_filter_threshold = is_lev
            ? static_cast<size_t>(std::max(0, static_cast<int>(config.threshold)))
            : 0;
        const bool three_prime_pam =
            (config.pam.position == PamPosition::THREE_PRIME);
        const size_t pattern_len = config.pattern.size();
        auto pam_sensitive_to_bulges = [&](Strand strand) {
            const bool plus = (strand == Strand::PLUS);
            const bool pam_after_window =
                (plus && three_prime_pam) || (!plus && !three_prime_pam);
            return !pam_after_window;
        };

        parallel_process_hits(result.hits,
            [&](SearchHit& hit) {
                // Start with canonical PAM (nominal pattern_len).
                extract_pam_bases(view, hit, pattern_len, config.pam);
                if (!pam_matches_pattern(hit.mismatch_info.pam_sequence,
                                         config.pam.pattern)
                    && pam_filter_threshold > 0
                    && pam_sensitive_to_bulges(hit.strand)) {
                    // Try shifted target lengths (±threshold).
                    const size_t lo = (pattern_len > pam_filter_threshold)
                        ? pattern_len - pam_filter_threshold : 1;
                    const size_t hi = pattern_len + pam_filter_threshold;
                    for (size_t tl = lo; tl <= hi; ++tl) {
                        std::string candidate = read_pam_sequence(
                            view, hit.genome_pos, hit.strand, tl, config.pam);
                        if (pam_matches_pattern(candidate, config.pam.pattern)) {
                            hit.mismatch_info.pam_sequence = std::move(candidate);
                            hit.mismatch_info.pam_type = classify_pam_for_dedup(
                                hit.mismatch_info.pam_sequence,
                                config.pam.position);
                            break;
                        }
                    }
                }
            },
            num_threads);

        std::vector<SearchHit> filtered_hits;
        filtered_hits.reserve(result.hits.size() / 10);

        for (auto& hit : result.hits) {
            if (pam_matches_pattern(hit.mismatch_info.pam_sequence, config.pam.pattern)) {
                filtered_hits.push_back(std::move(hit));
            }
        }

        result.hits = std::move(filtered_hits);
    }

    if (config.compute_mismatches) {
        parallel_process_hits(result.hits,
            [&](SearchHit& hit) {
                extract_mismatch_info(config.pattern, view, config.pam, hit);
            },
            num_threads);

        // v0.10.0 Lev PAM re-filter (single-pattern path). See batch path
        // for full rationale. Picks the LOWEST-distance alt-tl among
        // those with both (PAM matches at this tl) and (valid alignment
        // at this tl).
        if (config.pam.filtering_enabled() &&
            config.distance_mode == DistanceMode::LEVENSHTEIN) {
            const size_t pattern_len = config.pattern.size();
            const size_t pam_filter_threshold =
                static_cast<size_t>(std::max(0, static_cast<int>(config.threshold)));
            const bool three_prime_pam =
                (config.pam.position == PamPosition::THREE_PRIME);
            auto pam_sensitive_to_bulges = [&](Strand strand) {
                const bool plus = (strand == Strand::PLUS);
                const bool pam_after_window =
                    (plus && three_prime_pam) || (!plus && !three_prime_pam);
                return !pam_after_window;
            };

            std::vector<SearchHit> kept;
            kept.reserve(result.hits.size());
            for (auto& h : result.hits) {
                if (pam_matches_pattern(h.mismatch_info.pam_sequence,
                                        config.pam.pattern)) {
                    kept.push_back(std::move(h));
                    continue;
                }
                if (!pam_sensitive_to_bulges(h.strand) || pam_filter_threshold == 0) {
                    continue;
                }
                const size_t lo = (pattern_len > pam_filter_threshold)
                    ? pattern_len - pam_filter_threshold : 1;
                const size_t hi = pattern_len + pam_filter_threshold;
                SearchHit best_hit;
                bool have_best = false;
                for (size_t tl = lo; tl <= hi; ++tl) {
                    if (tl == pattern_len) continue;
                    std::string candidate = read_pam_sequence(
                        view, h.genome_pos, h.strand, tl, config.pam);
                    if (!pam_matches_pattern(candidate, config.pam.pattern)) continue;
                    SearchHit trial = h;
                    if (!repopulate_hit_at_tl(
                            config.pattern, view, config.pam, trial, tl,
                            config.threshold)) continue;
                    if (!have_best || trial.distance < best_hit.distance) {
                        best_hit = std::move(trial);
                        have_best = true;
                    }
                }
                if (have_best) kept.push_back(std::move(best_hit));
            }
            result.hits = std::move(kept);
        }
    }

    // Deduplicate overlapping alignment hits (Levenshtein only).
    // Hamming mode has no halo — each position is independent; skip dedup.
    // Levenshtein: Myers halo radius = threshold positions from the true match
    // end; use threshold+1 to safely collapse halos without dropping nearby
    // genuine distinct sites.
    if (!config.disable_deduplication && config.distance_mode == DistanceMode::LEVENSHTEIN) {
        result.hits = deduplicate_hits(std::move(result.hits),
                                       static_cast<size_t>(config.threshold) + 1);
    }

    // Compute CFD scores if requested (requires mismatch info) - PARALLEL
    if (config.compute_scores && config.compute_mismatches) {
        auto score_t0 = std::chrono::high_resolution_clock::now();
        parallel_process_hits(result.hits,
            [&](SearchHit& hit) {
                hit.scoring_info.cfd_score = cfd::compute_cfd_score(
                    config.pattern, hit.mismatch_info);
            },
            num_threads);
        result.scoring_time_ms = std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - score_t0).count();
    }

    // PAM filtering already done above in lazy path

    // Apply max_hits limit if specified
    if (config.max_hits > 0 && result.hits.size() > config.max_hits) {
        result.hits.resize(config.max_hits);
    }

    return result;
}

// Output formatting

std::string format_hits_tsv(const std::vector<SearchHit>& hits,
                            const std::string& pattern) {
    std::ostringstream oss;

    oss << "chrom\tstart\tend\tpattern\tdistance\tstrand\t"
        << "aligned_seq\tcigar\tpam_seq\t"
        << "alignment_ambiguous\tn_ambiguous_cells\t"
        << "cfd_score\n";

    size_t pattern_len = pattern.size();
    for (const auto& hit : hits) {
        size_t start = (hit.chrom_offset >= pattern_len - 1)
                       ? hit.chrom_offset - pattern_len + 1
                       : 0;
        size_t end = hit.chrom_offset + 1;

        oss << hit.chrom_name << '\t'
            << start << '\t'
            << end << '\t'
            << pattern << '\t'
            << static_cast<int>(hit.distance) << '\t'
            << static_cast<char>(hit.strand) << '\t'
            << (hit.mismatch_info.aligned_sequence.empty() ? "." : hit.mismatch_info.aligned_sequence) << '\t'
            << (hit.mismatch_info.cigar.empty() ? "." : hit.mismatch_info.cigar) << '\t'
            << (hit.mismatch_info.pam_sequence.empty() ? "." : hit.mismatch_info.pam_sequence) << '\t'
            << (hit.mismatch_info.alignment_is_ambiguous ? "true" : "false") << '\t'
            << static_cast<int>(hit.mismatch_info.n_ambiguous_cells) << '\t'
            << hit.scoring_info.cfd_score << '\n';
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

    oss << "\"cigar\":";
    if (mi.cigar.empty()) oss << "null";
    else oss << "\"" << json_escape(mi.cigar) << "\"";
    oss << ",";

    oss << "\"pam_seq\":";
    if (mi.pam_sequence.empty()) oss << "null";
    else oss << "\"" << json_escape(mi.pam_sequence) << "\"";
    oss << ",";

    oss << "\"alignment_ambiguous\":"  << (mi.alignment_is_ambiguous ? "true" : "false") << ","
        << "\"n_ambiguous_cells\":"    << static_cast<int>(mi.n_ambiguous_cells) << ","
        << "\"cfd_score\":"            << hit.scoring_info.cfd_score;
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
        single_config.pam = config.pam;
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

        // Thread count for parallel annotation.
        size_t num_threads = config.num_threads;
        if (num_threads == 0) {
            num_threads = std::thread::hardware_concurrency();
            if (num_threads == 0) num_threads = 1;
        }

        // Early PAM filter: when a PAM is set, drop non-matching sparse hits
        // *before* materializing heavyweight SearchHit structs. On human
        // Ham+NGG this skips ~98% of would-be SearchHit constructions (19M
        // sparse hits → ~400K keepers), saving several seconds of allocation
        // + string-construction overhead per 225-spacer run.
        //
        // Correctness note (v0.7.0): for Levenshtein mode, target-side bulges
        // can shorten the aligned target region below `pattern_len`, shifting
        // the PAM-adjacent position closer to the spacer. The true alignment
        // length is in [pattern_len - threshold, pattern_len + threshold].
        // Cases where PAM sits adjacent to the *align_begin* side of the
        // window (plus-strand 5' PAM, minus-strand 3' PAM) depend on this
        // length; the earlier v0.6.5 filter used only the nominal pattern_len
        // and silently dropped ~28% of legitimate indel-containing hits.
        // The fix below widens the filter's PAM search to cover every
        // candidate alignment length in the valid range. The filter still
        // conservatively drops only sparse hits that CANNOT host a valid PAM
        // under any bulge configuration. The authoritative PAM sequence is
        // re-derived from the full DP alignment in extract_mismatch_info()
        // (where text_consumed is known exactly).
        const bool do_pam_filter = config.pam.filtering_enabled();
        const size_t pattern_len = spacer.sequence.size();
        const bool is_lev = (config.distance_mode == DistanceMode::LEVENSHTEIN);
        const size_t pam_filter_threshold = is_lev
            ? static_cast<size_t>(std::max(0, static_cast<int>(config.threshold)))
            : 0;

        // Determine whether the PAM position depends on align_begin (which
        // shifts with target bulges) vs the fixed end_pos (which does not).
        // The PAM filter only needs to widen its check for the former cases.
        const bool three_prime_pam =
            (config.pam.position == PamPosition::THREE_PRIME);
        // PAM sits adjacent to align_begin (bulge-sensitive) on:
        //   plus-strand + 5'-PAM  or  minus-strand + 3'-PAM
        // PAM sits adjacent to end_pos (bulge-insensitive) on:
        //   plus-strand + 3'-PAM  or  minus-strand + 5'-PAM
        auto pam_sensitive_to_bulges = [&](Strand strand) {
            const bool plus = (strand == Strand::PLUS);
            const bool pam_after_window =
                (plus && three_prime_pam) || (!plus && !three_prime_pam);
            return !pam_after_window;
        };

        auto filter_pam_match = [&](size_t end_pos, Strand strand) -> bool {
            if (!pam_filter_threshold || !pam_sensitive_to_bulges(strand)) {
                // Fast path: PAM position independent of bulges, or Hamming mode.
                std::string pam_seq = read_pam_sequence(
                    view, end_pos, strand, pattern_len, config.pam);
                return pam_matches_pattern(pam_seq, config.pam.pattern);
            }
            // Lev mode + bulge-sensitive side: check every candidate target
            // length in [pattern_len - threshold, pattern_len + threshold].
            const size_t lo = (pattern_len > pam_filter_threshold)
                ? pattern_len - pam_filter_threshold : 1;
            const size_t hi = pattern_len + pam_filter_threshold;
            for (size_t target_len = lo; target_len <= hi; ++target_len) {
                std::string pam_seq = read_pam_sequence(
                    view, end_pos, strand, target_len, config.pam);
                if (pam_matches_pattern(pam_seq, config.pam.pattern)) return true;
            }
            return false;
        };

        auto materialize_and_filter = [&](const std::vector<GpuSparseHit>& sparse,
                                          Strand strand) {
            for (const auto& sparse_hit : sparse) {
                if (do_pam_filter) {
                    if (!filter_pam_match(sparse_hit.position, strand)) continue;
                    // Do NOT pre-populate pam_sequence here: for Lev with
                    // bulges, the correct PAM depends on the actual alignment
                    // length, which extract_mismatch_info derives below. In
                    // Hamming mode pam_sequence is re-read there too (cheap),
                    // so dropping the pre-populate is uniform.
                    result.hits.push_back(make_search_hit(
                        view, sparse_hit.position, sparse_hit.distance, strand));
                } else {
                    result.hits.push_back(make_search_hit(
                        view, sparse_hit.position, sparse_hit.distance, strand));
                }
            }
        };

        // Reserve: for PAM-filtered runs, most hits drop — don't over-reserve.
        size_t estimated_hits = 0;
        if (do_fwd) estimated_hits += gpu_results[fwd_idx].hits.size();
        if (do_rc)  estimated_hits += gpu_results[rc_idx].hits.size();
        if (do_pam_filter) estimated_hits = estimated_hits / 10 + 1024;
        result.hits.reserve(estimated_hits);

        if (do_fwd) materialize_and_filter(gpu_results[fwd_idx].hits, Strand::PLUS);
        if (do_rc)  materialize_and_filter(gpu_results[rc_idx].hits,  Strand::MINUS);

        // Sort hits by position
        std::sort(result.hits.begin(), result.hits.end());

        if (config.compute_mismatches) {
            parallel_process_hits(result.hits,
                [&](SearchHit& hit) {
                    extract_mismatch_info(spacer.sequence, view, config.pam, hit);
                },
                num_threads);

            // v0.10.0 Lev PAM post-filter — full re-derivation when
            // canonical PAM doesn't match (Option A).
            //
            // When the canonical (best-d) alignment's PAM doesn't match,
            // search every alt-tl with matching PAM and re-derive the
            // alignment from there. Pick the alt-tl with LOWEST distance
            // (among PAM-matching ones with valid alignments). Without
            // the lowest-d preference, the post-filter would pick an
            // arbitrary alt-tl and inflate the displayed distance — and
            // halo dedup would then collapse this hit with neighboring
            // canonical-match hits using stale (high) distance, picking
            // the wrong representative.
            if (do_pam_filter && is_lev && pam_filter_threshold > 0) {
                const size_t L = config.pam.effective_extract_length();
                std::vector<SearchHit> kept;
                kept.reserve(result.hits.size());
                for (auto& h : result.hits) {
                    // Fast path: canonical alignment already has matching PAM.
                    if (pam_matches_pattern(h.mismatch_info.pam_sequence,
                                            config.pam.pattern)) {
                        kept.push_back(std::move(h));
                        continue;
                    }
                    if (!pam_sensitive_to_bulges(h.strand)) {
                        // PAM bulge-insensitive: canonical PAM is
                        // authoritative, and it didn't match → drop.
                        continue;
                    }
                    const size_t lo = (pattern_len > pam_filter_threshold)
                        ? pattern_len - pam_filter_threshold : 1;
                    const size_t hi = pattern_len + pam_filter_threshold;
                    // Find the alt-tl with the LOWEST distance whose
                    // alignment is valid AND whose PAM matches. We do
                    // this in two passes so we don't waste time on the
                    // expensive repopulate when a closer tl is around.
                    SearchHit best_hit;
                    bool have_best = false;
                    for (size_t tl = lo; tl <= hi; ++tl) {
                        if (tl == pattern_len) continue;
                        std::string candidate = read_pam_sequence(
                            view, h.genome_pos, h.strand, tl, config.pam);
                        if (candidate.size() < L) continue;
                        if (!pam_matches_pattern(candidate, config.pam.pattern)) continue;
                        // Try the re-derivation. If alignment isn't valid
                        // at this tl (d > threshold), skip.
                        SearchHit trial = h;
                        if (!repopulate_hit_at_tl(
                                spacer.sequence, view, config.pam, trial, tl,
                                config.threshold)) continue;
                        if (!have_best || trial.distance < best_hit.distance) {
                            best_hit = std::move(trial);
                            have_best = true;
                        }
                    }
                    if (have_best) kept.push_back(std::move(best_hit));
                }
                result.hits = std::move(kept);
            }
        }

        // Deduplicate overlapping alignment hits (Levenshtein only).
        // Hamming mode has no halo — skip. Levenshtein: radius = threshold+1.
        if (!config.disable_deduplication && config.distance_mode == DistanceMode::LEVENSHTEIN) {
            result.hits = deduplicate_hits(std::move(result.hits),
                                           static_cast<size_t>(config.threshold) + 1);
        }

        // Compute CFD scores if requested - PARALLEL
        if (config.compute_scores && config.compute_mismatches) {
            auto score_t0 = std::chrono::high_resolution_clock::now();
            parallel_process_hits(result.hits,
                [&](SearchHit& hit) {
                    hit.scoring_info.cfd_score = cfd::compute_cfd_score(
                        spacer.sequence, hit.mismatch_info);
                },
                num_threads);
            batch_result.scoring_time_ms += std::chrono::duration<double, std::milli>(
                std::chrono::high_resolution_clock::now() - score_t0).count();
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
    std::atomic<size_t> next_index{0};       // Work queue index (atomic for work-stealing)
    std::atomic<size_t> completed{0};        // Progress counter
    std::atomic<size_t> total_hits{0};       // Accumulate total hits
    std::atomic<bool>   any_used_gpu{false}; // Track if any spacer used GPU
    std::atomic<int64_t> total_scoring_us{0}; // Accumulate CFD scoring time (microseconds)

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
            single_config.pam = config.pam;
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
            total_scoring_us += static_cast<int64_t>(
                spacer_result.result.scoring_time_ms * 1000.0);

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
    batch_result.scoring_time_ms = total_scoring_us.load() / 1000.0;

    auto end_time = std::chrono::high_resolution_clock::now();
    batch_result.total_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();

    return batch_result;
}

std::string format_batch_hits_tsv(const BatchSearchResult& result) {
    std::ostringstream oss;

    oss << "spacer\tchrom\tstart\tend\tpattern\tdistance\tstrand\t"
        << "aligned_seq\tcigar\tpam_seq\t"
        << "alignment_ambiguous\tn_ambiguous_cells\t"
        << "cfd_score\n";

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
                << (hit.mismatch_info.cigar.empty() ? "." : hit.mismatch_info.cigar) << '\t'
                << (hit.mismatch_info.pam_sequence.empty() ? "." : hit.mismatch_info.pam_sequence) << '\t'
                << (hit.mismatch_info.alignment_is_ambiguous ? "true" : "false") << '\t'
                << static_cast<int>(hit.mismatch_info.n_ambiguous_cells) << '\t'
                << hit.scoring_info.cfd_score << '\n';
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
            oss << "}";
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
