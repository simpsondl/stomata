#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <genome_loader.hpp>   // Genome, Chromosome

// Version

static constexpr const char* STOMATA_VERSION = "0.11.0";

// Data structures

// Strand orientation of a hit.
enum class Strand : char {
    PLUS  = '+',
    MINUS = '-'
};

// PAM position relative to the spacer.
// Cas9 / Cas12 / base editors disagree on which side of the spacer the PAM sits;
// Stomata is PAM-agnostic, so the user specifies the position along with the pattern.
enum class PamPosition : uint8_t {
    THREE_PRIME,  // PAM follows the spacer (e.g. SpCas9: 5'-[spacer]-NGG-3')
    FIVE_PRIME    // PAM precedes the spacer (e.g. Cas12a: 5'-TTTV-[spacer]-3')
};

// Pluggable PAM specification: IUPAC pattern + position + extraction length.
// Empty pattern means no PAM filter and no extraction.
// extract_length > 0 overrides the auto-derived length (= pattern.size()); use it
// when you want to annotate PAM context wider than the filter pattern.
struct PamSpec {
    std::string  pattern;
    PamPosition  position = PamPosition::THREE_PRIME;
    size_t       extract_length = 0;

    bool filtering_enabled() const { return !pattern.empty(); }
    size_t effective_extract_length() const {
        return extract_length > 0 ? extract_length : pattern.size();
    }
};

// Internal PAM classification used only by the (SpCas9-biased) halo deduplicator.
// TODO: superseded when the dedup logic is revisited — PAM priority here is
// tech-specific and contradicts the PAM-agnostic design. Kept dormant for
// non-3'-PAM searches (see classify_pam in search_pipeline.cpp).
enum class PamType : uint8_t {
    NGG,
    NAG,
    OTHER,
    INCOMPLETE
};

// Distance metric mode for pattern matching.
// Controls which algorithm and distance definition to use.
enum class DistanceMode : uint8_t {
    LEVENSHTEIN,  // Full edit distance (substitutions + insertions + deletions) - default
    HAMMING       // Substitution-only distance (no indels/bulges) - faster, simpler
};

// Alignment operation for traceback.
enum class AlignOp : uint8_t {
    MATCH,      // Match or mismatch (diagonal move)
    INS_TEXT,   // Insertion in text/genome (horizontal move) = DNA bulge
    INS_PATTERN // Insertion in pattern/RNA (vertical move) = RNA bulge
};

// Per-hit alignment annotation.
// Position-level detail is carried by the CIGAR string (M = match/mismatch,
// I = RNA bulge, D = DNA bulge); downstream callers parse the CIGAR for
// anything finer-grained.
struct MismatchInfo {
    std::string aligned_sequence;                 // Genome sequence at hit (spacer region)
    std::string pam_sequence;                     // Extracted PAM bases (length = PamSpec extract length)
    std::string cigar;                            // M / I / D run-length CIGAR

    bool     alignment_is_ambiguous = false;      // True if multiple optimal alignments exist
    uint16_t n_ambiguous_cells      = 0;          // Traceback cells where >1 op was co-optimal

    PamType  pam_type = PamType::OTHER;           // Internal dedup only; not exposed in output
};

// True when the CIGAR contains any insertion or deletion run.
bool cigar_has_indel(const std::string& cigar);

// Scoring information for a single hit.
// cfd_score is the per-mismatch SpCas9 CFD (Doench et al. 2016); used downstream
// by the spacer summary layer to compute per-spacer aggregate CFD.
struct ScoringInfo {
    double cfd_score = 0.0;  // Per-mismatch CFD (SpCas9; 0.0–1.0)
};

// A single hit from the search.
// Represents a position in the genome where the pattern aligns with edit
// distance <= threshold.
struct SearchHit {
    size_t      genome_pos = 0;     // Absolute position in the concatenated genome
    uint8_t     distance = 0;       // Edit distance at this position
    std::string chrom_name;         // Chromosome name
    size_t      chrom_offset = 0;   // Position within the chromosome (0-based)
    Strand      strand = Strand::PLUS;  // Strand orientation
    MismatchInfo mismatch_info;     // Detailed mismatch data (populated in post-processing)
    ScoringInfo  scoring_info;      // CFD scoring data (populated if compute_scores=true)

    // Default comparison for sorting (by genome position, then strand)
    bool operator<(const SearchHit& other) const {
        if (genome_pos != other.genome_pos) return genome_pos < other.genome_pos;
        return strand < other.strand;
    }
};

// Configuration for a genome search.
struct SearchConfig {
    std::string pattern;                          // The pattern to search for (nucleotide string)
    uint8_t     threshold = 3;                    // Maximum edit distance to report (default: 3 for CRISPR)
    DistanceMode distance_mode = DistanceMode::LEVENSHTEIN;  // Distance metric: Levenshtein (default) or Hamming
    bool        prefer_gpu = true;                // Use GPU if available (falls back to CPU)
    bool        search_both_strands = true;       // Search forward and reverse complement
    bool        compute_mismatches  = true;       // Extract detailed mismatch info (CIGAR, aligned seq, PAM)
    bool        compute_scores      = true;       // Compute CFD activity scores (requires mismatches)
    PamSpec     pam;                              // Pluggable PAM (empty = no filter); see PamSpec
    size_t      max_hits = 0;                     // Maximum hits to report (0 = unlimited)
    bool        forward_only = false;             // When true, skip reverse-complement search entirely
    bool        reverse_only = false;             // When true, skip forward search entirely (only RC pass runs)
    bool        disable_deduplication = false;    // DEBUG: Skip halo deduplication (default: false)

    // Optional search window in absolute (concatenated) genome coordinates [search_start, search_end).
    // If search_end == 0, the whole genome is searched (default).
    size_t      search_start = 0;
    size_t      search_end   = 0;
};

// Result of a genome search.
struct SearchResult {
    std::string           pattern;         // The pattern that was searched
    uint8_t               threshold;       // The threshold used
    std::vector<SearchHit> hits;           // All positions meeting the threshold
    size_t                total_positions; // Total genome positions scanned
    bool                  used_gpu;        // Whether GPU was used for the search
    double                scoring_time_ms = 0.0; // Time spent computing CFD scores
};

// Input normalization

// Replace U/u with T/t in a nucleotide sequence (RNA → DNA normalization).
// Allows users to provide RNA spacer sequences with uracil.
std::string normalize_uracil(const std::string& sequence);

// Filtering functions

// Filter a distance array and return positions where distance <= threshold.
// Returns a vector of positions (indices into the distances array).
std::vector<size_t> filter_by_threshold(const std::vector<uint8_t>& distances,
                                        uint8_t threshold);

// Chromosome resolution

// Given a list of chromosomes (with start positions) and an absolute genome
// position, return the chromosome name and offset within that chromosome.
// Throws std::out_of_range if pos is beyond the last chromosome.
std::tuple<std::string, size_t> resolve_chromosome(
    const std::vector<Chromosome>& chromosomes,
    size_t pos);

// Hit deduplication

// Deduplicate overlapping semi-global alignment hits (Levenshtein/Myers only).
// Myers' algorithm reports edit distance at every genome position, creating
// "halos" of low-distance positions around each real match. This function
// clusters consecutive hits on the same strand (within dedup_radius positions
// of each other) and keeps only the minimum-distance hit per cluster.
// Correct radius for Levenshtein: threshold+1 (halo extends threshold positions).
// For Hamming mode, skip this call entirely — each position is independent.
std::vector<SearchHit> deduplicate_hits(std::vector<SearchHit> hits,
                                         size_t dedup_radius);

// Search pipeline

// Run a full search of the pattern against the genome.
// 1. Validates input (pattern non-empty, valid characters)
// 2. Runs Myers algorithm (GPU if available and preferred, else CPU)
// 3. Filters results by threshold
// 4. Resolves chromosome names/offsets for each hit
// 5. Returns sorted hits
//
// Throws:
//   std::invalid_argument - if pattern is empty or contains invalid characters
SearchResult search_genome(const SearchConfig& config, GenomeView view);

// Convenience overload for Genome (creates view internally).
inline SearchResult search_genome(const SearchConfig& config, const Genome& genome) {
    return search_genome(config, make_view(genome));
}

// Output formatting

// Format hits as tab-separated values (TSV).
// Header: chrom  start  end  pattern  distance
// Positions are 0-based, end is exclusive.
std::string format_hits_tsv(const std::vector<SearchHit>& hits,
                            const std::string& pattern);

// Format hits as BED format.
// chrom  start  end  name  score  strand
// start is 0-based, end is exclusive, score is (m - distance).
std::string format_hits_bed(const std::vector<SearchHit>& hits,
                            const std::string& pattern);

// Format hits as a JSON array of hit objects (one object per hit).
// Fields mirror the TSV columns; positions are 0-based, end exclusive.
std::string format_hits_json(const std::vector<SearchHit>& hits,
                             const std::string& pattern);

// Batch processing

// A spacer entry from a batch input file.
struct SpacerEntry {
    std::string name;      // Identifier (e.g., "EMX1" or "spacer_1")
    std::string sequence;  // Nucleotide pattern (e.g., "GAGTCCGAGCAGAAGAAGAA")
    std::string source_file;   // Origin file for diagnostics; empty if constructed programmatically
    size_t      source_line = 0;  // 1-based line number in source_file (0 if unknown)
};

// Configuration for batch genome search.
struct BatchSearchConfig {
    std::vector<SpacerEntry> spacers;             // All spacers to search
    uint8_t     threshold = 3;                    // Maximum edit distance
    DistanceMode distance_mode = DistanceMode::LEVENSHTEIN;  // Distance metric: Levenshtein (default) or Hamming
    bool        prefer_gpu = true;                // Use GPU if available
    bool        search_both_strands = true;       // Search both strands
    bool        compute_mismatches  = true;       // Extract detailed mismatch info (CIGAR, aligned seq, PAM)
    bool        compute_scores      = true;       // Compute CFD activity scores
    PamSpec     pam;                              // Pluggable PAM (empty = no filter); see PamSpec
    size_t      max_hits_per_spacer = 0;          // Max hits per spacer (0 = unlimited)
    bool        verbose = false;                  // Print progress to stderr
    size_t      num_threads = 0;                  // Threads for parallel search (0 = auto)
    bool        forward_only = false;             // When true, skip reverse-complement search entirely
    bool        reverse_only = false;             // When true, skip forward search entirely (only RC pass runs)
    bool        disable_deduplication = false;    // DEBUG: Skip halo deduplication (default: false)

    // Optional search window in absolute genome coordinates [search_start, search_end).
    // If search_end == 0, the whole genome is searched (default).
    size_t      search_start = 0;
    size_t      search_end   = 0;
};

// Result for one spacer in batch mode.
struct BatchSpacerResult {
    std::string spacer_name;
    std::string spacer_sequence;
    SearchResult result;
};

// Result of a batch genome search.
struct BatchSearchResult {
    std::vector<BatchSpacerResult> spacer_results;
    size_t total_hits = 0;
    double total_time_ms = 0.0;   // Search + scoring (wall clock inside search_genome_batch)
    double scoring_time_ms = 0.0; // CFD scoring only (subset of total_time_ms)
    bool used_gpu = false;
    size_t spacers_skipped = 0;
};

// Parse a spacer file into a vector of SpacerEntry.
// Supports formats:
//   - Sequence only: "GAGTCCGAGCAGAAGAAGAA" (auto-named spacer_1, spacer_2, ...)
//   - Named: "EMX1\tGAGTCCGAGCAGAAGAAGAA"
//   - Full: "EMX1\tGAGTCCGAGCAGAAGAAGAA\tNGG" (extra columns ignored)
// Lines starting with # are comments; blank lines are ignored.
// Throws std::runtime_error if file cannot be opened.
// Throws std::invalid_argument if a line has invalid sequence.
std::vector<SpacerEntry> parse_spacer_file(const std::string& filepath);

// Run batch search: genome loaded once, all spacers searched.
// Invalid spacers are skipped with a warning to stderr.
BatchSearchResult search_genome_batch(const BatchSearchConfig& config,
                                      GenomeView view);

// Convenience overload for Genome (creates view internally).
inline BatchSearchResult search_genome_batch(const BatchSearchConfig& config,
                                             const Genome& genome) {
    return search_genome_batch(config, make_view(genome));
}

// Format batch results as TSV with spacer column first.
// Header: spacer  chrom  start  end  pattern  distance  strand  ...
std::string format_batch_hits_tsv(const BatchSearchResult& result);

// Format batch results as BED with spacer name in column 7.
// Columns: chrom  start  end  pattern  score  strand  spacer
std::string format_batch_hits_bed(const BatchSearchResult& result);

// Format batch results as a JSON array of hit objects (one object per hit).
// Each object carries the spacer name/sequence plus the same fields as
// format_hits_json.
std::string format_batch_hits_json(const BatchSearchResult& result);

// Summary mode

// Aggregated counts of off-targets by distance threshold.
struct SummaryEntry {
    std::string spacer_name;
    std::string spacer_sequence;
    size_t total_hits = 0;
    std::map<uint8_t, size_t> hits_by_distance;  // distance → count
};

// Result of a summary-mode search.
struct SummaryResult {
    std::vector<SummaryEntry> entries;
    size_t total_spacers_processed = 0;
    size_t spacers_skipped = 0;
    uint8_t threshold = 0;
};

// Generate summary from batch search results.
SummaryResult summarize_batch_results(const BatchSearchResult& batch_result,
                                       uint8_t threshold);

// Generate summary from a single search result.
SummaryResult summarize_single_result(const SearchResult& result,
                                       const std::string& pattern_name);

// Format summary as JSON.
std::string format_summary_json(const SummaryResult& summary);

// Format summary as TSV.
std::string format_summary_tsv(const SummaryResult& summary);
