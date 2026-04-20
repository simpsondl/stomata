#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <genome_loader.hpp>   // Genome, Chromosome

// Version

static constexpr const char* CYSTIDIA_VERSION = "0.5.0";

// Data structures

// Strand orientation of a hit.
enum class Strand : char {
    PLUS  = '+',
    MINUS = '-'
};

// PAM sequence type classification.
enum class PamType : uint8_t {
    NGG,         // NGG PAM (preferred SpCas9)
    NAG,         // NAG PAM (alternative SpCas9)
    OTHER,       // Other PAM sequence (not NGG/NAG)
    INCOMPLETE   // PAM at genome boundary (< 3 bases)
};

// PAM filter type for CRISPR off-target search.
// Controls which PAM sequences are retained in results.
enum class PamFilter : uint8_t {
    NONE,         // No filtering (report all hits regardless of PAM)
    NGG_ONLY,     // Only NGG PAMs
    NAG_ONLY,     // Only NAG PAMs
    NGG_OR_NAG,   // NGG or NAM PAMs (standard SpCas9, default)
    ANY           // Same as NONE (for clarity in user-facing options)
};

// Distance metric mode for pattern matching.
// Controls which algorithm and distance definition to use.
enum class DistanceMode : uint8_t {
    LEVENSHTEIN,  // Full edit distance (substitutions + insertions + deletions) - default
    HAMMING       // Substitution-only distance (no indels/bulges) - faster, simpler
};

// Type of alignment difference at a position.
// Used to annotate mismatches, insertions, and deletions.
enum class EditType : uint8_t {
    NONE       = 0,  // No difference (match)
    MISMATCH   = 1,  // Substitution (different bases)
    DNA_BULGE  = 2,  // Insertion in DNA/genome (deletion in pattern/RNA)
    RNA_BULGE  = 3   // Insertion in RNA/pattern (deletion in DNA/genome)
};

// Alignment operation for traceback
enum class AlignOp : uint8_t {
    MATCH,      // Match or mismatch (diagonal move)
    INS_TEXT,   // Insertion in text/genome (horizontal move) = DNA bulge
    INS_PATTERN // Insertion in pattern/RNA (vertical move) = RNA bulge
};

// Detailed annotation information for a single hit.
// Provides biologically-meaningful summary counts for CRISPR analysis.
// For a spacer of length m:
//   - DISTAL region: positions 1 to (m-8), distal from PAM
//   - SEED region:   positions (m-7) to m, proximal to PAM (8 nt)
struct MismatchInfo {
    std::string           aligned_sequence;         // Genome sequence at hit (spacer region)
    std::string           pam_sequence;             // 3 bases following spacer
    PamType               pam_type = PamType::OTHER;  // PAM classification (NGG/NAG/OTHER/INCOMPLETE)
    
    // Position-level detail (for advanced analysis)
    std::vector<uint8_t>  mismatch_positions;       // 1-based positions in spacer with edits
    std::vector<EditType> edit_types;               // Type of edit at each mismatch position
    
    // Biologically-meaningful summary counts
    uint8_t seed_mismatches      = 0;  // Count of mismatches in seed region
    uint8_t seed_dna_bulges      = 0;  // Count of DNA bulges in seed region
    uint8_t seed_rna_bulges      = 0;  // Count of RNA bulges in seed region
    uint8_t distal_mismatches    = 0;  // Count of mismatches in distal region
    uint8_t distal_dna_bulges    = 0;  // Count of DNA bulges in distal region
    uint8_t distal_rna_bulges    = 0;  // Count of RNA bulges in distal region
    
    // Alignment metadata
    bool        alignment_is_ambiguous  = false;  // True if multiple optimal alignments exist
    uint16_t    n_ambiguous_cells       = 0;      // Traceback cells where >1 op was co-optimal
    std::string cigar;                             // CIGAR-style string: M (match/mismatch), I (RNA bulge), D (DNA bulge)

    // Helper methods for common queries
    uint8_t total_mismatches() const { return seed_mismatches + distal_mismatches; }
    uint8_t total_dna_bulges() const { return seed_dna_bulges + distal_dna_bulges; }
    uint8_t total_rna_bulges() const { return seed_rna_bulges + distal_rna_bulges; }
    uint8_t total_seed_edits() const { return seed_mismatches + seed_dna_bulges + seed_rna_bulges; }
    uint8_t total_distal_edits() const { return distal_mismatches + distal_dna_bulges + distal_rna_bulges; }
    bool has_dna_bulge() const { return (seed_dna_bulges + distal_dna_bulges) > 0; }
    bool has_rna_bulge() const { return (seed_rna_bulges + distal_rna_bulges) > 0; }
};

// Scoring information for a single hit.
// Contains CFD (Cutting Frequency Determination) score and derived metrics.
// CFD scores predict off-target cleavage activity based on Doench et al. 2016.
struct ScoringInfo {
    double   cfd_score = 0.0;     // CFD score (0.0-1.0), higher = more likely to cleave
    uint8_t  risk_tier = 0;       // Risk classification: 0=low, 1=medium, 2=high

    // Helper method
    bool is_predicted_active(double threshold = 0.1) const {
        return cfd_score >= threshold;
    }
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
    uint8_t     threshold = 4;                    // Maximum edit distance to report (default: 4 for CRISPR)
    DistanceMode distance_mode = DistanceMode::LEVENSHTEIN;  // Distance metric: Levenshtein (default) or Hamming
    bool        prefer_gpu = true;                // Use GPU if available (falls back to CPU)
    bool        search_both_strands = true;       // Search forward and reverse complement
    bool        compute_mismatches  = true;       // Extract detailed mismatch positions
    bool        compute_scores      = true;       // Compute CFD activity scores (requires mismatches)
    PamFilter   pam_filter = PamFilter::NGG_OR_NAG;  // Filter hits by PAM type (default: NGG or NAG)
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

// Deduplicate overlapping semi-global alignment hits.
// Myers' algorithm reports edit distance at every genome position, creating
// "halos" of low-distance positions around each real match. This function
// clusters consecutive hits on the same strand (within pattern_len distance)
// and keeps only the minimum-distance hit per cluster.
std::vector<SearchHit> deduplicate_hits(std::vector<SearchHit> hits,
                                         size_t pattern_len);

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
    uint8_t     threshold = 4;                    // Maximum edit distance
    DistanceMode distance_mode = DistanceMode::LEVENSHTEIN;  // Distance metric: Levenshtein (default) or Hamming
    bool        prefer_gpu = true;                // Use GPU if available
    bool        search_both_strands = true;       // Search both strands
    bool        compute_mismatches  = true;       // Extract detailed mismatch info
    bool        compute_scores      = true;       // Compute CFD activity scores
    bool        compute_mit_score   = false;      // Compute MIT specificity scores (requires compute_scores=true)
    PamFilter   pam_filter = PamFilter::NONE;     // Filter by PAM type
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
    double mit_specificity_score = 0.0;  // MIT specificity score (0-100)
    uint8_t specificity_tier = 0;        // Quality tier: 0=poor, 1=fair, 2=good, 3=excellent
};

// Result of a batch genome search.
struct BatchSearchResult {
    std::vector<BatchSpacerResult> spacer_results;
    size_t total_hits = 0;
    double total_time_ms = 0.0;
    bool used_gpu = false;
    size_t spacers_skipped = 0;  // Count of spacers skipped due to validation errors
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
// format_hits_json, with added mit_specificity/specificity_tier columns.
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
