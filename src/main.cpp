#include <genome_loader.hpp>
#include <search_pipeline.hpp>
#include <spacer_summary.hpp>
#include <annotation_bed.hpp>
#include <gpu_engine.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

// Usage / Help

static void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]\n"
              << "\n"
              << "Stomata - Exhaustive CRISPR off-target search engine\n"
              << "\n"
              << "Required arguments (one of):\n"
              << "  --pattern PATTERN        Nucleotide pattern to search for (ACGTNacgtn)\n"
              << "  --spacer-file FILE       File with spacer sequences (one per line; use '-' for stdin)\n"
              << "  --index-genome FILE      Create .st index from FASTA (no search)\n"
              << "  --genome FILE            Path to genome FASTA or .st index file\n"
              << "\n"
              << "Optional arguments:\n"
              << "  --region REGION          Restrict search to CHR or CHR:START-END (1-based, inclusive-exclusive)\n"
              << "  --threshold N            Maximum edit distance to report (default: 3)\n"
              << "  --distance-mode MODE     Distance metric: levenshtein (default) or hamming\n"
              << "  --format FORMAT          Output format: tsv (default), bed, or json\n"
              << "  --output FILE            Write results to FILE instead of stdout\n"
              << "  --cpu-only               Force CPU computation (no GPU)\n"
              << "\n"
              << "Strand control:\n"
              << "  --strand MODE            Strand selection (default: both):\n"
              << "                             both   run both passes; report + and -\n"
              << "                             plus   run only the forward pass; only + hits\n"
              << "                             minus  run only the RC pass; only - hits\n"
              << "\n"
              << "PAM filtering (pluggable, cross-technology):\n"
              << "  --pam PATTERN            IUPAC PAM pattern to filter on (e.g. NGG, TTTV). Empty = no filter.\n"
              << "  --pam-position POS       PAM position relative to spacer: 3prime (default, Cas9) or 5prime (Cas12)\n"
              << "  --pam-extract-length N   Width of PAM bases to extract/annotate (default: length of --pam pattern)\n"
              << "  --compute-mismatches     Compute mismatch details (default: true)\n"
              << "  --no-compute-mismatches  Skip mismatch computation\n"
              << "\n"
              << "Activity scoring:\n"
              << "  --compute-scores         Compute CFD activity scores (default: true)\n"
              << "  --no-scores              Disable CFD scoring for faster output\n"
              << "\n"
              << "Spacer summary (per-spacer aggregate report):\n"
              << "  --spacer-summary FILE    Write per-spacer summary TSV to FILE\n"
              << "  --acfd-threshold FLOAT   Promiscuity threshold for aggregate CFD (default: 4.8)\n"
              << "  --intersect-bed LABEL:FILE  Count off-target overlaps with a BED file; repeatable.\n"
              << "                           LABEL appears as n_LABEL_overlaps in summary output.\n"
              << "                           Example: --intersect-bed CDS:/path/cds.bed\n"
              << "\n"
              << "Performance options:\n"
              << "  --threads N, -t N        Number of threads for batch processing (0 = auto)\n"
              << "\n"
              << "Input normalization:\n"
              << "  --treat-u-as-t           Treat U (uracil) as T (thymine) in spacers (default: on)\n"
              << "  --no-treat-u-as-t        Disable U→T normalization\n"
              << "\n"
              << "Output control:\n"
              << "  --max-hits N             Limit output to first N hits per spacer (0 = unlimited)\n"
              << "  --max-total-hits N       Batch mode: cap total hits across all spacers (0 = unlimited)\n"
              << "  --summary                Output aggregated counts by distance instead of per-hit details\n"
              << "  --summary-format FORMAT  Summary output format: json (default) or tsv\n"
              << "  --verbose                Show detailed statistics on stderr\n"
              << "  --quiet                  Suppress all non-error output on stderr\n"
              << "\n"
              << "Debug options:\n"
              << "  --no-deduplicate         Disable halo deduplication (emits every Myers end-position hit)\n"
              << "\n"
              << "Information:\n"
              << "  --help                   Show this help message and exit\n"
              << "  --version                Show version information and exit\n"
              << "  --quickstart             Run a self-contained smoke test (no genome/spacer needed)\n"
              << "                           and report whether Stomata is functional on this system.\n"
              << "\n"
              << "Spacer file format:\n"
              << "  # Lines starting with # are comments\n"
              << "  GAGTCCGAGCAGAAGAAGAA           # Just sequence (auto-named spacer_1, ...)\n"
              << "  EMX1    GAGTCCGAGCAGAAGAAGAA   # Name<TAB>Sequence\n"
              << "  FANCF   GGAATCCCTTCTGCAGCACC   NGG   # Additional columns ignored\n"
              << "\n"
              << "Examples:\n"
              << "  " << program_name << " --pattern ACGTACGTACGTACGTACGT --genome hg38.fa --threshold 3\n"
              << "  " << program_name << " --spacer-file guides.txt --genome hg38.fa --threshold 3\n"
              << "  " << program_name << " --pattern ACGTACGTACGTACGTACGT --genome hg38.fa --pam NGG --max-hits 1000\n"
              << "  " << program_name << " --index-genome hg38.fa                  # Creates hg38.fa.st\n"
              << "  " << program_name << " --genome hg38.fa.st --pattern ...    # Uses mmap for fast load\n"
              << "\n"
              << "Distance modes:\n"
              << "  levenshtein              Full edit distance (substitutions + indels) - default\n"
              << "  hamming                  Substitution-only distance (no indels) - faster\n"
              << "\n"
              << "Notes:\n"
              << "  - Pattern length must be 1-64 for GPU acceleration\n"
              << "  - Patterns longer than 64 bp will use CPU (multi-word Myers)\n"
              << "  - N in the pattern matches any nucleotide at zero cost\n"
              << "  - --pam accepts IUPAC codes (R, Y, S, W, K, M, B, D, H, V, N); supply one pattern per run\n"
              << "  - PAM filtering requires --compute-mismatches (default on)\n"
              << "  - Batch mode (--spacer-file) loads genome once for all spacers\n"
              << "  - .st index files are memory-mapped for instant loading\n"
              << "  - Hamming mode is faster but does not account for insertions/deletions\n"
              << "\n";
}

static void print_version() {
    std::cout << "stomata version " << STOMATA_VERSION << "\n"
              << "Built with GPU support: " << (gpu_available() ? "yes" : "no") << "\n";
}

// Argument parsing

// Strand selection. Each non-default mode skips one search pass entirely so the
// output naturally contains only hits on the requested strand:
//   BOTH  - run both passes (default)
//   PLUS  - run only the forward pass; all hits are on the + strand
//   MINUS - run only the reverse-complement pass; all hits are on the - strand
enum class StrandMode { BOTH, PLUS, MINUS };

struct Arguments {
    std::string pattern;
    std::string spacer_file;   // Batch mode: file with spacer sequences
    std::string genome_path;
    std::string index_genome;  // Index creation mode: input FASTA path
    std::string output_path;
    std::string format = "tsv";
    std::string pam_pattern;                   // IUPAC PAM pattern (empty = no filter)
    std::string pam_position = "3prime";       // 3prime (Cas9) or 5prime (Cas12)
    size_t      pam_extract_length = 0;        // 0 = use pam_pattern length
    std::string summary_format = "json";
    std::string distance_mode = "levenshtein";  // Distance metric: levenshtein or hamming
    int threshold = 3;
    size_t max_hits = 0;
    size_t max_total_hits = 0;  // Batch mode: global cap across all spacers (0 = unlimited)
    std::string region;  // Optional search window: "chr" or "chr:start-end" (1-based, inclusive-exclusive in output)
    bool cpu_only = false;
    bool verbose = false;
    bool quiet = false;
    bool show_help = false;
    bool show_version = false;
    bool quickstart = false;
    StrandMode strand_mode = StrandMode::BOTH;
    bool compute_mismatches = true;
    bool compute_scores = true;  // CFD activity scoring (default: on)
    bool treat_u_as_t = true;   // Normalize U→T in spacer sequences (default: on)
    bool summary_mode = false;  // Output aggregated counts instead of per-hit details
    bool disable_deduplication = false;  // Skip halo dedup (cross-tool validation).
    size_t num_threads = 0;  // 0 = auto (hardware_concurrency)
    std::string spacer_summary_path;   // If non-empty, write per-spacer summary TSV here
    double acfd_threshold = spacer::DEFAULT_ACFD_THRESHOLD;
    std::vector<std::pair<std::string, std::string>> intersect_beds;  // (label, path)
};

// Table-driven argument parser.
//
// Each OptionDef row binds a long flag (and optional alias) to an `apply`
// lambda.  Value-taking options receive the raw next-argv token; flag options
// receive nullptr.  The lambda returns false after printing an error to stderr
// to abort parsing.

namespace {

struct OptionDef {
    const char* name;
    const char* alias;     // nullptr if none
    bool takes_value;
    std::function<bool(Arguments&, const char*)> apply;
};

// Store a non-negative integer parsed from `v`.  Used for --threshold,
// --max-hits, --max-total-hits, --threads.
bool parse_nonneg_int(const char* flag, const char* v, long long& out) {
    try {
        out = std::stoll(v);
    } catch (const std::exception&) {
        std::cerr << "Error: " << flag << " requires a numeric argument\n";
        return false;
    }
    if (out < 0) {
        std::cerr << "Error: " << flag << " must be non-negative\n";
        return false;
    }
    return true;
}

std::string to_lower_copy(const char* v) {
    std::string s = v;
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Whitelist-check a string value; print error + allowed list on failure.
bool require_one_of(const char* flag, const std::string& got,
                    std::initializer_list<const char*> allowed,
                    const char* allowed_display) {
    for (const char* a : allowed) if (got == a) return true;
    std::cerr << "Error: " << flag << " must be " << allowed_display << "\n";
    return false;
}

const std::vector<OptionDef>& option_table() {
    static const std::vector<OptionDef> kOptions = {
        {"--help", "-h", false, [](Arguments& a, const char*){ a.show_help = true; return true; }},
        {"--version", "-v", false, [](Arguments& a, const char*){ a.show_version = true; return true; }},
        {"--quickstart", nullptr, false, [](Arguments& a, const char*){ a.quickstart = true; return true; }},

        {"--pattern",      nullptr, true, [](Arguments& a, const char* v){ a.pattern      = v; return true; }},
        {"--spacer-file",  nullptr, true, [](Arguments& a, const char* v){ a.spacer_file  = v; return true; }},
        {"--genome",       nullptr, true, [](Arguments& a, const char* v){ a.genome_path  = v; return true; }},
        {"--index-genome", nullptr, true, [](Arguments& a, const char* v){ a.index_genome = v; return true; }},
        {"--output",       nullptr, true, [](Arguments& a, const char* v){ a.output_path  = v; return true; }},
        {"--region",       nullptr, true, [](Arguments& a, const char* v){ a.region       = v; return true; }},

        {"--threshold", nullptr, true, [](Arguments& a, const char* v){
            long long n; if (!parse_nonneg_int("--threshold", v, n)) return false;
            a.threshold = static_cast<int>(n); return true;
        }},
        {"--max-hits", nullptr, true, [](Arguments& a, const char* v){
            long long n; if (!parse_nonneg_int("--max-hits", v, n)) return false;
            a.max_hits = static_cast<size_t>(n); return true;
        }},
        {"--max-total-hits", nullptr, true, [](Arguments& a, const char* v){
            long long n; if (!parse_nonneg_int("--max-total-hits", v, n)) return false;
            a.max_total_hits = static_cast<size_t>(n); return true;
        }},
        {"--threads", "-t", true, [](Arguments& a, const char* v){
            long long n; if (!parse_nonneg_int("--threads", v, n)) return false;
            a.num_threads = static_cast<size_t>(n); return true;
        }},

        {"--format", nullptr, true, [](Arguments& a, const char* v){
            a.format = v;
            return require_one_of("--format", a.format, {"tsv", "bed", "json"}, "'tsv', 'bed', or 'json'");
        }},
        {"--summary-format", nullptr, true, [](Arguments& a, const char* v){
            a.summary_format = v;
            return require_one_of("--summary-format", a.summary_format, {"json", "tsv"}, "'json' or 'tsv'");
        }},
        {"--distance-mode", nullptr, true, [](Arguments& a, const char* v){
            a.distance_mode = to_lower_copy(v);
            return require_one_of("--distance-mode", a.distance_mode, {"levenshtein", "hamming"}, "'levenshtein' or 'hamming'");
        }},
        {"--pam", nullptr, true, [](Arguments& a, const char* v){
            a.pam_pattern = v;
            for (char& c : a.pam_pattern) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
            for (char c : a.pam_pattern) {
                switch (c) {
                    case 'A': case 'C': case 'G': case 'T': case 'U':
                    case 'R': case 'Y': case 'S': case 'W': case 'K':
                    case 'M': case 'B': case 'D': case 'H': case 'V':
                    case 'N':
                        break;
                    default:
                        std::cerr << "Error: --pam has invalid IUPAC character '" << c << "'\n";
                        return false;
                }
            }
            return true;
        }},
        {"--pam-position", nullptr, true, [](Arguments& a, const char* v){
            a.pam_position = to_lower_copy(v);
            return require_one_of("--pam-position", a.pam_position,
                                  {"3prime", "5prime"}, "'3prime' or '5prime'");
        }},
        {"--pam-extract-length", nullptr, true, [](Arguments& a, const char* v){
            long long n; if (!parse_nonneg_int("--pam-extract-length", v, n)) return false;
            a.pam_extract_length = static_cast<size_t>(n);
            return true;
        }},
        {"--strand", nullptr, true, [](Arguments& a, const char* v){
            std::string s = to_lower_copy(v);
            if (s == "both")       a.strand_mode = StrandMode::BOTH;
            else if (s == "plus")  a.strand_mode = StrandMode::PLUS;
            else if (s == "minus") a.strand_mode = StrandMode::MINUS;
            else {
                std::cerr << "Error: --strand must be 'both', 'plus', or 'minus'\n";
                return false;
            }
            return true;
        }},

        {"--cpu-only", nullptr, false, [](Arguments& a, const char*){ a.cpu_only = true; return true; }},
        {"--verbose",  nullptr, false, [](Arguments& a, const char*){ a.verbose  = true; return true; }},
        {"--quiet",    nullptr, false, [](Arguments& a, const char*){ a.quiet    = true; return true; }},
        {"--summary",  nullptr, false, [](Arguments& a, const char*){ a.summary_mode = true; return true; }},

        {"--compute-mismatches",    nullptr, false, [](Arguments& a, const char*){ a.compute_mismatches = true;  return true; }},
        {"--no-compute-mismatches", nullptr, false, [](Arguments& a, const char*){ a.compute_mismatches = false; return true; }},
        {"--compute-scores",        nullptr, false, [](Arguments& a, const char*){ a.compute_scores = true;  return true; }},
        {"--no-scores",             nullptr, false, [](Arguments& a, const char*){ a.compute_scores = false; return true; }},
        {"--treat-u-as-t",          nullptr, false, [](Arguments& a, const char*){ a.treat_u_as_t = true;  return true; }},
        {"--no-treat-u-as-t",       nullptr, false, [](Arguments& a, const char*){ a.treat_u_as_t = false; return true; }},
        {"--no-deduplicate",        nullptr, false, [](Arguments& a, const char*){ a.disable_deduplication = true; return true; }},

        {"--spacer-summary", nullptr, true,  [](Arguments& a, const char* v){ a.spacer_summary_path = v; return true; }},
        {"--intersect-bed",  nullptr, true,  [](Arguments& a, const char* v) {
            std::string s = v;
            auto colon = s.find(':');
            if (colon == std::string::npos || colon == 0 || colon + 1 == s.size()) {
                std::cerr << "Error: --intersect-bed requires LABEL:FILE format\n";
                return false;
            }
            a.intersect_beds.emplace_back(s.substr(0, colon), s.substr(colon + 1));
            return true;
        }},
        {"--acfd-threshold", nullptr, true,  [](Arguments& a, const char* v){
            try { a.acfd_threshold = std::stod(v); }
            catch (const std::exception&) {
                std::cerr << "Error: --acfd-threshold requires a numeric argument\n";
                return false;
            }
            if (a.acfd_threshold < 0.0) {
                std::cerr << "Error: --acfd-threshold must be non-negative\n";
                return false;
            }
            return true;
        }},
    };
    return kOptions;
}

const OptionDef* find_option(const std::string& arg) {
    for (const auto& opt : option_table()) {
        if (arg == opt.name || (opt.alias && arg == opt.alias)) return &opt;
    }
    return nullptr;
}

}  // namespace

static bool parse_arguments(int argc, char** argv, Arguments& args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        const OptionDef* opt = find_option(arg);
        if (!opt) {
            std::cerr << "Error: Unknown argument: " << arg << "\n";
            return false;
        }

        const char* value = nullptr;
        if (opt->takes_value) {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << opt->name << " requires an argument\n";
                return false;
            }
            value = argv[++i];
        }

        if (!opt->apply(args, value)) return false;

        // Preserve prior short-circuit: --help / --version stop further
        // parsing (those code paths intentionally ignore everything else).
        // --quickstart is NOT short-circuited so it can compose with
        // --cpu-only / --threads / --verbose etc. for portability checks.
        if (args.show_help || args.show_version) return true;
    }
    return true;
}

static bool validate_arguments(const Arguments& args) {
    // Index creation mode: only --index-genome is required
    bool index_mode = !args.index_genome.empty();
    if (index_mode) {
        // In index mode, we don't need --pattern, --spacer-file, or --genome
        if (!args.pattern.empty() || !args.spacer_file.empty()) {
            std::cerr << "Error: --index-genome cannot be combined with --pattern or --spacer-file\n";
            return false;
        }
        if (!args.genome_path.empty()) {
            std::cerr << "Error: --index-genome cannot be combined with --genome\n";
            return false;
        }
        return true;
    }

    // Check that either --pattern or --spacer-file is provided (but not both)
    bool has_pattern = !args.pattern.empty();
    bool has_spacer_file = !args.spacer_file.empty();

    if (!has_pattern && !has_spacer_file) {
        std::cerr << "Error: --pattern or --spacer-file is required\n";
        return false;
    }

    if (has_pattern && has_spacer_file) {
        std::cerr << "Error: --pattern and --spacer-file are mutually exclusive\n";
        return false;
    }

    if (args.genome_path.empty()) {
        std::cerr << "Error: --genome is required\n";
        return false;
    }

    // Validate pattern characters (only if using --pattern)
    if (has_pattern) {
        for (char c : args.pattern) {
            switch (c) {
                case 'A': case 'a':
                case 'C': case 'c':
                case 'G': case 'g':
                case 'T': case 't':
                case 'N': case 'n':
                case 'U': case 'u':  // Accepted when treat_u_as_t is enabled
                    break;
                default:
                    std::cerr << "Error: Invalid character in pattern: '" << c << "'\n";
                    return false;
            }
        }
        // Check U usage when normalization is disabled
        if (!args.treat_u_as_t) {
            for (char c : args.pattern) {
                if (c == 'U' || c == 'u') {
                    std::cerr << "Error: Pattern contains 'U' but --no-treat-u-as-t is set. "
                              << "Use T instead of U, or enable U→T normalization.\n";
                    return false;
                }
            }
        }
    }

    // --pam + --no-compute-mismatches:
    //   Hamming mode — allowed. Hamming has no bulges; the canonical PAM is
    //     authoritative and the early filter is exact.
    //   Levenshtein mode — rejected. The widened early filter accepts hits
    //     whose PAM matches at any candidate target length; the principled
    //     re-filter that drops false positives (alignment-validity check)
    //     runs inside the compute_mismatches block. Without compute_mismatches
    //     the false positives are returned (see v0.8.0 fix).
    if (!args.pam_pattern.empty() && !args.compute_mismatches &&
        args.distance_mode == "levenshtein") {
        std::cerr << "Error: --pam requires --compute-mismatches in Levenshtein "
                     "mode (the alignment-validity re-filter cannot run otherwise). "
                     "Use --distance-mode hamming if you need --no-compute-mismatches "
                     "with --pam.\n";
        return false;
    }

    if ((args.summary_mode || !args.spacer_summary_path.empty()) &&
        (args.max_hits > 0 || args.max_total_hits > 0)) {
        std::cerr << "Error: hit caps cannot be combined with --summary or "
                     "--spacer-summary: truncated hits would undercount totals "
                     "and activity scores. Remove the caps for complete summaries.\n";
        return false;
    }

    return true;
}

// Subcommand handlers

// Self-contained smoke test for hardware portability checks. Bundles a tiny
// synthetic FASTA in the binary, runs Stomata on it under three configurations
// (Lev+NGG, Lev no-PAM, Hamming+NGG), and verifies hit counts/distances
// match the engineered ground truth. Exit code 0 on pass, non-zero on fail.
//
// Designed to be one-liner: `stomata --quickstart`.
//
// What it exercises end-to-end:
//   - FASTA parsing (load_fasta)
//   - GPU search dispatch (Myers and shift-add kernels)
//   - PAM filter (NGG)
//   - Lev's alignment-validity re-filter (v0.8.0 fix)
//   - Mismatch info / CIGAR construction
//   - CFD scoring
//
// What it does NOT exercise: index format (.st), batch mode, very long
// patterns (uint64 kernel branch), large genomes. Those are covered by
// the test suite (`ctest`).
static int run_quickstart(const Arguments& args) {
    // Synthetic genome: 3 deliberately placed sites for the test pattern
    // GAGTCCGAGCAGAAGAAGAA.
    //   Site 1 (offset  0): exact match,  PAM = AGG (NGG ✓)
    //   Site 2 (offset 50): 1 mismatch,   PAM = TGG (NGG ✓)
    //   Site 3 (offset 100): exact match, PAM = AAA (NGG ✗)
    // Expected:
    //   Lev+NGG / Hamming+NGG: 2 hits (sites 1, 2)
    //   Lev or Hamming, no PAM: 3 hits (sites 1, 2, 3)
    static const char* kPattern = "GAGTCCGAGCAGAAGAAGAA";
    static const char* kFasta =
        ">chr_quickstart synthetic-self-test\n"
        // Site 1: exact + AGG (NGG)        : positions 0..22
        "GAGTCCGAGCAGAAGAAGAAAGG"
        // padding                          : positions 23..49 (27 bp)
        "TTTTTTAAAAACCCCCAAAAAGTGTTC"
        // Site 2: 1 mm at pos 12 + TGG     : positions 50..72
        "AGAGTCCGAGCAGTAGAAGAATGG"
        // padding                          : positions 73..99 (26 bp)
        "AAAAAAGGGGGCCCCCAAAAAACCAA"
        // Site 3: exact + AAA (not NGG)    : positions 100..122
        "AGAGTCCGAGCAGAAGAAGAAAAA"
        // tail padding                     : positions 123..159
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    // Write to a tmp file (load_fasta needs a path).
    char tmpdir[] = "/tmp/stomata-quickstart-XXXXXX";
    if (!::mkdtemp(tmpdir)) {
        std::cerr << "FAIL: could not create tmpdir for quickstart\n";
        return 1;
    }
    const std::string fasta_path = std::string(tmpdir) + "/quickstart.fa";
    {
        std::ofstream f(fasta_path);
        if (!f) {
            std::cerr << "FAIL: could not write quickstart fasta to " << fasta_path << "\n";
            return 1;
        }
        f << kFasta;
    }

    auto cleanup = [&]() {
        std::remove(fasta_path.c_str());
        ::rmdir(tmpdir);
    };

    std::cout << "Stomata quickstart\n";
    std::cout << "  version:     " << STOMATA_VERSION << "\n";
    std::cout << "  GPU built-in: " << (gpu_available() ? "yes" : "no") << "\n\n";

    Genome genome;
    try {
        genome = load_fasta(fasta_path);
    } catch (const std::exception& e) {
        std::cerr << "FAIL [load_fasta]: " << e.what() << "\n";
        cleanup();
        return 1;
    }

    auto run_one = [&](const std::string& label,
                       DistanceMode mode,
                       const std::string& pam,
                       size_t expected_hits) -> bool {
        SearchConfig config;
        config.pattern = kPattern;
        config.threshold = 3;
        config.distance_mode = mode;
        if (!pam.empty()) {
            config.pam.pattern = pam;
            config.pam.position = PamPosition::THREE_PRIME;
            config.pam.extract_length = pam.size();
        }
        config.compute_mismatches = true;
        config.compute_scores = true;
        config.search_both_strands = true;
        config.prefer_gpu = !args.cpu_only;

        auto t0 = std::chrono::high_resolution_clock::now();
        SearchResult result;
        try {
            result = search_genome(config, make_view(genome));
        } catch (const std::exception& e) {
            std::cerr << "FAIL [" << label << "]: " << e.what() << "\n";
            return false;
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        const bool ok = (result.hits.size() == expected_hits);
        std::cout << "  " << (ok ? "[ok]   " : "[FAIL] ")
                  << label << ": "
                  << result.hits.size() << " hits "
                  << "(expected " << expected_hits << "), "
                  << std::fixed << std::setprecision(1) << ms << " ms"
                  << (result.used_gpu ? " [GPU]" : " [CPU]")
                  << "\n";

        if (!ok) {
            for (size_t i = 0; i < result.hits.size(); ++i) {
                const auto& h = result.hits[i];
                std::cout << "      hit " << i << ": "
                          << h.chrom_name << ":"
                          << h.chrom_offset << " "
                          << static_cast<char>(h.strand)
                          << " d=" << static_cast<int>(h.distance)
                          << " pam=" << (h.mismatch_info.pam_sequence.empty()
                                          ? "." : h.mismatch_info.pam_sequence)
                          << "\n";
            }
        }
        return ok;
    };

    bool all_ok = true;
    all_ok &= run_one("Hamming + NGG    ", DistanceMode::HAMMING,     "NGG", 2);
    all_ok &= run_one("Hamming PAM-agn  ", DistanceMode::HAMMING,     "",    3);
    all_ok &= run_one("Lev + NGG        ", DistanceMode::LEVENSHTEIN, "NGG", 2);
    all_ok &= run_one("Lev PAM-agnostic ", DistanceMode::LEVENSHTEIN, "",    3);

    cleanup();

    std::cout << "\n";
    if (all_ok) {
        std::cout << "Quickstart PASS — Stomata is functional on this system.\n";
        return 0;
    } else {
        std::cout << "Quickstart FAIL — see hit details above. "
                     "Likely causes: GPU compute capability not in the binary's "
                     "SASS targets, or a broken CUDA install. See INSTALL.md "
                     "Troubleshooting.\n";
        return 1;
    }
}

static int run_index_mode(const Arguments& args) {
    if (!args.quiet) {
        std::cerr << "Creating index from: " << args.index_genome << "\n";
    }

    auto load_start = std::chrono::high_resolution_clock::now();
    Genome genome = load_fasta(args.index_genome);
    auto load_end = std::chrono::high_resolution_clock::now();

    if (args.verbose) {
        auto load_ms = std::chrono::duration_cast<std::chrono::milliseconds>(load_end - load_start).count();
        std::cerr << "Parsed " << genome.total_bases << " bases from "
                  << genome.chromosomes.size() << " chromosome(s) in "
                  << load_ms << " ms\n";
    }

    std::string index_path = args.output_path.empty()
                             ? args.index_genome + ".st"
                             : args.output_path;

    auto write_start = std::chrono::high_resolution_clock::now();
    write_genome_index(genome, index_path);
    auto write_end = std::chrono::high_resolution_clock::now();

    if (args.verbose) {
        auto write_ms = std::chrono::duration_cast<std::chrono::milliseconds>(write_end - write_start).count();
        std::cerr << "Wrote index to: " << index_path << " in " << write_ms << " ms\n";
    } else if (!args.quiet) {
        std::cerr << "Index created: " << index_path << "\n";
    }

    return 0;
}

// Resolve optional --region to absolute genome coordinates.
// Accepts "chr" or "chr:start-end" (1-based, inclusive-exclusive end).
// On success, sets region_start / region_end (both zero means "whole genome").
// On parse or lookup failure, writes to stderr and returns false.
static bool resolve_region(const Arguments& args, const GenomeView& view,
                           size_t& region_start, size_t& region_end) {
    region_start = 0;
    region_end   = 0;
    if (args.region.empty()) return true;

    std::string chrom_name;
    size_t range_start_1b = 1;
    size_t range_end_1b = 0;
    size_t colon = args.region.find(':');
    if (colon == std::string::npos) {
        chrom_name = args.region;
    } else {
        chrom_name = args.region.substr(0, colon);
        std::string range_spec = args.region.substr(colon + 1);
        size_t dash = range_spec.find('-');
        if (dash == std::string::npos) {
            std::cerr << "Error: --region must be CHR or CHR:START-END (got '" << args.region << "')\n";
            return false;
        }
        try {
            range_start_1b = std::stoull(range_spec.substr(0, dash));
            range_end_1b   = std::stoull(range_spec.substr(dash + 1));
        } catch (const std::exception&) {
            std::cerr << "Error: --region range must be numeric (got '" << range_spec << "')\n";
            return false;
        }
        if (range_start_1b < 1 || range_end_1b <= range_start_1b) {
            std::cerr << "Error: --region range must be 1-based with END > START (got " << range_spec << ")\n";
            return false;
        }
    }

    const Chromosome* chrom = nullptr;
    for (const auto& c : *view.chromosomes) {
        if (c.name == chrom_name) { chrom = &c; break; }
    }
    if (!chrom) {
        std::cerr << "Error: --region chromosome '" << chrom_name << "' not found in genome.\n";
        std::cerr << "Available: ";
        for (const auto& c : *view.chromosomes) std::cerr << c.name << " ";
        std::cerr << "\n";
        return false;
    }

    if (range_end_1b == 0) {
        range_start_1b = 1;
        range_end_1b = chrom->length + 1;
    }
    if (range_end_1b > chrom->length + 1) {
        std::cerr << "Error: --region end " << range_end_1b << " exceeds " << chrom->name
                  << " length " << chrom->length << "\n";
        return false;
    }

    region_start = chrom->start + (range_start_1b - 1);
    region_end   = chrom->start + (range_end_1b - 1);

    if (args.verbose) {
        std::cerr << "Region: " << chrom_name << ":" << range_start_1b << "-" << range_end_1b
                  << " (absolute " << region_start << "-" << region_end << ")\n";
    }
    return true;
}

static PamSpec build_pam_spec(const Arguments& args) {
    PamSpec spec;
    spec.pattern = args.pam_pattern;
    spec.position = (args.pam_position == "5prime")
                    ? PamPosition::FIVE_PRIME
                    : PamPosition::THREE_PRIME;
    spec.extract_length = args.pam_extract_length;
    return spec;
}

static DistanceMode distance_mode_from_string(const std::string& s) {
    return (s == "hamming") ? DistanceMode::HAMMING : DistanceMode::LEVENSHTEIN;
}

// Map StrandMode to the search-pipeline gates. PLUS skips the RC pass,
// MINUS skips the forward pass; BOTH runs both. The output then naturally
// contains only hits on the requested strand(s) — no post-filter needed.
static bool forward_only_for(StrandMode mode) { return mode == StrandMode::PLUS; }
static bool reverse_only_for(StrandMode mode) { return mode == StrandMode::MINUS; }

static std::string run_batch_search(const Arguments& args, GenomeView view,
                                    size_t region_start, size_t region_end) {
    if (args.verbose) {
        std::cerr << "Parsing spacer file: " << args.spacer_file << "\n";
    }

    auto spacers = parse_spacer_file(args.spacer_file);

    if (args.treat_u_as_t) {
        for (auto& spacer : spacers) {
            spacer.sequence = normalize_uracil(spacer.sequence);
        }
    }

    if (args.verbose) {
        std::cerr << "Found " << spacers.size() << " spacer(s)\n";
    }

    BatchSearchConfig batch_config;
    batch_config.spacers = std::move(spacers);
    batch_config.threshold = static_cast<uint8_t>(args.threshold);
    batch_config.prefer_gpu = !args.cpu_only;
    batch_config.compute_mismatches = args.compute_mismatches;
    batch_config.compute_scores = args.compute_scores;
    batch_config.max_hits_per_spacer = args.max_hits;
    batch_config.forward_only = forward_only_for(args.strand_mode);
    batch_config.reverse_only = reverse_only_for(args.strand_mode);
    batch_config.verbose = args.verbose;
    batch_config.num_threads = args.num_threads;
    batch_config.disable_deduplication = args.disable_deduplication;
    batch_config.search_start = region_start;
    batch_config.search_end   = region_end;
    batch_config.distance_mode = distance_mode_from_string(args.distance_mode);
    batch_config.pam = build_pam_spec(args);

    BatchSearchResult batch_result = search_genome_batch(batch_config, view);

    // Apply global --max-total-hits cap (0 = unlimited). Hits are kept in
    // per-spacer order; we truncate from the tail across spacers.
    if (args.max_total_hits > 0 && batch_result.total_hits > args.max_total_hits) {
        size_t remaining = args.max_total_hits;
        size_t dropped = batch_result.total_hits - args.max_total_hits;
        for (auto& sr : batch_result.spacer_results) {
            if (remaining == 0) {
                sr.result.hits.clear();
            } else if (sr.result.hits.size() > remaining) {
                sr.result.hits.resize(remaining);
                remaining = 0;
            } else {
                remaining -= sr.result.hits.size();
            }
        }
        batch_result.total_hits = args.max_total_hits;
        if (!args.quiet) {
            std::cerr << "Note: --max-total-hits cap applied; dropped " << dropped
                      << " hits beyond the first " << args.max_total_hits << ".\n";
        }
    }

    if (args.verbose) {
        double kernel_ms = batch_result.total_time_ms - batch_result.scoring_time_ms;
        std::cerr << "Search (kernel):  " << static_cast<long>(kernel_ms)  << " ms"
                  << "  [GPU: " << (batch_result.used_gpu ? "yes" : "no") << "]\n";
        std::cerr << "Scoring (CFD):    " << static_cast<long>(batch_result.scoring_time_ms) << " ms\n";
        std::cerr << "Total hits: " << batch_result.total_hits << "\n";
    } else if (!args.quiet) {
        std::cerr << "Total hits: " << batch_result.total_hits << "\n";
    }

    // Write per-spacer summary if requested
    if (!args.spacer_summary_path.empty()) {
        // Load each named BED forest; keep alive for the duration of summarization.
        std::vector<std::unique_ptr<annot::BedForest>> forests;
        std::vector<spacer::AnnotationBed> annotations;
        std::vector<std::string> annotation_labels;
        for (const auto& [label, path] : args.intersect_beds) {
            forests.push_back(std::make_unique<annot::BedForest>(annot::load_bed(path)));
            annotations.push_back({label, forests.back().get()});
            annotation_labels.push_back(label);
        }

        spacer::SpacerSummaryConfig sum_config;
        sum_config.acfd_threshold = args.acfd_threshold;

        std::vector<spacer::SpacerSummaryRow> rows;
        rows.reserve(batch_result.spacer_results.size());
        for (const auto& sr : batch_result.spacer_results) {
            rows.push_back(spacer::compute_spacer_summary(
                sr.spacer_name, sr.spacer_sequence, sr.result.hits,
                sum_config, annotations));
        }

        std::ofstream sum_out(args.spacer_summary_path);
        if (!sum_out) {
            std::cerr << "Error: Could not open spacer summary file: "
                      << args.spacer_summary_path << "\n";
        } else {
            sum_out << spacer::format_spacer_summary_tsv(rows, annotation_labels);
            if (args.verbose)
                std::cerr << "Spacer summary written to: "
                          << args.spacer_summary_path << "\n";
        }
    }

    auto fmt_t0 = std::chrono::high_resolution_clock::now();
    std::string formatted;
    if (args.summary_mode) {
        auto summary = summarize_batch_results(
            batch_result, static_cast<uint8_t>(args.threshold));
        formatted = (args.summary_format == "tsv")
            ? format_summary_tsv(summary)
            : format_summary_json(summary);
    } else if (args.format == "tsv")  { formatted = format_batch_hits_tsv(batch_result); }
    else if (args.format == "json")   { formatted = format_batch_hits_json(batch_result); }
    else                              { formatted = format_batch_hits_bed(batch_result); }

    if (args.verbose) {
        auto fmt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - fmt_t0).count();
        std::cerr << "Format (TSV):     " << fmt_ms << " ms"
                  << "  [" << (formatted.size() / 1024 / 1024) << " MB]\n";
    }
    return formatted;
}

static std::string run_single_search(const Arguments& args, GenomeView view,
                                     size_t region_start, size_t region_end) {
    SearchConfig config;
    config.pattern = args.treat_u_as_t ? normalize_uracil(args.pattern) : args.pattern;
    config.threshold = static_cast<uint8_t>(args.threshold);
    config.prefer_gpu = !args.cpu_only;
    config.compute_mismatches = args.compute_mismatches;
    config.compute_scores = args.compute_scores;
    config.forward_only = forward_only_for(args.strand_mode);
    config.reverse_only = reverse_only_for(args.strand_mode);
    config.max_hits = args.max_hits;
    config.disable_deduplication = args.disable_deduplication;
    config.search_start = region_start;
    config.search_end   = region_end;
    config.distance_mode = distance_mode_from_string(args.distance_mode);
    config.pam = build_pam_spec(args);

    if (args.verbose) {
        const char* strand_str = (args.strand_mode == StrandMode::BOTH)  ? "both"
                               : (args.strand_mode == StrandMode::PLUS)  ? "plus"
                               :                                           "minus";
        std::cerr << "Searching for pattern: " << config.pattern << "\n";
        std::cerr << "Pattern length: " << config.pattern.size() << " bp\n";
        std::cerr << "Threshold: " << static_cast<int>(config.threshold) << "\n";
        std::cerr << "Strand: " << strand_str << "\n";
        std::cerr << "Compute mismatches: " << (config.compute_mismatches ? "yes" : "no") << "\n";
        std::cerr << "Compute CFD scores: " << (config.compute_scores ? "yes" : "no") << "\n";
        if (args.pam_pattern.empty()) {
            std::cerr << "PAM filter: (none)\n";
        } else {
            std::cerr << "PAM filter: " << args.pam_pattern
                      << " (" << args.pam_position << ")\n";
        }
        if (config.max_hits > 0) {
            std::cerr << "Max hits: " << config.max_hits << "\n";
        }
        std::cerr << "GPU available: " << (gpu_available() ? "yes" : "no") << "\n";
        std::cerr << "GPU preferred: " << (config.prefer_gpu ? "yes" : "no") << "\n";
    }

    auto search_start = std::chrono::high_resolution_clock::now();
    SearchResult result = search_genome(config, view);
    auto search_end = std::chrono::high_resolution_clock::now();

    if (args.verbose) {
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(search_end - search_start).count();
        long kernel_ms = total_ms - static_cast<long>(result.scoring_time_ms);
        std::cerr << "Search (kernel):  " << kernel_ms << " ms"
                  << "  [GPU: " << (result.used_gpu ? "yes" : "no") << "]\n";
        std::cerr << "Scoring (CFD):    " << static_cast<long>(result.scoring_time_ms) << " ms\n";
        std::cerr << "Positions scanned: " << result.total_positions << "\n";
        std::cerr << "Hits found: " << result.hits.size() << "\n";
    } else if (!args.quiet) {
        std::cerr << "Hits found: " << result.hits.size() << "\n";
    }

    // Write per-spacer summary if requested
    if (!args.spacer_summary_path.empty()) {
        std::vector<std::unique_ptr<annot::BedForest>> forests;
        std::vector<spacer::AnnotationBed> annotations;
        std::vector<std::string> annotation_labels;
        for (const auto& [label, path] : args.intersect_beds) {
            forests.push_back(std::make_unique<annot::BedForest>(annot::load_bed(path)));
            annotations.push_back({label, forests.back().get()});
            annotation_labels.push_back(label);
        }

        spacer::SpacerSummaryConfig sum_config;
        sum_config.acfd_threshold = args.acfd_threshold;

        std::vector<spacer::SpacerSummaryRow> rows;
        rows.push_back(spacer::compute_spacer_summary(
            config.pattern, config.pattern, result.hits, sum_config, annotations));

        std::ofstream sum_out(args.spacer_summary_path);
        if (!sum_out) {
            std::cerr << "Error: Could not open spacer summary file: "
                      << args.spacer_summary_path << "\n";
        } else {
            sum_out << spacer::format_spacer_summary_tsv(rows, annotation_labels);
            if (args.verbose)
                std::cerr << "Spacer summary written to: "
                          << args.spacer_summary_path << "\n";
        }
    }

    auto fmt_t0 = std::chrono::high_resolution_clock::now();
    std::string formatted;
    if (args.summary_mode) {
        auto summary = summarize_single_result(result, args.pattern);
        formatted = (args.summary_format == "tsv")
            ? format_summary_tsv(summary)
            : format_summary_json(summary);
    } else if (args.format == "tsv")  { formatted = format_hits_tsv(result.hits, config.pattern); }
    else if (args.format == "json")   { formatted = format_hits_json(result.hits, config.pattern); }
    else                              { formatted = format_hits_bed(result.hits, config.pattern); }

    if (args.verbose) {
        auto fmt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - fmt_t0).count();
        std::cerr << "Format (TSV):     " << fmt_ms << " ms"
                  << "  [" << (formatted.size() / 1024 / 1024) << " MB]\n";
    }
    return formatted;
}

// Open the genome (auto-detect .st for memory-mapped loading).
// genome_ptr and mapped_ptr are out-params; exactly one will be populated.
static GenomeView load_genome_auto(const Arguments& args,
                                   std::unique_ptr<Genome>& genome_ptr,
                                   std::unique_ptr<MappedGenome>& mapped_ptr) {
    if (args.verbose) {
        std::cerr << "Loading genome from: " << args.genome_path << "\n";
    }

    auto load_start = std::chrono::high_resolution_clock::now();
    GenomeView view;

    if (is_stomata_index(args.genome_path)) {
        mapped_ptr = std::make_unique<MappedGenome>(load_genome_index(args.genome_path));
        view = make_view(*mapped_ptr);
        if (args.verbose) {
            auto load_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - load_start).count();
            std::cerr << "Memory-mapped " << view.total_bases << " bases from "
                      << view.chromosomes->size() << " chromosome(s) in "
                      << load_ms << " ms\n";
        }
    } else {
        genome_ptr = std::make_unique<Genome>(load_fasta(args.genome_path));
        view = make_view(*genome_ptr);
        if (args.verbose) {
            auto load_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - load_start).count();
            std::cerr << "Loaded " << view.total_bases << " bases from "
                      << view.chromosomes->size() << " chromosome(s) in "
                      << load_ms << " ms\n";
        }
    }
    return view;
}

// Main

int main(int argc, char** argv) {
    Arguments args;

    if (!parse_arguments(argc, argv, args)) {
        std::cerr << "Run '" << argv[0] << " --help' for usage.\n";
        return 1;
    }

    if (args.show_help)    { print_usage(argv[0]); return 0; }
    if (args.show_version) { print_version();      return 0; }
    if (args.quickstart)   { return run_quickstart(args); }

    if (!validate_arguments(args)) {
        std::cerr << "Run '" << argv[0] << " --help' for usage.\n";
        return 1;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        if (!args.index_genome.empty()) {
            return run_index_mode(args);
        }

        std::unique_ptr<Genome> genome_ptr;
        std::unique_ptr<MappedGenome> mapped_ptr;
        GenomeView view = load_genome_auto(args, genome_ptr, mapped_ptr);

        size_t region_start = 0, region_end = 0;
        if (!resolve_region(args, view, region_start, region_end)) {
            return 1;
        }

        std::string output = args.spacer_file.empty()
            ? run_single_search(args, view, region_start, region_end)
            : run_batch_search(args, view, region_start, region_end);

        auto write_t0 = std::chrono::high_resolution_clock::now();
        if (args.output_path.empty()) {
            std::cout << output;
        } else {
            std::ofstream outfile(args.output_path);
            if (!outfile) {
                std::cerr << "Error: Could not open output file: " << args.output_path << "\n";
                return 1;
            }
            outfile << output;
        }

        if (args.verbose) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto write_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - write_t0).count();
            auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
            std::cerr << "Write (output):   " << write_ms << " ms"
                      << (args.output_path.empty() ? "  [stdout]" : ("  [" + args.output_path + "]")) << "\n";
            std::cerr << "Total:            " << total_ms << " ms\n";
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
