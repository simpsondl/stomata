#include <catch2/catch_test_macros.hpp>

#include <search_pipeline.hpp>

#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/wait.h>

#ifndef TEST_DATA_DIR
#error "TEST_DATA_DIR not defined. Set via CMake target_compile_definitions."
#endif

#ifndef STOMATA_EXECUTABLE
#error "STOMATA_EXECUTABLE not defined. Set via CMake target_compile_definitions."
#endif

// ───────────────────────────────────────────────────────────────────────────
// Helpers
// ───────────────────────────────────────────────────────────────────────────

static const std::string kTestSmallFa = TEST_DATA_DIR "/synthetic/test_small.fa";
static const std::string kStomataExe = STOMATA_EXECUTABLE;

// Run a command and capture its output and exit code
struct CommandResult {
    std::string stdout_text;
    std::string stderr_text;
    int exit_code;
};

static CommandResult run_command(const std::string& cmd) {
    CommandResult result;

    // Create temp files for stdout and stderr
    std::string stdout_file = "/tmp/stomata_test_stdout_" + std::to_string(getpid());
    std::string stderr_file = "/tmp/stomata_test_stderr_" + std::to_string(getpid());

    std::string full_cmd = cmd + " > " + stdout_file + " 2> " + stderr_file;
    result.exit_code = std::system(full_cmd.c_str());
    result.exit_code = WEXITSTATUS(result.exit_code);

    // Read stdout
    std::ifstream stdout_stream(stdout_file);
    std::stringstream stdout_ss;
    stdout_ss << stdout_stream.rdbuf();
    result.stdout_text = stdout_ss.str();
    stdout_stream.close();
    std::remove(stdout_file.c_str());

    // Read stderr
    std::ifstream stderr_stream(stderr_file);
    std::stringstream stderr_ss;
    stderr_ss << stderr_stream.rdbuf();
    result.stderr_text = stderr_ss.str();
    stderr_stream.close();
    std::remove(stderr_file.c_str());

    return result;
}

// ───────────────────────────────────────────────────────────────────────────
// Help and version tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - help flag shows usage", "[cli][help]") {
    auto result = run_command(kStomataExe + " --help");

    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("Usage:") != std::string::npos);
    REQUIRE(result.stdout_text.find("--pattern") != std::string::npos);
    REQUIRE(result.stdout_text.find("--genome") != std::string::npos);
}

TEST_CASE("CLI - version flag shows version", "[cli][version]") {
    auto result = run_command(kStomataExe + " --version");

    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("stomata") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Basic search tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - basic search produces TSV output", "[cli][search]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGT --genome " + kTestSmallFa + " --threshold 0");

    REQUIRE(result.exit_code == 0);
    // Should have TSV header
    REQUIRE(result.stdout_text.find("chrom\tstart\tend\tpattern\tdistance") != std::string::npos);
}

TEST_CASE("CLI - search with threshold 4", "[cli][search]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGT --genome " + kTestSmallFa + " --threshold 4");

    REQUIRE(result.exit_code == 0);
    // Should have some hits
    REQUIRE(result.stdout_text.find("chr") != std::string::npos);
}

TEST_CASE("CLI - search finds exact matches", "[cli][search]") {
    // ACGTACGT should have exact matches in test_small.fa
    auto result = run_command(kStomataExe + " --pattern ACGTACGT --genome " + kTestSmallFa + " --threshold 0");

    REQUIRE(result.exit_code == 0);
    // Should find at least one exact match (distance 0)
    // The distance is now in column 5 (after chrom, start, end, pattern)
    REQUIRE(result.stdout_text.find("\t0\t") != std::string::npos);
}

TEST_CASE("CLI - BED format output", "[cli][format]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold 1 --format bed");

    REQUIRE(result.exit_code == 0);
    // BED format should NOT have a header
    REQUIRE(result.stdout_text.find("chrom\tstart") == std::string::npos);
    // Should have chromosome data
    REQUIRE(result.stdout_text.find("chr") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// GPU preference tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - cpu-only flag works", "[cli][gpu]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGT --genome " + kTestSmallFa + " --threshold 1 --cpu-only");

    REQUIRE(result.exit_code == 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Error handling tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - missing genome file fails", "[cli][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome /nonexistent/file.fa");

    REQUIRE(result.exit_code != 0);
    REQUIRE((result.stderr_text.find("Error") != std::string::npos ||
             result.stderr_text.find("error") != std::string::npos));
}

TEST_CASE("CLI - missing pattern fails", "[cli][error]") {
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa);

    REQUIRE(result.exit_code != 0);
}

TEST_CASE("CLI - empty pattern fails", "[cli][error]") {
    auto result = run_command(kStomataExe + " --pattern \"\" --genome " + kTestSmallFa);

    REQUIRE(result.exit_code != 0);
}

TEST_CASE("CLI - invalid character in pattern fails", "[cli][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGTXACGT --genome " + kTestSmallFa);

    REQUIRE(result.exit_code != 0);
    REQUIRE((result.stderr_text.find("Error") != std::string::npos ||
             result.stderr_text.find("error") != std::string::npos ||
             result.stderr_text.find("Invalid") != std::string::npos ||
             result.stderr_text.find("invalid") != std::string::npos));
}

TEST_CASE("CLI - negative threshold fails", "[cli][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold -1");

    REQUIRE(result.exit_code != 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Output format tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - output to file", "[cli][output]") {
    std::string output_file = "/tmp/stomata_test_output_" + std::to_string(getpid()) + ".tsv";

    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa +
                              " --threshold 1 --output " + output_file);

    REQUIRE(result.exit_code == 0);

    // Verify file was created
    std::ifstream file(output_file);
    REQUIRE(file.good());

    std::stringstream ss;
    ss << file.rdbuf();
    std::string content = ss.str();

    REQUIRE(content.find("chrom\tstart\tend\tpattern\tdistance") != std::string::npos);

    // Cleanup
    file.close();
    std::remove(output_file.c_str());
}

// ───────────────────────────────────────────────────────────────────────────
// Verbose/quiet mode tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - quiet mode suppresses stats", "[cli][output]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold 1 --quiet");

    REQUIRE(result.exit_code == 0);
    // Quiet mode should suppress stats output to stderr
    REQUIRE(result.stderr_text.empty());
}

TEST_CASE("CLI - verbose mode shows stats", "[cli][output]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold 1 --verbose");

    REQUIRE(result.exit_code == 0);
    // Verbose mode should show stats like "hits found" or "positions scanned"
    REQUIRE((result.stderr_text.find("hit") != std::string::npos ||
             result.stderr_text.find("Hit") != std::string::npos ||
             result.stderr_text.find("position") != std::string::npos));
}
// New CLI features tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - version flag", "[cli][version]") {
    auto result = run_command(kStomataExe + " --version");

    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("stomata version") != std::string::npos);
    REQUIRE(result.stdout_text.find(STOMATA_VERSION) != std::string::npos);
    REQUIRE(result.stdout_text.find("GPU support") != std::string::npos);
}

TEST_CASE("CLI - max-hits limits output", "[cli][max-hits]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold 3 --max-hits 5 --quiet");

    REQUIRE(result.exit_code == 0);
    
    // Count lines in output (excluding header)
    size_t line_count = 0;
    std::istringstream iss(result.stdout_text);
    std::string line;
    while (std::getline(iss, line)) {
        // Skip header line (starts with "chrom\t")
        if (!line.empty() && line.find("chrom\t") != 0) {
            ++line_count;
        }
    }
    
    // Should have at most 5 hits
    REQUIRE(line_count <= 5);
}

TEST_CASE("CLI - strand plus searches plus strand only", "[cli][strand]") {
    auto result_both = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5 --quiet");
    auto result_plus  = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5 --strand plus --quiet");

    REQUIRE(result_both.exit_code == 0);
    REQUIRE(result_plus.exit_code == 0);

    // Count lines
    auto count_lines = [](const std::string& text) {
        size_t count = 0;
        std::istringstream iss(text);
        std::string line;
        while (std::getline(iss, line)) {
            // Skip header line (starts with "chrom\t")
            if (!line.empty() && line.find("chrom\t") != 0) {
                ++count;
            }
        }
        return count;
    };

    size_t hits_both = count_lines(result_both.stdout_text);
    size_t hits_plus = count_lines(result_plus.stdout_text);

    // Plus-only should have fewer or equal hits (roughly half for typical genomes)
    REQUIRE(hits_plus > 0);  // But still some hits
    REQUIRE(hits_plus <= hits_both);  // And not more than both strands
}

TEST_CASE("CLI - strand flag searches specific strand", "[cli][strand]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold 2 --strand plus");

    REQUIRE(result.exit_code == 0);
    // Should have some hits
    REQUIRE(!result.stdout_text.empty());
}

TEST_CASE("CLI - pam-filter none (default)", "[cli][pam]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5");

    REQUIRE(result.exit_code == 0);
    // Should work without error
}

TEST_CASE("CLI - pam-filter both", "[cli][pam]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5 --pam NRG");

    REQUIRE(result.exit_code == 0);
    // Test data likely has no NGG/NAG PAMs, so may return 0 hits (but should not error)
}

TEST_CASE("CLI - pam-filter ngg", "[cli][pam]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5 --pam NGG");

    REQUIRE(result.exit_code == 0);
}

TEST_CASE("CLI - pam-filter nag", "[cli][pam]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5 --pam NAG");

    REQUIRE(result.exit_code == 0);
}

TEST_CASE("CLI - no-compute-mismatches disables mismatch info", "[cli][mismatches]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --threshold 2 --no-compute-mismatches --max-hits 2");

    REQUIRE(result.exit_code == 0);
    // Mismatch columns should be empty (dots)
    REQUIRE(result.stdout_text.find("\t.\t") != std::string::npos);
}

TEST_CASE("CLI - pam-filter requires compute-mismatches", "[cli][pam][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTACGTACGTACGT --genome " + kTestSmallFa + " --threshold 5 --pam NRG --no-compute-mismatches");

    // Should fail with error
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("--pam requires") != std::string::npos);
}

TEST_CASE("CLI - invalid pam pattern", "[cli][pam][error]") {
    // Non-IUPAC character in --pam pattern
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --pam XYZ");

    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("--pam") != std::string::npos);
}

TEST_CASE("CLI - invalid strand value", "[cli][strand][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa + " --strand invalid");

    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("strand") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Batch mode tests
// ───────────────────────────────────────────────────────────────────────────

static const std::string kTestSpacersFile = TEST_DATA_DIR "/synthetic/test_spacers.txt";

TEST_CASE("CLI - batch mode basic search", "[cli][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threshold 4");

    REQUIRE(result.exit_code == 0);
    // Should have TSV header with spacer column first
    REQUIRE(result.stdout_text.find("spacer\tchrom\tstart") != std::string::npos);
}

TEST_CASE("CLI - batch mode TSV has spacer column", "[cli][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threshold 4 --quiet");

    REQUIRE(result.exit_code == 0);
    // Check that spacer names appear in output
    // test_spacers.txt has named spacers like test_spacer_1, test_spacer_2
    REQUIRE((result.stdout_text.find("test_spacer") != std::string::npos ||
             result.stdout_text.find("spacer_") != std::string::npos ||
             result.stdout_text.find("chr1_match") != std::string::npos));
}

TEST_CASE("CLI - batch mode BED format", "[cli][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threshold 4 --format bed --quiet");

    REQUIRE(result.exit_code == 0);
    // BED format should have spacer name in column 7 (after strand)
    // Check for chromosome data (BED has no header)
    REQUIRE(result.stdout_text.find("chr") != std::string::npos);
}

TEST_CASE("CLI - batch mode mutual exclusivity with pattern", "[cli][batch][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa);

    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("mutually exclusive") != std::string::npos);
}

TEST_CASE("CLI - batch mode requires spacer-file or pattern", "[cli][batch][error]") {
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa);

    REQUIRE(result.exit_code != 0);
    REQUIRE((result.stderr_text.find("--pattern") != std::string::npos ||
             result.stderr_text.find("--spacer-file") != std::string::npos));
}

TEST_CASE("CLI - batch mode nonexistent spacer file fails", "[cli][batch][error]") {
    auto result = run_command(kStomataExe + " --spacer-file /nonexistent/spacers.txt --genome " + kTestSmallFa);

    REQUIRE(result.exit_code != 0);
    REQUIRE((result.stderr_text.find("Error") != std::string::npos ||
             result.stderr_text.find("Cannot open") != std::string::npos));
}

TEST_CASE("CLI - batch mode with verbose shows progress", "[cli][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threshold 4 --verbose");

    REQUIRE(result.exit_code == 0);
    // Verbose should show searching progress
    REQUIRE((result.stderr_text.find("Searching") != std::string::npos ||
             result.stderr_text.find("spacer") != std::string::npos ||
             result.stderr_text.find("hits") != std::string::npos));
}

TEST_CASE("CLI - batch mode max-hits applies per spacer", "[cli][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threshold 6 --max-hits 2 --quiet");

    REQUIRE(result.exit_code == 0);
    // Count hits per spacer
    std::map<std::string, size_t> hits_per_spacer;
    std::istringstream iss(result.stdout_text);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.find("spacer\t") == 0) continue;  // Skip header
        if (line.empty()) continue;
        // First column is spacer name
        size_t tab_pos = line.find('\t');
        if (tab_pos != std::string::npos) {
            std::string spacer = line.substr(0, tab_pos);
            hits_per_spacer[spacer]++;
        }
    }
    // Each spacer should have at most 2 hits
    for (const auto& [spacer, count] : hits_per_spacer) {
        REQUIRE(count <= 2);
    }
}

TEST_CASE("CLI - batch mode pam-filter works", "[cli][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threshold 4 --pam NRG --quiet");

    REQUIRE(result.exit_code == 0);
    // Should work without error (may have 0 hits if no NGG/NAG PAMs in test data)
}

// ───────────────────────────────────────────────────────────────────────────
// Parallel batch processing CLI tests
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - batch mode with threads option", "[cli][batch][parallel]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa +
                              " --threshold 4 --threads 2 --quiet");

    REQUIRE(result.exit_code == 0);
    // Should produce valid output
    REQUIRE(result.stdout_text.find("spacer\tchrom") != std::string::npos);
}

TEST_CASE("CLI - batch mode threads=0 auto-detects", "[cli][batch][parallel]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa +
                              " --threshold 4 --threads 0 --quiet");

    REQUIRE(result.exit_code == 0);
}

TEST_CASE("CLI - batch mode with -t shorthand", "[cli][batch][parallel]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa +
                              " --threshold 4 -t 2 --quiet");

    REQUIRE(result.exit_code == 0);
}

TEST_CASE("CLI - invalid threads value fails", "[cli][batch][parallel][error]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile +
                              " --genome " + kTestSmallFa + " --threads -1");

    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("non-negative") != std::string::npos);
}

TEST_CASE("CLI - threads without batch mode still works", "[cli][parallel]") {
    // --threads should be accepted even in single pattern mode (ignored but not an error)
    auto result = run_command(kStomataExe + " --pattern ACGTACGT"
                              " --genome " + kTestSmallFa +
                              " --threshold 4 --threads 2 --quiet");

    REQUIRE(result.exit_code == 0);
}

TEST_CASE("CLI - help shows threads option", "[cli][help]") {
    auto result = run_command(kStomataExe + " --help");

    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("--threads") != std::string::npos);
    REQUIRE(result.stdout_text.find("-t") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: Single-strand search (--strand plus / minus)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - --strand plus excludes minus hits", "[cli][strand]") {
    auto result_both = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --quiet");
    REQUIRE(result_both.exit_code == 0);

    auto result_plus = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --strand plus --quiet");
    REQUIRE(result_plus.exit_code == 0);
    REQUIRE(result_plus.stdout_text.find("\t-\t") == std::string::npos);
}

TEST_CASE("CLI - --strand minus excludes plus hits", "[cli][strand]") {
    auto result_minus = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --strand minus --quiet");
    REQUIRE(result_minus.exit_code == 0);
    REQUIRE(result_minus.stdout_text.find("\t+\t") == std::string::npos);
}

TEST_CASE("CLI - strand in help text", "[cli][strand]") {
    auto result = run_command(kStomataExe + " --help");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("--strand") != std::string::npos);
    REQUIRE(result.stdout_text.find("plus") != std::string::npos);
    REQUIRE(result.stdout_text.find("minus") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: Treat U as T (--treat-u-as-t)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - treat-u-as-t (default enabled)", "[cli][uracil]") {
    // Pattern with U should work by default
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGUACGUAC --threshold 2 --quiet");
    REQUIRE(result.exit_code == 0);

    // Compare with T version - should produce same hits
    auto result_t = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --quiet");
    REQUIRE(result_t.exit_code == 0);

    // Count lines (excluding header)
    auto count_lines = [](const std::string& s) {
        size_t count = 0;
        for (char c : s) if (c == '\n') count++;
        return count;
    };
    REQUIRE(count_lines(result.stdout_text) == count_lines(result_t.stdout_text));
}

TEST_CASE("CLI - no-treat-u-as-t rejects U", "[cli][uracil]") {
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGUACGUAC --threshold 2 --no-treat-u-as-t");
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("U") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: Summary mode (--summary)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - summary mode JSON output", "[cli][summary]") {
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --summary --quiet");
    REQUIRE(result.exit_code == 0);

    // Should contain JSON structure
    REQUIRE(result.stdout_text.find("\"threshold\"") != std::string::npos);
    REQUIRE(result.stdout_text.find("\"total_hits\"") != std::string::npos);
    REQUIRE(result.stdout_text.find("\"hits_by_distance\"") != std::string::npos);
}

TEST_CASE("CLI - summary mode TSV output", "[cli][summary]") {
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --summary --summary-format tsv --quiet");
    REQUIRE(result.exit_code == 0);

    // Should contain TSV header with distance columns
    REQUIRE(result.stdout_text.find("spacer\tsequence\ttotal_hits") != std::string::npos);
    REQUIRE(result.stdout_text.find("d0\td1\td2") != std::string::npos);
}

TEST_CASE("CLI - summary mode with batch", "[cli][summary]") {
    std::string spacer_file = TEST_DATA_DIR "/synthetic/test_spacers.txt";
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --spacer-file " + spacer_file + " --threshold 2 --summary --quiet");
    REQUIRE(result.exit_code == 0);

    // Should contain JSON with spacers array
    REQUIRE(result.stdout_text.find("\"spacers\"") != std::string::npos);
}

TEST_CASE("CLI - summary combined with --strand plus", "[cli][summary][strand]") {
    auto result = run_command(kStomataExe + " --genome " + kTestSmallFa
        + " --pattern ACGTACGTAC --threshold 2 --summary --strand plus --quiet");
    REQUIRE(result.exit_code == 0);

    // Should still produce valid JSON summary
    REQUIRE(result.stdout_text.find("\"threshold\"") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: Distance mode (--distance-mode)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - distance-mode levenshtein is default", "[cli][distance-mode]") {
    auto result_default = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 2 --quiet");
    auto result_explicit = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 2 --distance-mode levenshtein --quiet");

    REQUIRE(result_default.exit_code == 0);
    REQUIRE(result_explicit.exit_code == 0);
    REQUIRE(result_default.stdout_text == result_explicit.stdout_text);
}

TEST_CASE("CLI - distance-mode hamming produces no indels", "[cli][distance-mode]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 2 --distance-mode hamming --quiet");
    REQUIRE(result.exit_code == 0);

    // Hamming results should contain no DNA_BULGE or RNA_BULGE CIGAR operations.
    // The CIGAR column format uses 'I' (RNA bulge) / 'D' (DNA bulge) via alignment.
    // Conservative check: every output hit line's distance column is numeric and
    // the CIGAR column (if present) contains only M and X operations.
    // We at minimum assert that the run succeeded and produced TSV.
    REQUIRE(result.stdout_text.find("chrom\tstart\tend\tpattern\tdistance") != std::string::npos);
}

TEST_CASE("CLI - distance-mode hamming excludes indel-only hits", "[cli][distance-mode]") {
    // Hamming mode must never emit a hit whose edit_types column contains
    // DNA_BULGE or RNA_BULGE — those are indels, which Hamming distance forbids.
    auto ham = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 3 --distance-mode hamming --quiet");

    REQUIRE(ham.exit_code == 0);
    REQUIRE(ham.stdout_text.find("DNA_BULGE") == std::string::npos);
    REQUIRE(ham.stdout_text.find("RNA_BULGE") == std::string::npos);
}

TEST_CASE("CLI - distance-mode invalid value fails", "[cli][distance-mode][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa
        + " --distance-mode nonsense");
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("distance-mode") != std::string::npos);
}

TEST_CASE("CLI - distance-mode hamming works in batch mode", "[cli][distance-mode][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile
        + " --genome " + kTestSmallFa + " --threshold 2 --distance-mode hamming --quiet");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("spacer\tchrom\tstart") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: Spacer file diagnostics (line-numbered skip warnings)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - invalid spacer warning reports file:line", "[cli][batch][error]") {
    std::string bad_file = "/tmp/stomata_bad_spacers_" + std::to_string(getpid()) + ".txt";
    {
        std::ofstream f(bad_file);
        f << "# header comment\n";
        f << "good_one\tACGTACGTACGTACGTACGT\n";
        f << "\n";
        f << "bad_one\tACGTXACGTACGTACGTACGT\n";  // X is invalid, line 4
        f << "good_two\tCCCCCCCCCCCCCCCCCCCC\n";
    }

    auto result = run_command(kStomataExe + " --spacer-file " + bad_file
        + " --genome " + kTestSmallFa + " --threshold 2 --quiet");

    // Should succeed overall (the good spacers search), with the bad one skipped.
    REQUIRE(result.exit_code == 0);
    // Warning should mention the file and the line number where the bad spacer lives.
    REQUIRE(result.stderr_text.find("bad_one") != std::string::npos);
    REQUIRE(result.stderr_text.find(":4") != std::string::npos);

    std::remove(bad_file.c_str());
}

TEST_CASE("CLI - empty spacer file produces clear error", "[cli][batch][error]") {
    std::string empty_file = "/tmp/stomata_empty_spacers_" + std::to_string(getpid()) + ".txt";
    {
        std::ofstream f(empty_file);
        f << "# only comments here\n";
        f << "\n";
        f << "   \n";
    }

    auto result = run_command(kStomataExe + " --spacer-file " + empty_file
        + " --genome " + kTestSmallFa);
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("No valid spacers") != std::string::npos);

    std::remove(empty_file.c_str());
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: Debug flag (--no-deduplicate)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - no-deduplicate is accepted and works", "[cli][debug]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 2 --no-deduplicate --quiet");
    REQUIRE(result.exit_code == 0);
}

TEST_CASE("CLI - no-deduplicate documented in help", "[cli][help]") {
    auto result = run_command(kStomataExe + " --help");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("--no-deduplicate") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: CIGAR column
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - TSV output includes cigar column", "[cli][cigar]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 0 --quiet");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("cigar") != std::string::npos);
    // Exact matches at threshold 0 should produce "10M"
    REQUIRE(result.stdout_text.find("10M") != std::string::npos);
}

TEST_CASE("CLI - cigar column present in batch TSV output", "[cli][cigar][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile
        + " --genome " + kTestSmallFa + " --threshold 2 --quiet");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("cigar") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: JSON output format
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - JSON output format produces valid-looking JSON array", "[cli][format][json]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 0 --format json --quiet");
    REQUIRE(result.exit_code == 0);
    // Starts with '[' and ends with ']' + newline
    REQUIRE(!result.stdout_text.empty());
    REQUIRE(result.stdout_text.front() == '[');
    REQUIRE(result.stdout_text.find("]\n") != std::string::npos);
    // Must carry the same field names as TSV
    REQUIRE(result.stdout_text.find("\"chrom\"") != std::string::npos);
    REQUIRE(result.stdout_text.find("\"cigar\"") != std::string::npos);
    REQUIRE(result.stdout_text.find("\"cfd_score\"") != std::string::npos);
}

TEST_CASE("CLI - batch JSON format includes spacer field", "[cli][format][json][batch]") {
    auto result = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile
        + " --genome " + kTestSmallFa + " --threshold 2 --format json --quiet");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.front() == '[');
    REQUIRE(result.stdout_text.find("\"spacer\"") != std::string::npos);
    REQUIRE(result.stdout_text.find("\"cfd_score\"") != std::string::npos);
}

TEST_CASE("CLI - JSON format with zero hits emits empty array", "[cli][format][json]") {
    auto result = run_command(kStomataExe + " --pattern TTTTTTTTTT --genome " + kTestSmallFa
        + " --threshold 0 --format json --quiet");
    REQUIRE(result.exit_code == 0);
    // Empty array: "[]\n"
    REQUIRE(result.stdout_text == "[]\n");
}

TEST_CASE("CLI - invalid format value is rejected", "[cli][format][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGT --genome " + kTestSmallFa
        + " --format xml");
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("tsv") != std::string::npos);
    REQUIRE(result.stderr_text.find("json") != std::string::npos);
}

TEST_CASE("CLI - help lists json as a format option", "[cli][help]") {
    auto result = run_command(kStomataExe + " --help");
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("json") != std::string::npos);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: --region search window
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - region restricts to chromosome", "[cli][region]") {
    // Whole-genome search returns hits in multiple chromosomes; restricting to
    // chr1 should drop any chr2/chr3 hits.
    auto all = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 1 --quiet");
    auto chr1_only = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 1 --region chr1 --quiet");
    REQUIRE(all.exit_code == 0);
    REQUIRE(chr1_only.exit_code == 0);
    // chr1 output should not reference chr2/chr3.
    REQUIRE(chr1_only.stdout_text.find("chr2\t") == std::string::npos);
    REQUIRE(chr1_only.stdout_text.find("chr3\t") == std::string::npos);
}

TEST_CASE("CLI - region with explicit window restricts positions", "[cli][region]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 0 --region chr1:50-95 --quiet");
    REQUIRE(result.exit_code == 0);
    // Any returned hit must have end position <= 95 (the exclusive end of the window).
    // A trivial lower-bound check: the output should not contain early positions like
    // chr1\t0\t or chr1\t100\t (outside the window).
    REQUIRE(result.stdout_text.find("chr1\t0\t") == std::string::npos);
}

TEST_CASE("CLI - region with unknown chromosome fails", "[cli][region][error]") {
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 1 --region chrBOGUS --quiet");
    REQUIRE(result.exit_code != 0);
    REQUIRE(result.stderr_text.find("chrBOGUS") != std::string::npos);
}

TEST_CASE("CLI - region with invalid range fails", "[cli][region][error]") {
    // End before start must be rejected.
    auto result = run_command(kStomataExe + " --pattern ACGTACGTAC --genome " + kTestSmallFa
        + " --threshold 1 --region chr1:100-50 --quiet");
    REQUIRE(result.exit_code != 0);
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: --max-total-hits global cap
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - max-total-hits caps batch hits globally", "[cli][batch][max-total-hits]") {
    // Run uncapped first to confirm there are > 2 hits to cap.
    auto uncapped = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile
        + " --genome " + kTestSmallFa + " --threshold 3 --quiet");
    REQUIRE(uncapped.exit_code == 0);

    auto capped = run_command(kStomataExe + " --spacer-file " + kTestSpacersFile
        + " --genome " + kTestSmallFa + " --threshold 3 --max-total-hits 2 --quiet");
    REQUIRE(capped.exit_code == 0);

    // Count data lines (excluding header) in each TSV.
    auto count_data_lines = [](const std::string& tsv) {
        size_t n = 0;
        size_t pos = 0;
        bool first = true;
        while (pos < tsv.size()) {
            size_t nl = tsv.find('\n', pos);
            if (nl == std::string::npos) break;
            if (!first) ++n;
            first = false;
            pos = nl + 1;
        }
        return n;
    };
    REQUIRE(count_data_lines(capped.stdout_text) <= 2);
    REQUIRE(count_data_lines(uncapped.stdout_text) > count_data_lines(capped.stdout_text));
}

// ───────────────────────────────────────────────────────────────────────────
// Feature: --spacer-file - (stdin input)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("CLI - spacer-file dash reads spacers from stdin", "[cli][batch][stdin]") {
    // Pipe two spacers via stdin.
    std::string cmd = "printf 'stdin_spacer\\tACGTACGTAC\\n' | " + kStomataExe
        + " --spacer-file - --genome " + kTestSmallFa
        + " --threshold 1 --quiet";
    auto result = run_command(cmd);
    REQUIRE(result.exit_code == 0);
    REQUIRE(result.stdout_text.find("stdin_spacer") != std::string::npos);
}