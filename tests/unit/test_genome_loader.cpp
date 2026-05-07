#include <catch2/catch_test_macros.hpp>
#include <genome_loader.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef TEST_DATA_DIR
#error "TEST_DATA_DIR not defined. Set via CMake target_compile_definitions."
#endif

namespace fs = std::filesystem;

// Path to the shared two-chromosome FASTA used by several tests below.
static constexpr const char* SMALL_FASTA    = TEST_DATA_DIR "/synthetic/test_small.fa";
static constexpr const char* SMALL_FASTA_GZ = TEST_DATA_DIR "/synthetic/test_small.fa.gz";

// Query a single bit from a packed 64-bit bit-vector.
// Bit layout: position p lives at bit (p % 64) of word (p / 64).
static bool get_bit(const std::vector<uint64_t>& bv, size_t pos) {
    return (bv[pos / 64] >> (pos % 64)) & 1;
}

// Write a temporary FASTA file and return its path.
// Used by parsing tests that need precise, small inputs without polluting
// the checked-in test-data directory.
static std::string write_temp_fasta(const std::string& content) {
    auto path = fs::temp_directory_path() / "stomata_test.fa";
    std::ofstream out(path);
    out << content;
    return path.string();
}

// ───────────────────────────────────────────────────────────────────────────
// FASTA Parsing
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("FASTA parsing - single sequence", "[genome_loader][fasta]") {
    std::string path = write_temp_fasta(
        ">seq1 A single short sequence\n"
        "ACGTACGT\n");

    auto entries = parse_fasta(path);

    REQUIRE(entries.size() == 1);
    REQUIRE(entries[0].name     == "seq1");
    REQUIRE(entries[0].sequence == "ACGTACGT");
}

TEST_CASE("FASTA parsing - multiple sequences (chromosomes)", "[genome_loader][fasta]") {
    // test_small.fa contains chr1 (100 bp) and chr2 (80 bp)
    auto entries = parse_fasta(SMALL_FASTA);

    REQUIRE(entries.size() == 2);

    REQUIRE(entries[0].name == "chr1");
    REQUIRE(entries[0].sequence.size() == 100);

    REQUIRE(entries[1].name == "chr2");
    REQUIRE(entries[1].sequence.size() == 80);
}

TEST_CASE("FASTA parsing - comment lines are skipped", "[genome_loader][fasta]") {
    std::string path = write_temp_fasta(
        "; File-level comment before any sequence\n"
        ">seq1 Test\n"
        "ACGT\n"
        "; Comment between sequence lines\n"
        "TGCA\n");

    auto entries = parse_fasta(path);

    REQUIRE(entries.size() == 1);
    // Comments must not appear in the sequence; the two lines concatenate
    REQUIRE(entries[0].sequence == "ACGTTGCA");
}

TEST_CASE("FASTA parsing - wrapped lines are concatenated", "[genome_loader][fasta]") {
    std::string path = write_temp_fasta(
        ">seq1 Wrapped\n"
        "AAAA\n"
        "CCCC\n"
        "GGGG\n"
        "TTTT\n");

    auto entries = parse_fasta(path);

    REQUIRE(entries.size() == 1);
    REQUIRE(entries[0].sequence == "AAAACCCCGGGGTTTT");
}

TEST_CASE("FASTA parsing - nonexistent file throws", "[genome_loader][fasta]") {
    REQUIRE_THROWS_AS(parse_fasta("/nonexistent/path/missing.fa"),
                      std::runtime_error);
}

// ───────────────────────────────────────────────────────────────────────────
// Nucleotide Encoding
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Nucleotide encoding - uppercase ACGT", "[genome_loader][encoding]") {
    // Bit layout per PROJECT_CONTEXT decision log:
    //   bit 0 = A  (0b0001)
    //   bit 1 = C  (0b0010)
    //   bit 2 = G  (0b0100)
    //   bit 3 = T  (0b1000)
    REQUIRE(encode_nucleotide('A') == 0b0001);
    REQUIRE(encode_nucleotide('C') == 0b0010);
    REQUIRE(encode_nucleotide('G') == 0b0100);
    REQUIRE(encode_nucleotide('T') == 0b1000);
}

TEST_CASE("Nucleotide encoding - lowercase acgt", "[genome_loader][encoding]") {
    REQUIRE(encode_nucleotide('a') == 0b0001);
    REQUIRE(encode_nucleotide('c') == 0b0010);
    REQUIRE(encode_nucleotide('g') == 0b0100);
    REQUIRE(encode_nucleotide('t') == 0b1000);
}

TEST_CASE("Nucleotide encoding - N in pattern is wildcard", "[genome_loader][encoding]") {
    REQUIRE(encode_nucleotide('N') == 0b1111);
    REQUIRE(encode_nucleotide('n') == 0b1111);
}

TEST_CASE("Genome nucleotide encoding - N is masked (no bits set)", "[genome_loader][encoding]") {
    REQUIRE(encode_genome_nucleotide('N') == 0b0000);
    REQUIRE(encode_genome_nucleotide('n') == 0b0000);
    REQUIRE(encode_genome_nucleotide('A') == 0b0001);
    REQUIRE(encode_genome_nucleotide('C') == 0b0010);
    REQUIRE(encode_genome_nucleotide('G') == 0b0100);
    REQUIRE(encode_genome_nucleotide('T') == 0b1000);
}

TEST_CASE("Nucleotide encoding - invalid characters throw", "[genome_loader][encoding]") {
    REQUIRE_THROWS_AS(encode_nucleotide('X'), std::invalid_argument);
    REQUIRE_THROWS_AS(encode_nucleotide('Y'), std::invalid_argument);
    REQUIRE_THROWS_AS(encode_nucleotide('Z'), std::invalid_argument);
    REQUIRE_THROWS_AS(encode_nucleotide('1'), std::invalid_argument);
    REQUIRE_THROWS_AS(encode_nucleotide(' '), std::invalid_argument);
}

TEST_CASE("Bit-vector encoding - ACGT positions are mutually exclusive",
          "[genome_loader][encoding]") {
    std::vector<FastaEntry> entries = {{"test", "ACGT"}};
    Genome genome = encode_genome(entries);

    // Each position must have exactly one bit set across the four vectors.
    SECTION("Position 0 is A") {
        REQUIRE( get_bit(genome.bv_A, 0));
        REQUIRE(!get_bit(genome.bv_C, 0));
        REQUIRE(!get_bit(genome.bv_G, 0));
        REQUIRE(!get_bit(genome.bv_T, 0));
    }

    SECTION("Position 1 is C") {
        REQUIRE(!get_bit(genome.bv_A, 1));
        REQUIRE( get_bit(genome.bv_C, 1));
        REQUIRE(!get_bit(genome.bv_G, 1));
        REQUIRE(!get_bit(genome.bv_T, 1));
    }

    SECTION("Position 2 is G") {
        REQUIRE(!get_bit(genome.bv_A, 2));
        REQUIRE(!get_bit(genome.bv_C, 2));
        REQUIRE( get_bit(genome.bv_G, 2));
        REQUIRE(!get_bit(genome.bv_T, 2));
    }

    SECTION("Position 3 is T") {
        REQUIRE(!get_bit(genome.bv_A, 3));
        REQUIRE(!get_bit(genome.bv_C, 3));
        REQUIRE(!get_bit(genome.bv_G, 3));
        REQUIRE( get_bit(genome.bv_T, 3));
    }
}

TEST_CASE("Bit-vector encoding - N in genome is masked (no bits set)",
          "[genome_loader][encoding]") {
    std::vector<FastaEntry> entries = {{"test", "ANA"}};
    Genome genome = encode_genome(entries);

    // Position 1 is N: no bits set (masked base, matches nothing)
    REQUIRE(!get_bit(genome.bv_A, 1));
    REQUIRE(!get_bit(genome.bv_C, 1));
    REQUIRE(!get_bit(genome.bv_G, 1));
    REQUIRE(!get_bit(genome.bv_T, 1));

    // Flanking A bases are still encoded correctly (only bv_A set)
    REQUIRE( get_bit(genome.bv_A, 0));
    REQUIRE(!get_bit(genome.bv_C, 0));
    REQUIRE( get_bit(genome.bv_A, 2));
    REQUIRE(!get_bit(genome.bv_C, 2));
}

TEST_CASE("Bit-vector encoding - lowercase input is accepted",
          "[genome_loader][encoding]") {
    std::vector<FastaEntry> entries = {{"test", "acgt"}};
    Genome genome = encode_genome(entries);

    REQUIRE( get_bit(genome.bv_A, 0));
    REQUIRE( get_bit(genome.bv_C, 1));
    REQUIRE( get_bit(genome.bv_G, 2));
    REQUIRE( get_bit(genome.bv_T, 3));
}

TEST_CASE("Bit-vector encoding - invalid character in sequence throws",
          "[genome_loader][encoding]") {
    std::vector<FastaEntry> entries = {{"bad", "ACXGT"}};
    REQUIRE_THROWS_AS(encode_genome(entries), std::invalid_argument);
}

// ───────────────────────────────────────────────────────────────────────────
// Chromosome Metadata
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("Chromosome metadata - names are extracted from FASTA headers",
          "[genome_loader][metadata]") {
    std::vector<FastaEntry> entries = {
        {"chr1", std::string(10, 'A')},
        {"chrX", std::string(20, 'C')}
    };
    Genome genome = encode_genome(entries);

    REQUIRE(genome.chromosomes.size() == 2);
    REQUIRE(genome.chromosomes[0].name == "chr1");
    REQUIRE(genome.chromosomes[1].name == "chrX");
}

TEST_CASE("Chromosome metadata - start/length positions are contiguous",
          "[genome_loader][metadata]") {
    std::vector<FastaEntry> entries = {
        {"chr1", std::string(10, 'A')},   // 10 bases
        {"chr2", std::string(25, 'C')},   // 25 bases
        {"chr3", std::string( 5, 'G')}    //  5 bases
    };
    Genome genome = encode_genome(entries);

    REQUIRE(genome.chromosomes[0].start  == 0);
    REQUIRE(genome.chromosomes[0].length == 10);

    REQUIRE(genome.chromosomes[1].start  == 10);
    REQUIRE(genome.chromosomes[1].length == 25);

    REQUIRE(genome.chromosomes[2].start  == 35);
    REQUIRE(genome.chromosomes[2].length == 5);

    REQUIRE(genome.total_bases == 40);
}

TEST_CASE("Chromosome metadata - chromosomes of different lengths",
          "[genome_loader][metadata]") {
    // Lengths deliberately span multiple 64-bit word boundaries
    std::vector<FastaEntry> entries = {
        {"tiny",   std::string(  4, 'A')},   //   4 bases  – well under one word
        {"medium", std::string( 64, 'C')},   //  64 bases  – exactly one word
        {"large",  std::string(200, 'G')}    // 200 bases  – spans multiple words
    };
    Genome genome = encode_genome(entries);

    REQUIRE(genome.chromosomes[0].length ==   4);
    REQUIRE(genome.chromosomes[1].length ==  64);
    REQUIRE(genome.chromosomes[2].length == 200);
    REQUIRE(genome.total_bases           == 268);

    // Bit-vectors must be large enough for the full genome: ceil(268/64) = 5 words
    size_t expected_words = (268 + 63) / 64;
    REQUIRE(genome.bv_A.size() == expected_words);
    REQUIRE(genome.bv_C.size() == expected_words);
    REQUIRE(genome.bv_G.size() == expected_words);
    REQUIRE(genome.bv_T.size() == expected_words);
}

// ───────────────────────────────────────────────────────────────────────────
// Gzipped FASTA
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("gzip FASTA parsing - matches plain FASTA output",
          "[genome_loader][fasta][gzip]") {
    auto plain = parse_fasta(SMALL_FASTA);
    auto gz    = parse_fasta(SMALL_FASTA_GZ);

    REQUIRE(gz.size() == plain.size());
    for (size_t i = 0; i < plain.size(); ++i) {
        REQUIRE(gz[i].name     == plain[i].name);
        REQUIRE(gz[i].sequence == plain[i].sequence);
    }
}

TEST_CASE("gzip FASTA parsing - nonexistent .gz file throws",
          "[genome_loader][fasta][gzip]") {
    REQUIRE_THROWS_AS(parse_fasta("/nonexistent/path/missing.fa.gz"),
                      std::runtime_error);
}

TEST_CASE("gzip FASTA parsing - corrupt .gz file throws",
          "[genome_loader][fasta][gzip]") {
    auto path = fs::temp_directory_path() / "stomata_test_corrupt.fa.gz";
    {
        std::ofstream out(path, std::ios::binary);
        out.write("this is not valid gzip data", 27);
    }
    REQUIRE_THROWS_AS(parse_fasta(path.string()), std::runtime_error);
}

TEST_CASE("gzip load_fasta - end-to-end produces same Genome as plain",
          "[genome_loader][integration][gzip]") {
    Genome plain = load_fasta(SMALL_FASTA);
    Genome gz    = load_fasta(SMALL_FASTA_GZ);

    REQUIRE(gz.total_bases        == plain.total_bases);
    REQUIRE(gz.chromosomes.size() == plain.chromosomes.size());
    for (size_t i = 0; i < plain.chromosomes.size(); ++i) {
        REQUIRE(gz.chromosomes[i].name   == plain.chromosomes[i].name);
        REQUIRE(gz.chromosomes[i].start  == plain.chromosomes[i].start);
        REQUIRE(gz.chromosomes[i].length == plain.chromosomes[i].length);
    }
    REQUIRE(gz.bv_A == plain.bv_A);
    REQUIRE(gz.bv_C == plain.bv_C);
    REQUIRE(gz.bv_G == plain.bv_G);
    REQUIRE(gz.bv_T == plain.bv_T);
}

// ───────────────────────────────────────────────────────────────────────────
// Integration: load_fasta end-to-end on the synthetic test file
// ───────────────────────────────────────────────────────────────────────────
// test_small.fa layout (see tests/data/synthetic/test_small.fa):
//   chr1 – 100 bp: (ACGT)×14 + NNNN + (ACGT)×10
//   chr2 –  80 bp: (TGCA)×20
//   Total: 180 bp → ceil(180/64) = 3 uint64_t words per bit-vector

TEST_CASE("load_fasta - parses and encodes test_small.fa end-to-end",
          "[genome_loader][integration]") {
    Genome genome = load_fasta(SMALL_FASTA);

    // ── overall dimensions ──
    REQUIRE(genome.total_bases      == 180);
    REQUIRE(genome.chromosomes.size() == 2);

    REQUIRE(genome.chromosomes[0].name   == "chr1");
    REQUIRE(genome.chromosomes[0].start  == 0);
    REQUIRE(genome.chromosomes[0].length == 100);

    REQUIRE(genome.chromosomes[1].name   == "chr2");
    REQUIRE(genome.chromosomes[1].start  == 100);
    REQUIRE(genome.chromosomes[1].length == 80);

    // ── bit-vector word count ──
    REQUIRE(genome.bv_A.size() == 3);
    REQUIRE(genome.bv_C.size() == 3);
    REQUIRE(genome.bv_G.size() == 3);
    REQUIRE(genome.bv_T.size() == 3);

    SECTION("chr1 begins with repeating ACGT") {
        REQUIRE( get_bit(genome.bv_A, 0));  // A
        REQUIRE( get_bit(genome.bv_C, 1));  // C
        REQUIRE( get_bit(genome.bv_G, 2));  // G
        REQUIRE( get_bit(genome.bv_T, 3));  // T
    }

    SECTION("chr1 N bases at positions 56-59 are masked (no bits set)") {
        for (size_t i = 56; i < 60; ++i) {
            REQUIRE(!get_bit(genome.bv_A, i));
            REQUIRE(!get_bit(genome.bv_C, i));
            REQUIRE(!get_bit(genome.bv_G, i));
            REQUIRE(!get_bit(genome.bv_T, i));
        }
    }

    SECTION("sequence resumes with ACGT immediately after the N run") {
        // Position 60 is the first base of the second line of chr1
        REQUIRE( get_bit(genome.bv_A, 60));  // A
        REQUIRE( get_bit(genome.bv_C, 61));  // C
        REQUIRE( get_bit(genome.bv_G, 62));  // G
        REQUIRE( get_bit(genome.bv_T, 63));  // T
    }

    SECTION("chr2 starts with repeating TGCA at genome offset 100") {
        REQUIRE( get_bit(genome.bv_T, 100));  // T
        REQUIRE( get_bit(genome.bv_G, 101));  // G
        REQUIRE( get_bit(genome.bv_C, 102));  // C
        REQUIRE( get_bit(genome.bv_A, 103));  // A
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Stomata Index (.st format)
// ───────────────────────────────────────────────────────────────────────────

TEST_CASE("is_stomata_index - detects .st extension", "[genome_loader][index]") {
    REQUIRE(is_stomata_index("genome.st"));
    REQUIRE(is_stomata_index("/path/to/hg38.fa.st"));
    REQUIRE(is_stomata_index("a.st"));
    REQUIRE(is_stomata_index(".st"));  // Hidden file with .st extension
    REQUIRE(!is_stomata_index("genome.fa"));
    REQUIRE(!is_stomata_index("genome.fa.gz"));
    REQUIRE(!is_stomata_index("stomata"));  // Too short, no dot prefix
}

TEST_CASE("Index write/read round-trip - preserves genome data",
          "[genome_loader][index][integration]") {
    // Load original genome
    Genome original = load_fasta(SMALL_FASTA);

    // Write to temporary index file
    auto index_path = fs::temp_directory_path() / "stomata_test_roundtrip.st";
    write_genome_index(original, index_path.string());

    // Read back via memory mapping
    MappedGenome mapped = load_genome_index(index_path.string());

    // Verify dimensions match
    REQUIRE(mapped.total_bases == original.total_bases);
    REQUIRE(mapped.num_words == (original.total_bases + 63) / 64);
    REQUIRE(mapped.chromosomes.size() == original.chromosomes.size());

    // Verify chromosome metadata
    for (size_t i = 0; i < original.chromosomes.size(); ++i) {
        REQUIRE(mapped.chromosomes[i].name == original.chromosomes[i].name);
        REQUIRE(mapped.chromosomes[i].start == original.chromosomes[i].start);
        REQUIRE(mapped.chromosomes[i].length == original.chromosomes[i].length);
    }

    // Verify bit-vectors are identical (zero-copy from mmap)
    for (size_t w = 0; w < mapped.num_words; ++w) {
        REQUIRE(mapped.bv_A[w] == original.bv_A[w]);
        REQUIRE(mapped.bv_C[w] == original.bv_C[w]);
        REQUIRE(mapped.bv_G[w] == original.bv_G[w]);
        REQUIRE(mapped.bv_T[w] == original.bv_T[w]);
    }

    // Clean up
    fs::remove(index_path);
}

TEST_CASE("GenomeView - works with both Genome and MappedGenome",
          "[genome_loader][index]") {
    // Load original genome
    Genome genome = load_fasta(SMALL_FASTA);

    // Write to index and read back
    auto index_path = fs::temp_directory_path() / "stomata_test_view.st";
    write_genome_index(genome, index_path.string());
    MappedGenome mapped = load_genome_index(index_path.string());

    // Create views from both
    GenomeView view_genome = make_view(genome);
    GenomeView view_mapped = make_view(mapped);

    // Both views should provide identical data
    REQUIRE(view_genome.total_bases == view_mapped.total_bases);
    REQUIRE(view_genome.num_words == view_mapped.num_words);
    REQUIRE(view_genome.chromosomes->size() == view_mapped.chromosomes->size());

    // Verify bit-vector pointers return same data
    for (size_t w = 0; w < view_genome.num_words; ++w) {
        REQUIRE(view_genome.bv_A[w] == view_mapped.bv_A[w]);
        REQUIRE(view_genome.bv_C[w] == view_mapped.bv_C[w]);
        REQUIRE(view_genome.bv_G[w] == view_mapped.bv_G[w]);
        REQUIRE(view_genome.bv_T[w] == view_mapped.bv_T[w]);
    }

    // Clean up
    fs::remove(index_path);
}

TEST_CASE("extract_genome_slice - works with GenomeView from MappedGenome",
          "[genome_loader][index]") {
    // Load and create index
    Genome genome = load_fasta(SMALL_FASTA);
    auto index_path = fs::temp_directory_path() / "stomata_test_slice.st";
    write_genome_index(genome, index_path.string());
    MappedGenome mapped = load_genome_index(index_path.string());

    // Extract slices from both
    GenomeView view_genome = make_view(genome);
    GenomeView view_mapped = make_view(mapped);

    // Compare slices at various positions
    REQUIRE(extract_genome_slice(view_genome, 0, 20) ==
            extract_genome_slice(view_mapped, 0, 20));
    REQUIRE(extract_genome_slice(view_genome, 50, 30) ==
            extract_genome_slice(view_mapped, 50, 30));
    REQUIRE(extract_genome_slice(view_genome, 100, 40) ==
            extract_genome_slice(view_mapped, 100, 40));

    // Clean up
    fs::remove(index_path);
}

TEST_CASE("Index loading - invalid magic throws", "[genome_loader][index]") {
    auto path = fs::temp_directory_path() / "stomata_test_bad_magic.st";
    {
        std::ofstream out(path, std::ios::binary);
        // Write invalid magic bytes
        out.write("NOTSTOMA", 8);
        // Pad to minimum header size
        char padding[56] = {0};
        out.write(padding, 56);
    }

    REQUIRE_THROWS_AS(load_genome_index(path.string()), std::runtime_error);

    fs::remove(path);
}

TEST_CASE("Index loading - truncated file throws", "[genome_loader][index]") {
    auto path = fs::temp_directory_path() / "stomata_test_truncated.st";
    {
        std::ofstream out(path, std::ios::binary);
        // Write only partial header (less than 64 bytes)
        out.write("STOMATA\x00\x01\x00", 8);
    }

    REQUIRE_THROWS_AS(load_genome_index(path.string()), std::runtime_error);

    fs::remove(path);
}

TEST_CASE("Index loading - nonexistent file throws", "[genome_loader][index]") {
    REQUIRE_THROWS_AS(load_genome_index("/nonexistent/path/missing.st"),
                      std::runtime_error);
}

TEST_CASE("MappedGenome move semantics", "[genome_loader][index]") {
    // Load and create index
    Genome genome = load_fasta(SMALL_FASTA);
    auto index_path = fs::temp_directory_path() / "stomata_test_move.st";
    write_genome_index(genome, index_path.string());

    // Load initial mapped genome
    MappedGenome mapped1 = load_genome_index(index_path.string());
    REQUIRE(mapped1.total_bases == genome.total_bases);
    REQUIRE(mapped1.bv_A != nullptr);

    // Move construct
    MappedGenome mapped2(std::move(mapped1));
    REQUIRE(mapped2.total_bases == genome.total_bases);
    REQUIRE(mapped2.bv_A != nullptr);
    REQUIRE(mapped1.bv_A == nullptr);  // Moved-from object is cleared

    // Move assign
    MappedGenome mapped3;
    mapped3 = std::move(mapped2);
    REQUIRE(mapped3.total_bases == genome.total_bases);
    REQUIRE(mapped3.bv_A != nullptr);
    REQUIRE(mapped2.bv_A == nullptr);  // Moved-from object is cleared

    // Clean up
    fs::remove(index_path);
}
