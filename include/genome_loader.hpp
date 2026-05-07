#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct FastaEntry {
    std::string name;
    std::string sequence;
};

struct Chromosome {
    std::string name;
    size_t start;
    size_t length;
};

struct Genome {
    std::vector<uint64_t> bv_A;
    std::vector<uint64_t> bv_C;
    std::vector<uint64_t> bv_G;
    std::vector<uint64_t> bv_T;
    std::vector<Chromosome> chromosomes;
    size_t total_bases;
};

// Memory-mapped genome index (.st format)

// Magic bytes: "STOMATA1" (8 bytes)
constexpr char STOMATA_INDEX_MAGIC[8] = {'S', 'T', 'O', 'M', 'A', 'T', 'A', '1'};
constexpr uint64_t STOMATA_INDEX_VERSION = 1;

// Memory-mapped genome (points into mmap'd region, zero-copy access)
struct MappedGenome {
    const uint64_t* bv_A;      // Pointer into mmap'd file
    const uint64_t* bv_C;
    const uint64_t* bv_G;
    const uint64_t* bv_T;
    std::vector<Chromosome> chromosomes;  // Parsed from file header
    size_t total_bases;
    size_t num_words;          // = (total_bases + 63) / 64

    // RAII handle for mmap
    void* mapped_data;
    size_t mapped_size;
    int fd;

    MappedGenome();
    ~MappedGenome();

    // Move-only (no copy)
    MappedGenome(MappedGenome&& other) noexcept;
    MappedGenome& operator=(MappedGenome&& other) noexcept;
    MappedGenome(const MappedGenome&) = delete;
    MappedGenome& operator=(const MappedGenome&) = delete;
};

// Lightweight view that works with either Genome or MappedGenome.
// Does not own any data; the underlying Genome/MappedGenome must outlive the view.
struct GenomeView {
    const uint64_t* bv_A;
    const uint64_t* bv_C;
    const uint64_t* bv_G;
    const uint64_t* bv_T;
    const std::vector<Chromosome>* chromosomes;
    size_t total_bases;
    size_t num_words;  // = (total_bases + 63) / 64
};

// Create a GenomeView from either a Genome or MappedGenome.
GenomeView make_view(const Genome& genome);
GenomeView make_view(const MappedGenome& genome);

// Check if a filepath has the .st extension.
bool is_stomata_index(const std::string& filepath);

// Write a Genome to an .st index file for fast memory-mapped loading.
// The file format includes a header, chromosome table, and contiguous bit-vectors.
// Throws std::runtime_error if the file cannot be written.
void write_genome_index(const Genome& genome, const std::string& filepath);

// Load an .st index file via memory mapping.
// Returns a MappedGenome with zero-copy access to bit-vectors.
// Throws std::runtime_error if the file cannot be opened or is invalid.
MappedGenome load_genome_index(const std::string& filepath);

// Nucleotide encoding

// Returns a 4-bit mask indicating which nucleotides this character can represent.
// Bit layout: bit0=A, bit1=C, bit2=G, bit3=T.
// For pattern/PAM encoding: Returns 0b1111 for N/n (wildcard: matches all four nucleotides).
// For reference genome encoding, use encode_genome_nucleotide instead.
// Throws std::invalid_argument for anything else.
uint8_t encode_nucleotide(char c);

// Returns a 4-bit mask for encoding reference genome nucleotides.
// N/n is treated as a masked base (returns 0b0000, matches nothing).
// Valid nucleotides ACGT return their single-bit representation.
// Throws std::invalid_argument for anything else.
uint8_t encode_genome_nucleotide(char c);

// Parse a FASTA file into one FastaEntry per sequence.
// Name is the first whitespace-delimited token after '>'.
// Lines starting with ';' are treated as comments and skipped.
// Throws std::runtime_error if the file cannot be opened.
std::vector<FastaEntry> parse_fasta(const std::string& filepath);

// Encode parsed FASTA entries into a packed bit-vector Genome.
// Each nucleotide position sets exactly one bit in the corresponding vector;
// N sets no bits (masked base: does not match any nucleotide in pattern).
// Chromosomes are laid out contiguously.
// Throws std::invalid_argument if any sequence contains an invalid character.
Genome encode_genome(const std::vector<FastaEntry>& entries);

// Convenience: parse and encode in one step.
Genome load_fasta(const std::string& filepath);

// Sequence utilities

// Generate the reverse complement of a nucleotide sequence.
// Valid characters: ACGTNacgtn. Case is preserved (A->T, a->t).
// Throws std::invalid_argument for invalid characters.
std::string reverse_complement(const std::string& seq);

// Extract a contiguous slice of the genome as a nucleotide string.
// Returns uppercase ACGTN characters.
// start + length must be <= genome.total_bases.
// Throws std::out_of_range if the slice extends beyond the genome.
std::string extract_genome_slice(GenomeView view, size_t start, size_t length);

// Overload for Genome (convenience wrapper).
inline std::string extract_genome_slice(const Genome& genome, size_t start, size_t length) {
    return extract_genome_slice(make_view(genome), start, length);
}
