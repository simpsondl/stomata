#include <genome_loader.hpp>

#include <array>
#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

// For memory mapping
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// Nucleotide encoding

// For pattern/PAM sequences: N is a wildcard (matches ACGT)
uint8_t encode_nucleotide(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0b0001;
        case 'C': return 0b0010;
        case 'G': return 0b0100;
        case 'T': return 0b1000;
        case 'N': return 0b1111;  // wildcard
        default:
            throw std::invalid_argument(
                std::string("Invalid nucleotide character: '") + c + "'");
    }
}

// For reference genomes: N is a masked base (matches nothing)
uint8_t encode_genome_nucleotide(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0b0001;
        case 'C': return 0b0010;
        case 'G': return 0b0100;
        case 'T': return 0b1000;
        case 'N': return 0b0000;  // masked: no bits set, matches nothing
        default:
            throw std::invalid_argument(
                std::string("Invalid nucleotide character: '") + c + "'");
    }
}

// gzip decompression

static std::string decompress_gzip(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + filepath);
    }

    std::vector<char> compressed((std::istreambuf_iterator<char>(file)),
                                  std::istreambuf_iterator<char>());

    z_stream strm = {};
    // 15 + 16: base window bits + gzip format flag
    if (inflateInit2(&strm, 15 + 16) != Z_OK) {
        throw std::runtime_error("Failed to initialize zlib decompression");
    }

    strm.next_in  = reinterpret_cast<Bytef*>(compressed.data());
    strm.avail_in = static_cast<uInt>(compressed.size());

    std::string decompressed;
    std::array<char, 4096> buffer;
    int ret;
    do {
        strm.next_out  = reinterpret_cast<Bytef*>(buffer.data());
        strm.avail_out = static_cast<uInt>(buffer.size());
        ret = inflate(&strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
            inflateEnd(&strm);
            throw std::runtime_error("gzip decompression failed for: " + filepath);
        }
        decompressed.append(buffer.data(), buffer.size() - strm.avail_out);
    } while (ret != Z_STREAM_END);

    inflateEnd(&strm);
    return decompressed;
}

// FASTA parsing

static std::vector<FastaEntry> parse_fasta_stream(std::istream& in) {
    std::vector<FastaEntry> entries;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == ';') {
            continue;   // skip blank lines and comments
        }

        if (line[0] == '>') {
            // Header line: name is the first whitespace-delimited token after '>'
            entries.emplace_back();
            std::string header = line.substr(1);
            size_t space = header.find_first_of(" \t");
            entries.back().name = header.substr(0, space);
        } else if (!entries.empty()) {
            // Sequence line: append to the current entry
            entries.back().sequence += line;
        }
    }

    return entries;
}

std::vector<FastaEntry> parse_fasta(const std::string& filepath) {
    // Gzipped path: decompress into memory, then parse the string
    if (filepath.size() >= 3 &&
        filepath.compare(filepath.size() - 3, 3, ".gz") == 0) {
        std::string content = decompress_gzip(filepath);
        std::istringstream stream(std::move(content));
        return parse_fasta_stream(stream);
    }

    // Plain path
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + filepath);
    }
    return parse_fasta_stream(file);
}

// Bit-vector encoding

Genome encode_genome(const std::vector<FastaEntry>& entries) {
    size_t total = 0;
    for (const auto& entry : entries) {
        total += entry.sequence.size();
    }

    size_t num_words = (total + 63) / 64;

    Genome genome;
    genome.total_bases = total;
    genome.bv_A.assign(num_words, 0);
    genome.bv_C.assign(num_words, 0);
    genome.bv_G.assign(num_words, 0);
    genome.bv_T.assign(num_words, 0);

    size_t pos = 0;
    for (const auto& entry : entries) {
        Chromosome chrom;
        chrom.name   = entry.name;
        chrom.start  = pos;
        chrom.length = entry.sequence.size();

        for (char c : entry.sequence) {
            uint8_t  mask = encode_genome_nucleotide(c);   // N -> 0b0000 (masked)
            size_t   word = pos / 64;
            uint64_t bit  = 1ULL << (pos % 64);

            if (mask & 0b0001) genome.bv_A[word] |= bit;
            if (mask & 0b0010) genome.bv_C[word] |= bit;
            if (mask & 0b0100) genome.bv_G[word] |= bit;
            if (mask & 0b1000) genome.bv_T[word] |= bit;

            ++pos;
        }

        genome.chromosomes.push_back(std::move(chrom));
    }

    return genome;
}

// Convenience loader

Genome load_fasta(const std::string& filepath) {
    return encode_genome(parse_fasta(filepath));
}

// Sequence utilities

std::string reverse_complement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.size());

    // Process in reverse order
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        char c = *it;
        switch (c) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
            case 'a': rc += 't'; break;
            case 't': rc += 'a'; break;
            case 'c': rc += 'g'; break;
            case 'g': rc += 'c'; break;
            case 'N': rc += 'N'; break;
            case 'n': rc += 'n'; break;
            default:
                throw std::invalid_argument(
                    std::string("Invalid nucleotide for reverse complement: '") + c + "'");
        }
    }
    return rc;
}

std::string extract_genome_slice(GenomeView view, size_t start, size_t length) {
    if (start + length > view.total_bases) {
        throw std::out_of_range(
            "Genome slice [" + std::to_string(start) + ", " +
            std::to_string(start + length) + ") exceeds genome bounds (" +
            std::to_string(view.total_bases) + ")");
    }

    std::string result;
    result.reserve(length);

    for (size_t i = 0; i < length; ++i) {
        size_t pos = start + i;
        size_t word = pos / 64;
        uint64_t bit = 1ULL << (pos % 64);

        bool a = (view.bv_A[word] & bit) != 0;
        bool c = (view.bv_C[word] & bit) != 0;
        bool g = (view.bv_G[word] & bit) != 0;
        bool t = (view.bv_T[word] & bit) != 0;

        // Decode: all four bits set = N, otherwise exactly one bit set
        if (a && c && g && t) {
            result += 'N';
        } else if (a) {
            result += 'A';
        } else if (c) {
            result += 'C';
        } else if (g) {
            result += 'G';
        } else if (t) {
            result += 'T';
        } else {
            // Should not happen if genome was encoded correctly
            result += 'N';
        }
    }

    return result;
}

// GenomeView creation

GenomeView make_view(const Genome& genome) {
    GenomeView view;
    view.bv_A = genome.bv_A.data();
    view.bv_C = genome.bv_C.data();
    view.bv_G = genome.bv_G.data();
    view.bv_T = genome.bv_T.data();
    view.chromosomes = &genome.chromosomes;
    view.total_bases = genome.total_bases;
    view.num_words = (genome.total_bases + 63) / 64;
    return view;
}

GenomeView make_view(const MappedGenome& genome) {
    GenomeView view;
    view.bv_A = genome.bv_A;
    view.bv_C = genome.bv_C;
    view.bv_G = genome.bv_G;
    view.bv_T = genome.bv_T;
    view.chromosomes = &genome.chromosomes;
    view.total_bases = genome.total_bases;
    view.num_words = genome.num_words;
    return view;
}

// MappedGenome RAII

MappedGenome::MappedGenome()
    : bv_A(nullptr), bv_C(nullptr), bv_G(nullptr), bv_T(nullptr),
      total_bases(0), num_words(0),
      mapped_data(nullptr), mapped_size(0), fd(-1) {}

MappedGenome::~MappedGenome() {
    if (mapped_data != nullptr && mapped_data != MAP_FAILED) {
        munmap(mapped_data, mapped_size);
    }
    if (fd >= 0) {
        close(fd);
    }
}

MappedGenome::MappedGenome(MappedGenome&& other) noexcept
    : bv_A(other.bv_A), bv_C(other.bv_C), bv_G(other.bv_G), bv_T(other.bv_T),
      chromosomes(std::move(other.chromosomes)),
      total_bases(other.total_bases), num_words(other.num_words),
      mapped_data(other.mapped_data), mapped_size(other.mapped_size), fd(other.fd) {
    // Clear other's resources so its destructor doesn't free them
    other.bv_A = nullptr;
    other.bv_C = nullptr;
    other.bv_G = nullptr;
    other.bv_T = nullptr;
    other.mapped_data = nullptr;
    other.mapped_size = 0;
    other.fd = -1;
}

MappedGenome& MappedGenome::operator=(MappedGenome&& other) noexcept {
    if (this != &other) {
        // Clean up existing resources
        if (mapped_data != nullptr && mapped_data != MAP_FAILED) {
            munmap(mapped_data, mapped_size);
        }
        if (fd >= 0) {
            close(fd);
        }

        // Move from other
        bv_A = other.bv_A;
        bv_C = other.bv_C;
        bv_G = other.bv_G;
        bv_T = other.bv_T;
        chromosomes = std::move(other.chromosomes);
        total_bases = other.total_bases;
        num_words = other.num_words;
        mapped_data = other.mapped_data;
        mapped_size = other.mapped_size;
        fd = other.fd;

        // Clear other
        other.bv_A = nullptr;
        other.bv_C = nullptr;
        other.bv_G = nullptr;
        other.bv_T = nullptr;
        other.mapped_data = nullptr;
        other.mapped_size = 0;
        other.fd = -1;
    }
    return *this;
}

// Index file format utilities

bool is_stomata_index(const std::string& filepath) {
    return filepath.size() >= 3 &&
           filepath.compare(filepath.size() - 3, 3, ".st") == 0;
}

// Header layout (64 bytes total):
//   0-7:   Magic bytes "STOMATA1"
//   8-15:  Version (uint64_t)
//  16-23:  total_bases (uint64_t)
//  24-31:  num_chromosomes (uint64_t)
//  32-39:  chromosome_table_offset (uint64_t)
//  40-47:  bitvector_offset (uint64_t)
//  48-63:  Reserved (16 bytes)

constexpr size_t HEADER_SIZE = 64;
constexpr size_t ALIGNMENT = 64;  // Align bit-vectors to 64-byte boundary

void write_genome_index(const Genome& genome, const std::string& filepath) {
    std::ofstream file(filepath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create index file: " + filepath);
    }

    size_t num_words = (genome.total_bases + 63) / 64;
    size_t num_chromosomes = genome.chromosomes.size();

    // Calculate chromosome table size
    size_t chrom_table_size = 0;
    for (const auto& chrom : genome.chromosomes) {
        chrom_table_size += 4;                    // name_length (uint32_t)
        chrom_table_size += chrom.name.size();    // name (NOT null-terminated)
        chrom_table_size += 8;                    // start (uint64_t)
        chrom_table_size += 8;                    // length (uint64_t)
    }

    // Offsets
    uint64_t chrom_table_offset = HEADER_SIZE;
    uint64_t bitvector_offset = chrom_table_offset + chrom_table_size;
    // Align bitvector_offset to 64-byte boundary
    bitvector_offset = (bitvector_offset + ALIGNMENT - 1) & ~(ALIGNMENT - 1);

    // Write header
    file.write(STOMATA_INDEX_MAGIC, 8);
    file.write(reinterpret_cast<const char*>(&STOMATA_INDEX_VERSION), 8);
    file.write(reinterpret_cast<const char*>(&genome.total_bases), 8);
    file.write(reinterpret_cast<const char*>(&num_chromosomes), 8);
    file.write(reinterpret_cast<const char*>(&chrom_table_offset), 8);
    file.write(reinterpret_cast<const char*>(&bitvector_offset), 8);
    // Reserved padding (16 bytes)
    char reserved[16] = {0};
    file.write(reserved, 16);

    // Write chromosome table
    for (const auto& chrom : genome.chromosomes) {
        uint32_t name_len = static_cast<uint32_t>(chrom.name.size());
        file.write(reinterpret_cast<const char*>(&name_len), 4);
        file.write(chrom.name.data(), chrom.name.size());
        file.write(reinterpret_cast<const char*>(&chrom.start), 8);
        file.write(reinterpret_cast<const char*>(&chrom.length), 8);
    }

    // Padding to align bit-vectors
    size_t current_pos = chrom_table_offset + chrom_table_size;
    size_t padding = bitvector_offset - current_pos;
    if (padding > 0) {
        std::vector<char> pad(padding, 0);
        file.write(pad.data(), padding);
    }

    // Write bit-vectors contiguously
    file.write(reinterpret_cast<const char*>(genome.bv_A.data()), num_words * 8);
    file.write(reinterpret_cast<const char*>(genome.bv_C.data()), num_words * 8);
    file.write(reinterpret_cast<const char*>(genome.bv_G.data()), num_words * 8);
    file.write(reinterpret_cast<const char*>(genome.bv_T.data()), num_words * 8);

    if (!file.good()) {
        throw std::runtime_error("Error writing index file: " + filepath);
    }
}

MappedGenome load_genome_index(const std::string& filepath) {
    MappedGenome result;

    // Open file
    result.fd = open(filepath.c_str(), O_RDONLY);
    if (result.fd < 0) {
        throw std::runtime_error("Cannot open index file: " + filepath);
    }

    // Get file size
    struct stat st;
    if (fstat(result.fd, &st) < 0) {
        close(result.fd);
        throw std::runtime_error("Cannot stat index file: " + filepath);
    }
    result.mapped_size = static_cast<size_t>(st.st_size);

    // Memory map the file
    result.mapped_data = mmap(nullptr, result.mapped_size, PROT_READ, MAP_PRIVATE,
                               result.fd, 0);
    if (result.mapped_data == MAP_FAILED) {
        close(result.fd);
        throw std::runtime_error("Cannot mmap index file: " + filepath);
    }

    const char* data = static_cast<const char*>(result.mapped_data);

    // Validate header
    if (result.mapped_size < HEADER_SIZE) {
        throw std::runtime_error("Index file too small: " + filepath);
    }

    // Check magic bytes
    if (std::memcmp(data, STOMATA_INDEX_MAGIC, 8) != 0) {
        throw std::runtime_error("Invalid index file magic: " + filepath);
    }

    // Read header fields
    uint64_t version;
    std::memcpy(&version, data + 8, 8);
    if (version != STOMATA_INDEX_VERSION) {
        throw std::runtime_error("Unsupported index version " +
                                  std::to_string(version) + " in: " + filepath);
    }

    std::memcpy(&result.total_bases, data + 16, 8);

    uint64_t num_chromosomes;
    std::memcpy(&num_chromosomes, data + 24, 8);

    uint64_t chrom_table_offset;
    std::memcpy(&chrom_table_offset, data + 32, 8);

    uint64_t bitvector_offset;
    std::memcpy(&bitvector_offset, data + 40, 8);

    result.num_words = (result.total_bases + 63) / 64;

    // Parse chromosome table
    const char* chrom_ptr = data + chrom_table_offset;
    result.chromosomes.reserve(num_chromosomes);
    for (uint64_t i = 0; i < num_chromosomes; ++i) {
        Chromosome chrom;

        uint32_t name_len;
        std::memcpy(&name_len, chrom_ptr, 4);
        chrom_ptr += 4;

        chrom.name.assign(chrom_ptr, name_len);
        chrom_ptr += name_len;

        std::memcpy(&chrom.start, chrom_ptr, 8);
        chrom_ptr += 8;

        std::memcpy(&chrom.length, chrom_ptr, 8);
        chrom_ptr += 8;

        result.chromosomes.push_back(std::move(chrom));
    }

    // Set bit-vector pointers (zero-copy access into mmap'd region)
    const uint64_t* bv_base = reinterpret_cast<const uint64_t*>(data + bitvector_offset);
    result.bv_A = bv_base;
    result.bv_C = bv_base + result.num_words;
    result.bv_G = bv_base + 2 * result.num_words;
    result.bv_T = bv_base + 3 * result.num_words;

    return result;
}
