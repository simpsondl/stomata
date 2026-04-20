#include <shift_add_cpu.hpp>
#include <bio_utils.hpp>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Bit-parallel bitap for Hamming distance (Baeza-Yates–Gonnet).
//
// R[d] is a bitmask over pattern positions: bit i is set iff pattern[0..i]
// matches some suffix of text[0..j] with exactly d substitutions. At each
// text position we advance R[d] with:
//   R_d[j] = ((R_d[j-1] << 1) & char_mask)          // match
//          | (R_{d-1}[j-1] << 1)                     // substitution
//          | 1                                       // start new match
// The bit at position (m-1) of R[d] signals a full-pattern match ending at j.
// For j < m-1 we inspect bit j — the start is forced to text[0] so the result
// matches the naive prefix-vs-prefix distance.
//
// Single-word implementation: pattern length capped at 64. This covers every
// CRISPR spacer we see in practice (18–24 nt).

namespace {

// Advance the bitap state for one text position and return the minimum
// distance at which the pattern matches.
inline uint8_t bitap_step(uint64_t* R, int max_errors, uint64_t c_mask,
                          size_t check_bit, uint8_t default_dist) {
    uint64_t old_R = R[0];
    R[0] = ((old_R << 1) | 1) & c_mask;

    for (int d = 1; d <= max_errors; ++d) {
        uint64_t temp = R[d];
        R[d] = ((temp << 1) & c_mask) | (old_R << 1) | 1;
        old_R = temp;
    }

    const uint64_t check_mask = 1ULL << check_bit;
    for (int d = 0; d <= max_errors; ++d) {
        if (R[d] & check_mask) return static_cast<uint8_t>(d);
    }
    return default_dist;
}

// Bitap driver. `nuc_at(j)` returns 0..3 for A/C/G/T or -1 for N (which is
// treated as a non-match against any pattern base).
template<typename NucFn>
ShiftAddResult run_bitap(size_t m, size_t n,
                         const std::vector<uint64_t> PM[4],
                         NucFn nuc_at) {
    ShiftAddResult result;
    result.pattern_len = m;
    result.text_len    = n;
    result.distances.assign(n, static_cast<uint8_t>(m));
    if (m == 0 || n == 0) return result;

    const int max_errors = static_cast<int>(m);
    uint64_t R[MAX_SINGLE_WORD_PATTERN + 1] = {0};

    for (size_t j = 0; j < n; ++j) {
        const int nuc = nuc_at(j);
        const uint64_t c_mask = (nuc >= 0) ? PM[nuc][0] : 0;
        const size_t check_bit = (j < m - 1) ? j : (m - 1);
        result.distances[j] = bitap_step(R, max_errors, c_mask, check_bit,
                                         static_cast<uint8_t>(m));
    }
    return result;
}

// Nucleotide lookup for text strings: A/C/G/T (case-insensitive) → 0..3,
// anything else (N or invalid — callers have pre-validated) → -1.
inline int nuc_from_char(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}

// Decode a nucleotide from the genome bit-vectors at absolute position idx.
// Returns 0..3 for A/C/G/T, or -1 for N / invalid encoding.
inline int decode_nucleotide(const GenomeView& view, size_t idx) {
    const size_t word_idx = idx / MYERS_WORD_BITS;
    const size_t bit_pos  = idx % MYERS_WORD_BITS;
    const uint64_t bit_mask = 1ULL << bit_pos;

    uint8_t enc = 0;
    if (view.bv_A[word_idx] & bit_mask) enc |= 0b0001;
    if (view.bv_C[word_idx] & bit_mask) enc |= 0b0010;
    if (view.bv_G[word_idx] & bit_mask) enc |= 0b0100;
    if (view.bv_T[word_idx] & bit_mask) enc |= 0b1000;

    switch (enc) {
        case 0b0001: return 0;
        case 0b0010: return 1;
        case 0b0100: return 2;
        case 0b1000: return 3;
        default:     return -1;  // N or masked
    }
}

void require_single_word(const std::string& pattern, const char* fn) {
    if (pattern.size() > MAX_SINGLE_WORD_PATTERN) {
        throw std::invalid_argument(
            std::string(fn) + ": pattern length > " +
            std::to_string(MAX_SINGLE_WORD_PATTERN) + " not yet supported");
    }
}

}  // namespace

ShiftAddResult shift_add_hamming_distance(const std::string& pattern,
                                          const std::string& text) {
    if (pattern.empty()) {
        throw std::invalid_argument("shift_add_hamming_distance: pattern is empty");
    }
    if (text.empty()) {
        throw std::invalid_argument("shift_add_hamming_distance: text is empty");
    }

    std::string p_upper = validate_and_upper(pattern, "pattern");
    std::string t_upper = validate_and_upper(text, "text");

    require_single_word(p_upper, "shift_add_hamming_distance");

    std::vector<uint64_t> PM[4];
    build_pattern_masks(p_upper, PM);

    const size_t n = t_upper.size();
    return run_bitap(p_upper.size(), n, PM,
                     [&](size_t j) { return nuc_from_char(t_upper[j]); });
}

ShiftAddResult shift_add_hamming_distance_bv(const std::string& pattern,
                                             GenomeView         view,
                                             size_t             start,
                                             size_t             length) {
    if (pattern.empty()) {
        throw std::invalid_argument("shift_add_hamming_distance_bv: pattern is empty");
    }
    if (length == 0) {
        throw std::invalid_argument("shift_add_hamming_distance_bv: length is zero");
    }
    if (start + length > view.total_bases) {
        throw std::invalid_argument(
            "shift_add_hamming_distance_bv: slice [" + std::to_string(start) +
            ", " + std::to_string(start + length) + ") exceeds genome size " +
            std::to_string(view.total_bases));
    }

    std::string p_upper = validate_and_upper(pattern, "pattern");
    require_single_word(p_upper, "shift_add_hamming_distance_bv");

    std::vector<uint64_t> PM[4];
    build_pattern_masks(p_upper, PM);

    return run_bitap(p_upper.size(), length, PM,
                     [&](size_t j) { return decode_nucleotide(view, start + j); });
}
