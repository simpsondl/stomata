#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Bit-width of the single-word Myers / shift-add bitap lane.
// Kernels and inner loops are written against this width and a pattern
// length limit of MAX_SINGLE_WORD_PATTERN.
constexpr size_t MYERS_WORD_BITS = 64;
constexpr size_t MAX_SINGLE_WORD_PATTERN = 64;

// Returns pattern upper-cased, rejecting anything outside ACGTN.
// 'label' is prepended to error messages to identify which argument failed.
std::string validate_and_upper(const std::string& s, const char* label);

// Nucleotide index used by shift-add and verifier paths. Returns -1 for N.
int nuc_index(char c);

// Number of 64-bit words needed to hold 'bits' bits.
constexpr size_t num_words(size_t bits) {
    return (bits + MYERS_WORD_BITS - 1) / MYERS_WORD_BITS;
}

// Per-nucleotide pattern-match bit-vectors PM[0..3] indexed A=0 C=1 G=2 T=3.
// Bit i of PM[k] is set iff pattern[i] == nucleotide k OR pattern[i] == 'N'.
// N in the pattern is a wildcard (sets all four vectors).
void build_pattern_masks(const std::string& pattern,
                         std::vector<uint64_t> PM[4]);
