#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <genome_loader.hpp>   // Genome struct

// Result of a shift-add bit-parallel Hamming distance sweep.
// distances[i] is the minimum Hamming distance (substitutions only) between
// the pattern and the substring text[i - m + 1 .. i] (i.e., the best alignment
// ending at position i).  The first m-1 entries are undefined / not meaningful
// because a full-length alignment cannot end before position m-1.
//
// NOTE: Hamming distance ignores insertions and deletions, so it will typically
// be >= the Levenshtein distance computed by Myers' algorithm. Use this when
// you only care about substitutions (no indels/bulges).
struct ShiftAddResult {
    std::vector<uint8_t> distances;   // one entry per text position
    size_t pattern_len;               // m
    size_t text_len;                  // n
};

// Compute semi-global Hamming distances (substitution-only) between a short
// pattern and a long text using a bit-parallel bitap (Baeza-Yates–Gonnet).
// Both strings must contain only [ACGTNacgtn].
// N in pattern is a wildcard (matches any nucleotide at zero cost).
// N in text is masked (counts as a mismatch against any pattern base).
//
// Substitutions only — insertions and deletions are not modeled, so the
// reported distance is typically >= the Levenshtein distance from Myers.
//
// Throws std::invalid_argument if either string is empty or contains an
// invalid character.
ShiftAddResult shift_add_hamming_distance(const std::string& pattern,
                                          const std::string& text);

// Bit-vector overload: sweep the pattern over a contiguous slice of a
// pre-encoded Genome via GenomeView.  'start' and 'length' define the slice
// within the genome's bit-vectors (must satisfy start + length <= view.total_bases).
// This is the overload that the GPU kernel will mirror.
//
// Throws std::invalid_argument if pattern is empty or contains an invalid
// character, or if the slice is out of range.
ShiftAddResult shift_add_hamming_distance_bv(const std::string& pattern,
                                             GenomeView         view,
                                             size_t             start,
                                             size_t             length);

// Convenience overload for Genome (creates view internally).
inline ShiftAddResult shift_add_hamming_distance_bv(const std::string& pattern,
                                                    const Genome&      genome,
                                                    size_t             start,
                                                    size_t             length) {
    return shift_add_hamming_distance_bv(pattern, make_view(genome), start, length);
}
