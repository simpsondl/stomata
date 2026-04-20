#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <genome_loader.hpp>   // Genome struct

// Result of a Myers bit-parallel edit-distance sweep.
// distances[i] is the minimum edit distance between the pattern and the
// substring text[i - m + 1 .. i]  (i.e. the best alignment *ending* at
// position i).  The first m-1 entries are undefined / not meaningful because
// a full-length alignment cannot end before position m-1.
struct MysersResult {
    std::vector<uint8_t> distances;   // one entry per text position
    size_t pattern_len;               // m
    size_t text_len;                  // n
};

// Compute semi-global edit distances between a short pattern and a long text
// using Myers' O(nm/w) bit-parallel algorithm (w = 64 on this platform).
// Both strings must contain only [ACGTNacgtn].
// N in pattern is a wildcard (matches any nucleotide).
// N in text is a wildcard (matches any nucleotide).
// Note: For reference genome searches, use myers_edit_distances_bv which
// correctly handles masked Ns in the genome (encoded as 0b0000).
// Throws std::invalid_argument if either string is empty or contains an
// invalid character.
MysersResult myers_edit_distances(const std::string& pattern,
                                  const std::string& text);

// Bit-vector overload: sweep the pattern over a contiguous slice of a
// pre-encoded Genome via GenomeView.  'start' and 'length' define the slice
// within the genome's bit-vectors (must satisfy start + length <= view.total_bases).
// This is the overload that the GPU kernel will mirror.
// Throws std::invalid_argument if pattern is empty or contains an invalid
// character, or if the slice is out of range.
MysersResult myers_edit_distances_bv(const std::string& pattern,
                                     GenomeView         view,
                                     size_t             start,
                                     size_t             length);

// Convenience overload for Genome (creates view internally).
inline MysersResult myers_edit_distances_bv(const std::string& pattern,
                                            const Genome&      genome,
                                            size_t             start,
                                            size_t             length) {
    return myers_edit_distances_bv(pattern, make_view(genome), start, length);
}
