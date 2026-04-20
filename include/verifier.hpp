#pragma once

#include <cstdint>
#include <string>
#include <vector>

// Result of a brute-force O(mn) semi-global edit-distance computation.
// distances[j] is the minimum edit distance between the full pattern and
// any substring of the text that ends at position j.  The vector has one
// entry per text position; the first m-1 entries may be large (no full-
// length window fits) but are still well-defined.
struct VerifierResult {
    std::vector<uint8_t> distances;   // one entry per text position
    size_t pattern_len;               // m
    size_t text_len;                  // n
};

// Compute semi-global edit distances using a textbook O(mn) DP recurrence.
// This is the ground-truth reference implementation used to validate Myers'
// bit-parallel algorithm.
//
// Semi-global semantics: free gaps at the start of the text (the first row
// of the DP matrix is all zeroes).  Equivalently, distances[j] is the
// minimum edit distance between the pattern and *any* substring of the text
// that ends at position j.
//
// N is treated as a wildcard that matches any nucleotide (including another
// N) at zero cost.  Both upper- and lower-case input is accepted.
//
// Throws std::invalid_argument if either string is empty or contains a
// character outside [ACGTNacgtn].
VerifierResult verify_edit_distances(const std::string& pattern,
                                     const std::string& text);
