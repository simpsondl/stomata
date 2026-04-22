#pragma once

#include <map>
#include <string>
#include <vector>
#include <cstddef>

namespace annot {

struct BedInterval {
    size_t start;  // 0-based, inclusive
    size_t end;    // 0-based, exclusive
};

// Per-chromosome sorted interval list. Sorted by start for binary-search overlap.
using BedForest = std::map<std::string, std::vector<BedInterval>>;

// Load a BED file into a BedForest.
// Only chrom/start/end columns are used; track/browser lines and # comments are skipped.
// Throws std::runtime_error if the file cannot be opened.
BedForest load_bed(const std::string& path);

// Returns true if [query_start, query_end) overlaps any interval on chrom.
bool overlaps_any(const std::string& chrom, size_t query_start, size_t query_end,
                  const BedForest& forest);

} // namespace annot
