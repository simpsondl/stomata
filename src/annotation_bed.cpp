#include "annotation_bed.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace annot {

BedForest load_bed(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open BED file: " + path);
    }

    BedForest forest;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.size() >= 5 && line.substr(0, 5) == "track") continue;
        if (line.size() >= 7 && line.substr(0, 7) == "browser") continue;

        std::istringstream iss(line);
        std::string chrom;
        size_t start, end;
        if (!(iss >> chrom >> start >> end)) continue;

        forest[chrom].push_back({start, end});
    }

    for (auto& [chrom, intervals] : forest) {
        std::sort(intervals.begin(), intervals.end(),
                  [](const BedInterval& a, const BedInterval& b) {
                      return a.start < b.start;
                  });
    }

    return forest;
}

bool overlaps_any(const std::string& chrom, size_t query_start, size_t query_end,
                  const BedForest& forest) {
    auto it = forest.find(chrom);
    if (it == forest.end()) return false;

    const auto& intervals = it->second;
    if (intervals.empty()) return false;

    // All intervals with start < query_end could overlap.
    // Find the first with start >= query_end (upper bound on start).
    auto upper = std::lower_bound(
        intervals.begin(), intervals.end(), query_end,
        [](const BedInterval& iv, size_t val) { return iv.start < val; });

    // Among [begin, upper), check whether any has end > query_start.
    for (auto jt = intervals.begin(); jt != upper; ++jt) {
        if (jt->end > query_start) return true;
    }
    return false;
}

} // namespace annot
