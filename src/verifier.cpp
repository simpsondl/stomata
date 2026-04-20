#include <verifier.hpp>
#include <bio_utils.hpp>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

VerifierResult verify_edit_distances(const std::string& pattern,
                                     const std::string& text) {
    if (pattern.empty())
        throw std::invalid_argument("verify_edit_distances: pattern must not be empty");
    if (text.empty())
        throw std::invalid_argument("verify_edit_distances: text must not be empty");

    std::string pat = validate_and_upper(pattern, "pattern");
    std::string txt = validate_and_upper(text,    "text");

    size_t m = pat.size();
    size_t n = txt.size();

    // Two-row DP.  prev = row i-1, curr = row i.
    //   dp[0][j] = 0   for all j   (semi-global: free start in text)
    //   dp[i][0] = i               (aligning pattern[0..i-1] vs empty text)
    //   dp[i][j] = min( dp[i-1][j]   + 1,          // deletion from pattern
    //                   dp[i][j-1]   + 1,          // insertion in text
    //                   dp[i-1][j-1] + cost )      // match (0) / sub (1)
    // uint16_t cells avoid any truncation during the +1 additions.
    std::vector<uint16_t> prev(n + 1, 0);   // row 0: all zeroes
    std::vector<uint16_t> curr(n + 1, 0);

    for (size_t i = 1; i <= m; ++i) {
        curr[0] = static_cast<uint16_t>(i);   // gap in text = i deletions

        for (size_t j = 1; j <= n; ++j) {
            // Text N is masked: does NOT match ANY pattern nucleotide (including pattern N).
            // Pattern N is wildcard: matches any text nucleotide A/C/G/T (but not text N).
            uint16_t cost;
            if (txt[j-1] == 'N') {
                cost = 1;  // Text masked N doesn't match anything (even pattern N)
            } else if (pat[i-1] == 'N') {
                cost = 0;  // Pattern wildcard matches any A/C/G/T in text
            } else {
                cost = (pat[i-1] == txt[j-1]) ? 0 : 1;
            }

            curr[j] = std::min({
                static_cast<uint16_t>(prev[j]   + 1),   // deletion
                static_cast<uint16_t>(curr[j-1] + 1),   // insertion
                static_cast<uint16_t>(prev[j-1] + cost) // match / sub
            });
        }

        std::swap(prev, curr);
    }

    // After the loop prev holds row m (the final row).
    VerifierResult result;
    result.pattern_len = m;
    result.text_len    = n;
    result.distances.resize(n);
    for (size_t j = 0; j < n; ++j) {
        result.distances[j] = static_cast<uint8_t>(prev[j + 1]);
    }
    return result;
}
