#include <myers_cpu.hpp>
#include <bio_utils.hpp>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Myers bit-parallel core (single-word, m <= 64)

static MysersResult myers_single_word(const std::string& pattern,
                                      const std::string& text,
                                      const std::vector<uint64_t> PM[4]) {
    size_t m = pattern.size();
    size_t n = text.size();

    // Mask to isolate the lower m bits
    uint64_t top_bit = 1ULL << (m - 1);
    uint64_t mask    = (m == 64) ? ~0ULL : ((1ULL << m) - 1);

    uint64_t Pv = mask;   // all vertical deltas positive initially
    uint64_t Mv = 0;      // no negative vertical deltas
    int      score = static_cast<int>(m);   // initial score = m (all deletions)

    MysersResult result;
    result.pattern_len = m;
    result.text_len    = n;
    result.distances.resize(n);

    for (size_t j = 0; j < n; ++j) {
        int idx = nuc_index(text[j]);
        // N in reference genome does NOT match any pattern position (masked).
        // This is critical for genome sequence analysis: N = unknown base.
        // Only user-provided patterns can have N as wildcard (handled in build_pattern_masks).
        uint64_t Eq = (idx >= 0) ? PM[idx][0] : 0;

        // Myers' bit-parallel update equations
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq | Mv;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;

        // Update score based on the top bit of Ph and Mh
        if (Ph & top_bit) ++score;
        if (Mh & top_bit) --score;

        // Shift Ph and Mh left by 1.  Semi-global (free start in text):
        // H[0][j] = 0, so inject 0 at bit 0.  (Injecting 1 would be global.)
        Ph = (Ph << 1) & mask;
        Mh = (Mh << 1) & mask;

        // Update vertical delta vectors; mask Pv because ~(Xv|Ph) sets upper bits.
        Pv = (Mh | ~(Xv | Ph)) & mask;
        Mv = Ph & Xv;

        result.distances[j] = static_cast<uint8_t>(score);
    }

    return result;
}

// Myers bit-parallel core (multi-word, m > 64)

static MysersResult myers_multi_word(const std::string& pattern,
                                           const std::string& text,
                                           const std::vector<uint64_t> PM[4]) {
    size_t m = pattern.size();
    size_t n = text.size();
    size_t w = num_words(m);

    size_t top_bits  = m % 64;
    uint64_t top_mask = (top_bits == 0) ? ~0ULL : ((1ULL << top_bits) - 1);
    uint64_t top_bit  = 1ULL << ((top_bits == 0 ? 64 : top_bits) - 1);

    std::vector<uint64_t> Pv(w, ~0ULL);
    std::vector<uint64_t> Mv(w, 0ULL);
    Pv[w - 1] &= top_mask;

    int score = static_cast<int>(m);

    MysersResult result;
    result.pattern_len = m;
    result.text_len    = n;
    result.distances.resize(n);

    for (size_t j = 0; j < n; ++j) {
        int idx = nuc_index(text[j]);
        // N in reference genome does NOT match any pattern position (masked).
        // This is critical for genome sequence analysis: N = unknown base.
        // Only user-provided patterns can have N as wildcard (handled in build_pattern_masks).

        uint64_t carry_add = 0;
        uint64_t carry_Ph  = 0;   // semi-global: H[0][j]=0, inject 0 at bit 0
        uint64_t carry_Mh  = 0;
        uint64_t Ph_top    = 0;   // will hold the top bit of Ph (pre-shift) for score
        uint64_t Mh_top    = 0;

        for (size_t k = 0; k < w; ++k) {
            // N in reference genome does NOT match any pattern position (masked).
            uint64_t Eq = (idx >= 0) ? PM[idx][k] : 0;
            uint64_t pv = Pv[k];
            uint64_t mv = Mv[k];

            uint64_t Xv = Eq | mv;

            // (Eq & Pv) + Pv  with carry propagation
            uint64_t eq_and_pv = Eq & pv;
            uint64_t sum       = eq_and_pv + pv;
            uint64_t c1        = (sum < eq_and_pv) ? 1ULL : 0ULL;
            uint64_t sum2      = sum + carry_add;
            uint64_t c2        = (sum2 < sum) ? 1ULL : 0ULL;
            carry_add          = c1 | c2;

            uint64_t Xh = (sum2 ^ pv) | Eq | mv;
            uint64_t Ph = mv | ~(Xh | pv);
            uint64_t Mh = pv & Xh;

            // Record top bit of Ph / Mh in the last word (before shift)
            if (k == w - 1) {
                Ph_top = (Ph >> ((top_bits == 0 ? 64 : top_bits) - 1)) & 1;
                Mh_top = (Mh >> ((top_bits == 0 ? 64 : top_bits) - 1)) & 1;
            }

            // Shift Ph left by 1 with carry
            uint64_t Ph_shifted = (Ph << 1) | carry_Ph;
            carry_Ph            = (Ph >> 63) & 1;

            // Shift Mh left by 1 with carry
            uint64_t Mh_shifted = (Mh << 1) | carry_Mh;
            carry_Mh            = (Mh >> 63) & 1;

            // Update vertical deltas using shifted horizontal deltas
            Pv[k] = Mh_shifted | ~(Xv | Ph_shifted);
            Mv[k] = Ph_shifted & Xv;

            // Mask top word
            if (k == w - 1) {
                Pv[k] &= top_mask;
                Mv[k] &= top_mask;
            }
        }

        // Update score
        if (Ph_top) ++score;
        if (Mh_top) --score;

        result.distances[j] = static_cast<uint8_t>(score);
    }

    return result;
}

// Bit-vector text extraction helper
// Note: extract_genome_slice is now defined in genome_loader.cpp
// Public API

MysersResult myers_edit_distances(const std::string& pattern,
                                  const std::string& text) {
    if (pattern.empty())
        throw std::invalid_argument("myers_edit_distances: pattern must not be empty");
    if (text.empty())
        throw std::invalid_argument("myers_edit_distances: text must not be empty");

    std::string pat = validate_and_upper(pattern, "pattern");
    std::string txt = validate_and_upper(text,    "text");

    std::vector<uint64_t> PM[4];
    build_pattern_masks(pat, PM);

    if (pat.size() <= 64)
        return myers_single_word(pat, txt, PM);
    else
        return myers_multi_word(pat, txt, PM);
}

MysersResult myers_edit_distances_bv(const std::string& pattern,
                                     GenomeView         view,
                                     size_t             start,
                                     size_t             length) {
    if (pattern.empty())
        throw std::invalid_argument("myers_edit_distances_bv: pattern must not be empty");
    if (start + length > view.total_bases)
        throw std::invalid_argument("myers_edit_distances_bv: slice [start, start+length) is out of range");
    if (length == 0)
        throw std::invalid_argument("myers_edit_distances_bv: length must be > 0");

    std::string pat = validate_and_upper(pattern, "pattern");
    size_t m = pat.size();

    std::vector<uint64_t> PM[4];
    build_pattern_masks(pat, PM);

    // For now, only implement single-word version
    if (m > 64) {
        // Fall back to string-based version for multi-word patterns
        std::string text = extract_genome_slice(view, start, length);
        return myers_edit_distances(pattern, text);
    }

    // Single-word implementation working directly on bit-vectors
    uint64_t top_bit = 1ULL << (m - 1);
    uint64_t mask    = (m == 64) ? ~0ULL : ((1ULL << m) - 1);

    uint64_t Pv = mask;
    uint64_t Mv = 0;
    int      score = static_cast<int>(m);

    MysersResult result;
    result.pattern_len = m;
    result.text_len    = length;
    result.distances.resize(length);

    for (size_t j = 0; j < length; ++j) {
        size_t abs_pos = start + j;
        size_t word    = abs_pos / 64;
        uint64_t bit   = 1ULL << (abs_pos % 64);

        // Extract nucleotide bits from genome view
        bool bA = view.bv_A[word] & bit;
        bool bC = view.bv_C[word] & bit;
        bool bG = view.bv_G[word] & bit;
        bool bT = view.bv_T[word] & bit;

        // Compute Eq: which pattern positions match this text character?
        // If all bits are 0 (masked N), Eq = 0 (matches nothing)
        uint64_t Eq;
        if (!bA && !bC && !bG && !bT) {
            Eq = 0;  // Masked N: matches no pattern positions
        } else {
            // Normal nucleotide: OR together the PM vectors for set bits
            Eq = 0;
            if (bA) Eq |= PM[0][0];
            if (bC) Eq |= PM[1][0];
            if (bG) Eq |= PM[2][0];
            if (bT) Eq |= PM[3][0];
        }

        // Myers' bit-parallel update equations
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq | Mv;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;

        // Update score based on the top bit of Ph and Mh
        if (Ph & top_bit) ++score;
        if (Mh & top_bit) --score;

        // Shift Ph and Mh left by 1.  Semi-global (free start in text):
        // H[0][j] = 0, so inject 0 at bit 0.
        Ph = (Ph << 1) & mask;
        Mh = (Mh << 1) & mask;

        // Update vertical delta vectors; mask Pv because ~(Xv|Ph) sets upper bits.
        Pv = (Mh | ~(Xv | Ph)) & mask;
        Mv = Ph & Xv;

        result.distances[j] = static_cast<uint8_t>(score);
    }

    return result;
}
