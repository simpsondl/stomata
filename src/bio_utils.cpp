#include <bio_utils.hpp>

#include <cctype>
#include <stdexcept>

std::string validate_and_upper(const std::string& s, const char* label) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        switch (u) {
            case 'A': case 'C': case 'G': case 'T': case 'N':
                out += u;
                break;
            default:
                throw std::invalid_argument(
                    std::string(label) + " contains invalid character: '" + c + "'");
        }
    }
    return out;
}

int nuc_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;   // N
    }
}

void build_pattern_masks(const std::string& pattern,
                         std::vector<uint64_t> PM[4]) {
    size_t w = num_words(pattern.size());
    for (int k = 0; k < 4; ++k) PM[k].assign(w, 0);

    for (size_t i = 0; i < pattern.size(); ++i) {
        size_t   word = i / MYERS_WORD_BITS;
        uint64_t bit  = 1ULL << (i % MYERS_WORD_BITS);
        switch (pattern[i]) {
            case 'A': PM[0][word] |= bit; break;
            case 'C': PM[1][word] |= bit; break;
            case 'G': PM[2][word] |= bit; break;
            case 'T': PM[3][word] |= bit; break;
            case 'N':
                PM[0][word] |= bit;
                PM[1][word] |= bit;
                PM[2][word] |= bit;
                PM[3][word] |= bit;
                break;
        }
    }
}
