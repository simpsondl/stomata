#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <genome_loader.hpp>   // Genome struct

// Default chunk stride used by the GPU kernel.  One thread processes
// 'stride' output positions plus an overlap warm-up prefix.  Exported here
// so tests can place exact matches precisely at chunk boundaries.
static constexpr size_t GPU_DEFAULT_STRIDE = 65536;

// Result of a GPU Myers bit-parallel edit-distance sweep.
// Semantics are identical to MysersResult (myers_cpu.hpp):
// distances[i] is the minimum semi-global edit distance between the pattern
// and any substring of the genome slice ending at position (start + i).
struct GpuMysersResult {
    std::vector<uint8_t> distances;   // one entry per position in [start, start+length)
    size_t pattern_len;               // m
    size_t text_len;                  // length (== distances.size())
};

// A single sparse hit from GPU output (only positions with distance <= threshold).
struct GpuSparseHit {
    size_t  position;   // Absolute position in genome
    uint8_t distance;   // Edit distance at this position
};

// Sparse result for a single pattern from multi-pattern GPU search.
struct GpuSparseResult {
    std::vector<GpuSparseHit> hits;
    size_t pattern_len;
    size_t text_len;    // Total genome length searched
};

// Returns true if at least one CUDA-capable device is available.
bool gpu_available();

// Myers' bit-parallel edit distance on the GPU.
// The genome is tiled into GPU_DEFAULT_STRIDE chunks; each chunk runs on a
// separate thread with a 2*m warm-up overlap, which guarantees bit-identical
// results to myers_edit_distances_bv() (semi-global, free start in text).
// N in the pattern is a wildcard; N in the genome matches nothing.
// Constraints: pattern length ∈ [1, 64]; pattern chars ∈ [ACGTNacgtn].
GpuMysersResult myers_gpu(const std::string& pattern,
                          GenomeView         view,
                          size_t             start,
                          size_t             length);

// Convenience overload for Genome (creates view internally).
inline GpuMysersResult myers_gpu(const std::string& pattern,
                                 const Genome&      genome,
                                 size_t             start,
                                 size_t             length) {
    return myers_gpu(pattern, make_view(genome), start, length);
}

// Batched Myers GPU search: returns positions with edit distance <= threshold.
// The distances vector is expanded to full length (non-hit positions = 255) for
// backwards-compatible dense access. For sparse output, use myers_gpu_batch_sparse.
// Constraints: same as myers_gpu() — each pattern length ∈ [1, 64].
std::vector<GpuMysersResult> myers_gpu_batch_threshold(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold);

// Sparse batched Myers GPU search: returns only hits (no dense distance arrays).
// Recommended API for multi-pattern CRISPR off-target scans.
// The genome is transferred to the GPU once and up to 32 patterns run per kernel
// launch; per-pattern results come back in the same order as `patterns`.
std::vector<GpuSparseResult> myers_gpu_batch_sparse(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold);

inline std::vector<GpuSparseResult> myers_gpu_batch_sparse(
    const std::vector<std::string>& patterns,
    const Genome&                   genome,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    return myers_gpu_batch_sparse(patterns, make_view(genome), start, length, threshold);
}

inline std::vector<GpuMysersResult> myers_gpu_batch_threshold(
    const std::vector<std::string>& patterns,
    const Genome&                   genome,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    return myers_gpu_batch_threshold(patterns, make_view(genome), start, length, threshold);
}

// GPU Shift-Add (substitution-only Hamming distance).
// distances[i] is the minimum Hamming distance between the pattern and any
// substring of the genome slice ending at position (start + i).
struct GpuShiftAddResult {
    std::vector<uint8_t> distances;   // one entry per position in [start, start+length)
    size_t pattern_len;               // m
    size_t text_len;                  // length (== distances.size())
};

// Shift-add bit-parallel Hamming distance on the GPU.
// Semi-global alignment with free start in text. N in pattern is a wildcard;
// N in genome (masked position) matches nothing.
// Constraints: pattern length ∈ [1, 64]; pattern chars ∈ [ACGTNacgtn].
GpuShiftAddResult shift_add_gpu(const std::string& pattern,
                                GenomeView         view,
                                size_t             start,
                                size_t             length);

// Convenience overload for Genome (creates view internally).
inline GpuShiftAddResult shift_add_gpu(const std::string& pattern,
                                       const Genome&      genome,
                                       size_t             start,
                                       size_t             length) {
    return shift_add_gpu(pattern, make_view(genome), start, length);
}

// Batched shift-add with threshold: dense output (non-hit positions = 255)
// for backwards compatibility. Use shift_add_gpu_batch_sparse for sparse hits.
std::vector<GpuShiftAddResult> shift_add_gpu_batch_threshold(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold);

// Convenience overload for Genome.
inline std::vector<GpuShiftAddResult> shift_add_gpu_batch_threshold(
    const std::vector<std::string>& patterns,
    const Genome&                   genome,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    return shift_add_gpu_batch_threshold(patterns, make_view(genome), start, length, threshold);
}

// Sparse batched shift-add: returns only hits. Use when only substitutions matter.
std::vector<GpuSparseResult> shift_add_gpu_batch_sparse(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold);

// Convenience overload for Genome.
inline std::vector<GpuSparseResult> shift_add_gpu_batch_sparse(
    const std::vector<std::string>& patterns,
    const Genome&                   genome,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    return shift_add_gpu_batch_sparse(patterns, make_view(genome), start, length, threshold);
}

