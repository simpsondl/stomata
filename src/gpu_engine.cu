#include <gpu_engine.hpp>
#include <bio_utils.hpp>

#include <cuda_runtime.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <atomic>

// Maximum patterns per multi-pattern kernel launch.
// Sized to fit pattern masks in shared memory: 32 × 4 × 8 = 1 KB.
static constexpr int MAX_PATTERNS_PER_KERNEL = 32;

// Maximum hits per pattern retained on-device before host compaction.
// Caps memory for promiscuous patterns that saturate the output buffer.
static constexpr size_t MAX_HITS_PER_PATTERN = 10'000'000;

// Sparse hit output entry (position + distance + pattern index).
struct SparseHit {
    uint32_t position;    // Relative position in genome slice
    uint8_t  distance;
    uint8_t  pattern_idx; // Which pattern in the batch matched (0-31)
    uint16_t padding;
};

// CUDA error helper

static void check_cuda(cudaError_t err, const char* msg) {
    if (err != cudaSuccess)
        throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
}

// RAII wrapper for device memory

struct DevicePtr {
    void* ptr = nullptr;
    DevicePtr() = default;
    ~DevicePtr() { if (ptr) cudaFree(ptr); }
    DevicePtr(const DevicePtr&) = delete;
    DevicePtr& operator=(const DevicePtr&) = delete;
};

// Kernel argument struct (passed by value; must stay under 4 096 bytes)

struct KernelArgs {
    const uint64_t* bv_A;         // device pointer — genome bit-vector A
    const uint64_t* bv_C;         //                          C
    const uint64_t* bv_G;         //                          G
    const uint64_t* bv_T;         //                          T
    uint8_t*        distances;    // device pointer — output array

    size_t genome_start;          // absolute start position in original genome
    size_t genome_len;            // number of positions to sweep (= length)
    size_t bv_word_offset;        // first uint64 word copied to device (for indexing)

    uint64_t PM[4];               // pattern-match masks  [A C G T]
    int      m;                   // pattern length  1 ≤ m ≤ 64

    size_t stride;                // non-overlapping output positions per thread
    size_t overlap;               // warm-up prefix length  (= 2 * m)
};

// Multi-pattern kernel arguments (sparse output)

struct MultiPatternArgs {
    const uint64_t* bv_A;
    const uint64_t* bv_C;
    const uint64_t* bv_G;
    const uint64_t* bv_T;

    size_t genome_start;
    size_t genome_len;
    size_t bv_word_offset;

    // Pattern masks: PM[pattern_idx][nucleotide] — flattened to [4 * pattern_idx + nuc]
    uint64_t PM[MAX_PATTERNS_PER_KERNEL * 4];
    int      pattern_lengths[MAX_PATTERNS_PER_KERNEL];
    int      num_patterns;
    int      threshold;           // Only output hits with distance <= threshold

    // Sparse output
    SparseHit* hits;              // Output buffer for sparse hits
    uint32_t*  hit_counts;        // Per-pattern hit counts (atomically incremented)
    size_t     max_hits_per_pattern;

    size_t stride;
    size_t overlap;
};

// Single-pattern Myers kernel: one thread per genome chunk.
// The 2*m warm-up overlap covers the maximum alignment context for a
// semi-global match, so every position in the output range sees identical
// Myers state to a single-thread sweep. Output is dense (one score per position).
__global__ void myers_kernel(KernelArgs args) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    size_t num_chunks = (args.genome_len + args.stride - 1) / args.stride;
    if (tid >= num_chunks) return;

    // Output range (relative to genome_start)
    size_t out_start = tid * args.stride;
    size_t out_end   = out_start + args.stride;
    if (out_end > args.genome_len) out_end = args.genome_len;

    // Processing range includes warm-up overlap (chunk 0 has none)
    size_t proc_start = (tid == 0) ? 0 : (out_start - args.overlap);
    size_t proc_end   = out_end;

    // Myers single-word state
    int      m       = args.m;
    uint64_t mask    = (m == 64) ? ~0ULL : ((1ULL << m) - 1);
    uint64_t top_bit = 1ULL << (m - 1);

    uint64_t Pv    = mask;
    uint64_t Mv    = 0ULL;
    int      score = m;

    // Cached word index and loaded bit-vector words — avoids redundant
    // global-memory loads for 64 consecutive positions in the same word.
    size_t   cached_word = ~(size_t)0;   // invalid sentinel
    uint64_t wA = 0, wC = 0, wG = 0, wT = 0;

    for (size_t rel = proc_start; rel < proc_end; ++rel) {
        size_t abs_pos = args.genome_start + rel;
        size_t  word   = abs_pos / 64;
        uint64_t bit   = 1ULL << (abs_pos % 64);

        // Load four bit-vector words once per 64-position block
        if (word != cached_word) {
            size_t dev_idx = word - args.bv_word_offset;
            wA = args.bv_A[dev_idx];
            wC = args.bv_C[dev_idx];
            wG = args.bv_G[dev_idx];
            wT = args.bv_T[dev_idx];
            cached_word = word;
        }

        // Decode nucleotide → Eq mask
        bool bA = wA & bit;
        bool bC = wC & bit;
        bool bG = wG & bit;
        bool bT = wT & bit;

        uint64_t Eq;
        // If all bits are 0 (masked N in reference genome), Eq = 0 (matches nothing).
        // This is critical: N in the genome sequence means unknown base and should NOT match anything.
        // N as wildcard only appears in user-provided patterns (handled in build_pattern_masks on host).
        if (!bA && !bC && !bG && !bT) {
            Eq = 0;  // Masked N: matches no pattern positions
        } else if (bA) {
            Eq = args.PM[0];
        } else if (bC) {
            Eq = args.PM[1];
        } else if (bG) {
            Eq = args.PM[2];
        } else {
            Eq = args.PM[3];    // T
        }

        // Myers bit-parallel update (single-word)
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq | Mv;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;

        // Score update from top-row horizontal deltas
        if (Ph & top_bit) ++score;
        if (Mh & top_bit) --score;

        // Shift horizontal deltas left; inject 0 at bit 0 (semi-global)
        Ph = (Ph << 1) & mask;
        Mh = (Mh << 1) & mask;

        // Update vertical deltas
        Pv = (Mh | ~(Xv | Ph)) & mask;
        Mv = Ph & Xv;

        // Write output only for positions inside this thread's output range
        if (rel >= out_start)
            args.distances[rel] = static_cast<uint8_t>(score);
    }
}

// Multi-pattern Myers kernel: up to 32 patterns per genome pass.
// Each thread iterates the genome chunk once and updates Myers state for
// every pattern in parallel. Output is sparse: hits (distance <= threshold)
// are pushed via per-pattern atomic counters, avoiding N dense output arrays.
__global__ void myers_multi_pattern_kernel(MultiPatternArgs args) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    size_t num_chunks = (args.genome_len + args.stride - 1) / args.stride;
    if (tid >= num_chunks) return;

    // Output range (relative to genome_start)
    size_t out_start = tid * args.stride;
    size_t out_end   = out_start + args.stride;
    if (out_end > args.genome_len) out_end = args.genome_len;

    // Processing range includes warm-up overlap (chunk 0 has none)
    // Use maximum pattern length for overlap to cover all patterns
    size_t proc_start = (tid == 0) ? 0 : (out_start - args.overlap);
    size_t proc_end   = out_end;

    // Per-pattern Myers state — stored in registers/local memory
    uint64_t Pv[MAX_PATTERNS_PER_KERNEL];
    uint64_t Mv[MAX_PATTERNS_PER_KERNEL];
    int      score[MAX_PATTERNS_PER_KERNEL];
    uint64_t mask[MAX_PATTERNS_PER_KERNEL];
    uint64_t top_bit[MAX_PATTERNS_PER_KERNEL];

    // Initialize all patterns
    for (int p = 0; p < args.num_patterns; ++p) {
        int m = args.pattern_lengths[p];
        mask[p]    = (m == 64) ? ~0ULL : ((1ULL << m) - 1);
        top_bit[p] = 1ULL << (m - 1);
        Pv[p]      = mask[p];
        Mv[p]      = 0ULL;
        score[p]   = m;
    }

    // Cached word index and loaded bit-vector words
    size_t   cached_word = ~(size_t)0;
    uint64_t wA = 0, wC = 0, wG = 0, wT = 0;

    for (size_t rel = proc_start; rel < proc_end; ++rel) {
        size_t abs_pos = args.genome_start + rel;
        size_t  word   = abs_pos / 64;
        uint64_t bit   = 1ULL << (abs_pos % 64);

        // Load bit-vector words once per 64-position block
        if (word != cached_word) {
            size_t dev_idx = word - args.bv_word_offset;
            wA = args.bv_A[dev_idx];
            wC = args.bv_C[dev_idx];
            wG = args.bv_G[dev_idx];
            wT = args.bv_T[dev_idx];
            cached_word = word;
        }

        // Decode nucleotide
        bool bA = wA & bit;
        bool bC = wC & bit;
        bool bG = wG & bit;
        bool bT = wT & bit;

        // Determine which nucleotide (or masked N)
        int nuc_idx = -1;  // -1 = masked N
        if (bA)      nuc_idx = 0;
        else if (bC) nuc_idx = 1;
        else if (bG) nuc_idx = 2;
        else if (bT) nuc_idx = 3;
        // else: masked N, nuc_idx stays -1

        // Update Myers state for ALL patterns
        for (int p = 0; p < args.num_patterns; ++p) {
            // Get Eq mask for this pattern and nucleotide
            uint64_t Eq = (nuc_idx >= 0) ? args.PM[p * 4 + nuc_idx] : 0ULL;

            // Myers bit-parallel update
            uint64_t Xv = Eq | Mv[p];
            uint64_t Xh = (((Eq & Pv[p]) + Pv[p]) ^ Pv[p]) | Eq | Mv[p];
            uint64_t Ph = Mv[p] | ~(Xh | Pv[p]);
            uint64_t Mh = Pv[p] & Xh;

            // Score update
            if (Ph & top_bit[p]) ++score[p];
            if (Mh & top_bit[p]) --score[p];

            // Shift horizontal deltas
            Ph = (Ph << 1) & mask[p];
            Mh = (Mh << 1) & mask[p];

            // Update vertical deltas
            Pv[p] = (Mh | ~(Xv | Ph)) & mask[p];
            Mv[p] = Ph & Xv;

            // Output sparse hit if within output range and passes threshold
            if (rel >= out_start && score[p] <= args.threshold) {
                // Atomically increment hit count and get write index
                uint32_t idx = atomicAdd(&args.hit_counts[p], 1);
                
                // Only write if within buffer capacity
                if (idx < args.max_hits_per_pattern) {
                    size_t write_pos = p * args.max_hits_per_pattern + idx;
                    args.hits[write_pos].position = static_cast<uint32_t>(rel);
                    args.hits[write_pos].distance = static_cast<uint8_t>(score[p]);
                    args.hits[write_pos].pattern_idx = static_cast<uint8_t>(p);
                }
            }
        }
    }
}

// Public API

bool gpu_available() {
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0);
}

GpuMysersResult myers_gpu(const std::string& pattern,
                          GenomeView         view,
                          size_t             start,
                          size_t             length) {
    // ── Input validation ──────────────────────────────────────────────────
    if (pattern.empty())
        throw std::invalid_argument("myers_gpu: pattern must not be empty");
    if (pattern.size() > 64)
        throw std::invalid_argument("myers_gpu: pattern length must be ≤ 64 (single-word)");
    if (length == 0)
        throw std::invalid_argument("myers_gpu: length must be > 0");
    if (start + length > view.total_bases)
        throw std::invalid_argument("myers_gpu: slice [start, start+length) is out of range");

    std::string pat = validate_and_upper(pattern, "pattern");
    size_t m = pat.size();

    // ── Build pattern-match masks (host) ──────────────────────────────────
    uint64_t PM[4] = {0, 0, 0, 0};   // A  C  G  T
    for (size_t i = 0; i < m; ++i) {
        uint64_t bit = 1ULL << i;
        switch (pat[i]) {
            case 'A': PM[0] |= bit; break;
            case 'C': PM[1] |= bit; break;
            case 'G': PM[2] |= bit; break;
            case 'T': PM[3] |= bit; break;
            case 'N': PM[0] |= bit; PM[1] |= bit;
                      PM[2] |= bit; PM[3] |= bit; break;
        }
    }

    // ── Determine which bit-vector words are needed ───────────────────────
    size_t first_word = start / 64;
    size_t last_word  = (start + length - 1) / 64 + 1;   // exclusive
    size_t num_words  = last_word - first_word;
    size_t bv_bytes   = num_words * sizeof(uint64_t);

    // ── Allocate device memory ────────────────────────────────────────────
    DevicePtr d_bvA, d_bvC, d_bvG, d_bvT, d_dist;

    check_cuda(cudaMalloc(&d_bvA.ptr, bv_bytes), "cudaMalloc bv_A");
    check_cuda(cudaMalloc(&d_bvC.ptr, bv_bytes), "cudaMalloc bv_C");
    check_cuda(cudaMalloc(&d_bvG.ptr, bv_bytes), "cudaMalloc bv_G");
    check_cuda(cudaMalloc(&d_bvT.ptr, bv_bytes), "cudaMalloc bv_T");

    size_t out_bytes = length * sizeof(uint8_t);
    check_cuda(cudaMalloc(&d_dist.ptr, out_bytes), "cudaMalloc distances");

    // ── Copy bit-vector slice to device ───────────────────────────────────
    check_cuda(cudaMemcpy(d_bvA.ptr, view.bv_A + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_A");
    check_cuda(cudaMemcpy(d_bvC.ptr, view.bv_C + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_C");
    check_cuda(cudaMemcpy(d_bvG.ptr, view.bv_G + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_G");
    check_cuda(cudaMemcpy(d_bvT.ptr, view.bv_T + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_T");

    // ── Sentinel-fill output (catches unwritten positions during debug) ───
    check_cuda(cudaMemset(d_dist.ptr, 0xFF, out_bytes), "cudaMemset distances");

    // ── Prepare kernel arguments ──────────────────────────────────────────
    size_t overlap    = 2 * m;
    size_t stride     = GPU_DEFAULT_STRIDE;
    size_t num_chunks = (length + stride - 1) / stride;

    KernelArgs args;
    args.bv_A            = static_cast<const uint64_t*>(d_bvA.ptr);
    args.bv_C            = static_cast<const uint64_t*>(d_bvC.ptr);
    args.bv_G            = static_cast<const uint64_t*>(d_bvG.ptr);
    args.bv_T            = static_cast<const uint64_t*>(d_bvT.ptr);
    args.distances       = static_cast<uint8_t*>(d_dist.ptr);
    args.genome_start    = start;
    args.genome_len      = length;
    args.bv_word_offset  = first_word;
    args.PM[0]           = PM[0];
    args.PM[1]           = PM[1];
    args.PM[2]           = PM[2];
    args.PM[3]           = PM[3];
    args.m               = static_cast<int>(m);
    args.stride          = stride;
    args.overlap         = overlap;

    // ── Launch kernel ─────────────────────────────────────────────────────
    dim3 block(256);
    dim3 grid((num_chunks + 255) / 256);

    myers_kernel<<<grid, block>>>(args);

    check_cuda(cudaGetLastError(),        "kernel launch");
    check_cuda(cudaDeviceSynchronize(),   "cudaDeviceSynchronize");

    // ── Copy results back to host ─────────────────────────────────────────
    GpuMysersResult result;
    result.pattern_len = m;
    result.text_len    = length;
    result.distances.resize(length);
    check_cuda(cudaMemcpy(result.distances.data(), d_dist.ptr,
                          out_bytes, cudaMemcpyDeviceToHost), "memcpy distances");

    return result;
}

// Batched GPU search with sparse output (multi-pattern kernel)

std::vector<GpuMysersResult> myers_gpu_batch_threshold(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    
    if (patterns.empty()) {
        return {};
    }
    
    // ── Validate inputs ───────────────────────────────────────────────────
    if (length == 0)
        throw std::invalid_argument("myers_gpu_batch: length must be > 0");
    if (start + length > view.total_bases)
        throw std::invalid_argument("myers_gpu_batch: slice out of range");
    if (threshold < 0)
        throw std::invalid_argument("myers_gpu_batch: threshold must be >= 0");
    
    // Validate and normalize all patterns upfront
    std::vector<std::string> normalized_patterns;
    normalized_patterns.reserve(patterns.size());
    size_t max_pattern_len = 0;
    
    for (size_t i = 0; i < patterns.size(); ++i) {
        if (patterns[i].empty())
            throw std::invalid_argument("myers_gpu_batch: pattern " + 
                std::to_string(i) + " is empty");
        if (patterns[i].size() > 64)
            throw std::invalid_argument("myers_gpu_batch: pattern " + 
                std::to_string(i) + " exceeds 64 bp");
        
        normalized_patterns.push_back(validate_and_upper(patterns[i], "pattern"));
        max_pattern_len = std::max(max_pattern_len, patterns[i].size());
    }
    
    // ── Allocate and copy genome bit-vectors to GPU (once) ────────────────
    size_t first_word = start / 64;
    size_t last_word  = (start + length - 1) / 64 + 1;
    size_t num_words  = last_word - first_word;
    size_t bv_bytes   = num_words * sizeof(uint64_t);
    
    DevicePtr d_bvA, d_bvC, d_bvG, d_bvT;
    check_cuda(cudaMalloc(&d_bvA.ptr, bv_bytes), "cudaMalloc bv_A");
    check_cuda(cudaMalloc(&d_bvC.ptr, bv_bytes), "cudaMalloc bv_C");
    check_cuda(cudaMalloc(&d_bvG.ptr, bv_bytes), "cudaMalloc bv_G");
    check_cuda(cudaMalloc(&d_bvT.ptr, bv_bytes), "cudaMalloc bv_T");
    
    check_cuda(cudaMemcpy(d_bvA.ptr, view.bv_A + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_A");
    check_cuda(cudaMemcpy(d_bvC.ptr, view.bv_C + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_C");
    check_cuda(cudaMemcpy(d_bvG.ptr, view.bv_G + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_G");
    check_cuda(cudaMemcpy(d_bvT.ptr, view.bv_T + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_T");
    
    // ── Prepare results vector ────────────────────────────────────────────
    std::vector<GpuMysersResult> results(patterns.size());
    for (size_t i = 0; i < patterns.size(); ++i) {
        results[i].pattern_len = normalized_patterns[i].size();
        results[i].text_len    = length;
        // Note: distances vector stays empty until we populate from sparse hits
    }
    
    // ── Process patterns in batches of MAX_PATTERNS_PER_KERNEL ────────────
    size_t stride = GPU_DEFAULT_STRIDE;
    size_t overlap = 2 * max_pattern_len;
    size_t num_chunks = (length + stride - 1) / stride;
    
    // Allocate sparse output buffers
    size_t max_hits = MAX_HITS_PER_PATTERN;
    
    for (size_t batch_start = 0; batch_start < patterns.size(); 
         batch_start += MAX_PATTERNS_PER_KERNEL) {
        
        size_t batch_end = std::min(batch_start + MAX_PATTERNS_PER_KERNEL, 
                                    patterns.size());
        int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);
        
        // Allocate sparse hit buffer and counters for this batch
        size_t hits_buffer_size = num_patterns_in_batch * max_hits * sizeof(SparseHit);
        size_t counters_size = num_patterns_in_batch * sizeof(uint32_t);
        
        DevicePtr d_hits, d_counters;
        check_cuda(cudaMalloc(&d_hits.ptr, hits_buffer_size), "cudaMalloc hits");
        check_cuda(cudaMalloc(&d_counters.ptr, counters_size), "cudaMalloc counters");
        
        // Zero the counters
        check_cuda(cudaMemset(d_counters.ptr, 0, counters_size), "memset counters");
        
        // Prepare multi-pattern kernel arguments
        MultiPatternArgs args;
        args.bv_A            = static_cast<const uint64_t*>(d_bvA.ptr);
        args.bv_C            = static_cast<const uint64_t*>(d_bvC.ptr);
        args.bv_G            = static_cast<const uint64_t*>(d_bvG.ptr);
        args.bv_T            = static_cast<const uint64_t*>(d_bvT.ptr);
        args.genome_start    = start;
        args.genome_len      = length;
        args.bv_word_offset  = first_word;
        args.num_patterns    = num_patterns_in_batch;
        args.threshold       = threshold;
        args.hits            = static_cast<SparseHit*>(d_hits.ptr);
        args.hit_counts      = static_cast<uint32_t*>(d_counters.ptr);
        args.max_hits_per_pattern = max_hits;
        args.stride          = stride;
        args.overlap         = overlap;
        
        // Build pattern masks for all patterns in batch
        std::memset(args.PM, 0, sizeof(args.PM));
        std::memset(args.pattern_lengths, 0, sizeof(args.pattern_lengths));
        
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            const std::string& pat = normalized_patterns[batch_start + p];
            args.pattern_lengths[p] = static_cast<int>(pat.size());
            
            for (size_t j = 0; j < pat.size(); ++j) {
                uint64_t bit = 1ULL << j;
                switch (pat[j]) {
                    case 'A': args.PM[p * 4 + 0] |= bit; break;
                    case 'C': args.PM[p * 4 + 1] |= bit; break;
                    case 'G': args.PM[p * 4 + 2] |= bit; break;
                    case 'T': args.PM[p * 4 + 3] |= bit; break;
                    case 'N': args.PM[p * 4 + 0] |= bit; 
                              args.PM[p * 4 + 1] |= bit;
                              args.PM[p * 4 + 2] |= bit; 
                              args.PM[p * 4 + 3] |= bit; break;
                }
            }
        }
        
        // Launch multi-pattern kernel
        dim3 block(256);
        dim3 grid((num_chunks + 255) / 256);
        myers_multi_pattern_kernel<<<grid, block>>>(args);
        check_cuda(cudaGetLastError(), "multi-pattern kernel launch");
        check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
        
        // Read back hit counts
        std::vector<uint32_t> hit_counts(num_patterns_in_batch);
        check_cuda(cudaMemcpy(hit_counts.data(), d_counters.ptr,
                              counters_size, cudaMemcpyDeviceToHost), "memcpy counters");
        
        // Read back sparse hits for each pattern
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            size_t result_idx = batch_start + p;
            uint32_t count = std::min(hit_counts[p], static_cast<uint32_t>(max_hits));
            
            if (count == 0) {
                // No hits for this pattern — leave distances vector empty
                continue;
            }
            
            // Read sparse hits for this pattern
            std::vector<SparseHit> sparse_hits(count);
            size_t read_offset = p * max_hits * sizeof(SparseHit);
            check_cuda(cudaMemcpy(sparse_hits.data(), 
                                  static_cast<char*>(d_hits.ptr) + read_offset,
                                  count * sizeof(SparseHit), 
                                  cudaMemcpyDeviceToHost), "memcpy sparse hits");
            
            // Expand sparse hits to full distances array
            // For backwards compatibility, we populate the full array
            // In the future, we could return sparse results directly
            results[result_idx].distances.resize(length, 255);  // 255 = no hit
            
            for (const auto& hit : sparse_hits) {
                if (hit.position < length) {
                    results[result_idx].distances[hit.position] = hit.distance;
                }
            }
        }
    }
    
    return results;
}

// TRUE SPARSE Multi-Pattern GPU Search (no 3GB array expansion)
// This is the high-performance path for batch off-target searches.
// Returns sparse results directly - caller doesn't need to allocate 3GB arrays.
// ───────────────────────────────────────────────────────────────────────────

std::vector<GpuSparseResult> myers_gpu_batch_sparse(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    
    if (patterns.empty()) {
        return {};
    }
    
    // ── Validate inputs ───────────────────────────────────────────────────
    if (length == 0)
        throw std::invalid_argument("myers_gpu_batch_sparse: length must be > 0");
    if (start + length > view.total_bases)
        throw std::invalid_argument("myers_gpu_batch_sparse: slice out of range");
    if (threshold < 0)
        throw std::invalid_argument("myers_gpu_batch_sparse: threshold must be >= 0");
    
    // Validate and normalize all patterns upfront
    std::vector<std::string> normalized_patterns;
    normalized_patterns.reserve(patterns.size());
    size_t max_pattern_len = 0;
    
    for (size_t i = 0; i < patterns.size(); ++i) {
        if (patterns[i].empty())
            throw std::invalid_argument("myers_gpu_batch_sparse: pattern " + 
                std::to_string(i) + " is empty");
        if (patterns[i].size() > 64)
            throw std::invalid_argument("myers_gpu_batch_sparse: pattern " + 
                std::to_string(i) + " exceeds 64 bp");
        
        normalized_patterns.push_back(validate_and_upper(patterns[i], "pattern"));
        max_pattern_len = std::max(max_pattern_len, patterns[i].size());
    }
    
    // ── Allocate and copy genome bit-vectors to GPU (once) ────────────────
    size_t first_word = start / 64;
    size_t last_word  = (start + length - 1) / 64 + 1;
    size_t num_words  = last_word - first_word;
    size_t bv_bytes   = num_words * sizeof(uint64_t);
    
    DevicePtr d_bvA, d_bvC, d_bvG, d_bvT;
    check_cuda(cudaMalloc(&d_bvA.ptr, bv_bytes), "cudaMalloc bv_A");
    check_cuda(cudaMalloc(&d_bvC.ptr, bv_bytes), "cudaMalloc bv_C");
    check_cuda(cudaMalloc(&d_bvG.ptr, bv_bytes), "cudaMalloc bv_G");
    check_cuda(cudaMalloc(&d_bvT.ptr, bv_bytes), "cudaMalloc bv_T");
    
    check_cuda(cudaMemcpy(d_bvA.ptr, view.bv_A + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_A");
    check_cuda(cudaMemcpy(d_bvC.ptr, view.bv_C + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_C");
    check_cuda(cudaMemcpy(d_bvG.ptr, view.bv_G + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_G");
    check_cuda(cudaMemcpy(d_bvT.ptr, view.bv_T + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_T");
    
    // ── Prepare sparse results ────────────────────────────────────────────
    std::vector<GpuSparseResult> results(patterns.size());
    for (size_t i = 0; i < patterns.size(); ++i) {
        results[i].pattern_len = normalized_patterns[i].size();
        results[i].text_len    = length;
    }
    
    // ── Process patterns in batches of MAX_PATTERNS_PER_KERNEL ────────────
    size_t stride = GPU_DEFAULT_STRIDE;
    size_t overlap = 2 * max_pattern_len;
    size_t num_chunks = (length + stride - 1) / stride;
    
    size_t max_hits = MAX_HITS_PER_PATTERN;
    
    for (size_t batch_start = 0; batch_start < patterns.size(); 
         batch_start += MAX_PATTERNS_PER_KERNEL) {
        
        size_t batch_end = std::min(batch_start + MAX_PATTERNS_PER_KERNEL, 
                                    patterns.size());
        int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);
        
        // Allocate sparse hit buffer and counters for this batch
        size_t hits_buffer_size = num_patterns_in_batch * max_hits * sizeof(SparseHit);
        size_t counters_size = num_patterns_in_batch * sizeof(uint32_t);
        
        DevicePtr d_hits, d_counters;
        check_cuda(cudaMalloc(&d_hits.ptr, hits_buffer_size), "cudaMalloc hits");
        check_cuda(cudaMalloc(&d_counters.ptr, counters_size), "cudaMalloc counters");
        
        // Zero the counters
        check_cuda(cudaMemset(d_counters.ptr, 0, counters_size), "memset counters");
        
        // Prepare multi-pattern kernel arguments
        MultiPatternArgs args;
        args.bv_A            = static_cast<const uint64_t*>(d_bvA.ptr);
        args.bv_C            = static_cast<const uint64_t*>(d_bvC.ptr);
        args.bv_G            = static_cast<const uint64_t*>(d_bvG.ptr);
        args.bv_T            = static_cast<const uint64_t*>(d_bvT.ptr);
        args.genome_start    = start;
        args.genome_len      = length;
        args.bv_word_offset  = first_word;
        args.num_patterns    = num_patterns_in_batch;
        args.threshold       = threshold;
        args.hits            = static_cast<SparseHit*>(d_hits.ptr);
        args.hit_counts      = static_cast<uint32_t*>(d_counters.ptr);
        args.max_hits_per_pattern = max_hits;
        args.stride          = stride;
        args.overlap         = overlap;
        
        // Build pattern masks for all patterns in batch
        std::memset(args.PM, 0, sizeof(args.PM));
        std::memset(args.pattern_lengths, 0, sizeof(args.pattern_lengths));
        
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            const std::string& pat = normalized_patterns[batch_start + p];
            args.pattern_lengths[p] = static_cast<int>(pat.size());
            
            for (size_t j = 0; j < pat.size(); ++j) {
                uint64_t bit = 1ULL << j;
                switch (pat[j]) {
                    case 'A': args.PM[p * 4 + 0] |= bit; break;
                    case 'C': args.PM[p * 4 + 1] |= bit; break;
                    case 'G': args.PM[p * 4 + 2] |= bit; break;
                    case 'T': args.PM[p * 4 + 3] |= bit; break;
                    case 'N': args.PM[p * 4 + 0] |= bit; 
                              args.PM[p * 4 + 1] |= bit;
                              args.PM[p * 4 + 2] |= bit; 
                              args.PM[p * 4 + 3] |= bit; break;
                }
            }
        }
        
        // Launch multi-pattern kernel
        dim3 block(256);
        dim3 grid((num_chunks + 255) / 256);
        myers_multi_pattern_kernel<<<grid, block>>>(args);
        check_cuda(cudaGetLastError(), "sparse multi-pattern kernel launch");
        check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
        
        // Read back hit counts
        std::vector<uint32_t> hit_counts(num_patterns_in_batch);
        check_cuda(cudaMemcpy(hit_counts.data(), d_counters.ptr,
                              counters_size, cudaMemcpyDeviceToHost), "memcpy counters");
        
        // Read back sparse hits and convert to GpuSparseHit format
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            size_t result_idx = batch_start + p;
            uint32_t count = std::min(hit_counts[p], static_cast<uint32_t>(max_hits));
            
            if (count == 0) {
                // No hits for this pattern
                continue;
            }
            
            // Read sparse hits for this pattern
            std::vector<SparseHit> sparse_hits(count);
            size_t read_offset = p * max_hits * sizeof(SparseHit);
            check_cuda(cudaMemcpy(sparse_hits.data(), 
                                  static_cast<char*>(d_hits.ptr) + read_offset,
                                  count * sizeof(SparseHit), 
                                  cudaMemcpyDeviceToHost), "memcpy sparse hits");
            
            // Convert to GpuSparseHit format (absolute positions)
            results[result_idx].hits.reserve(count);
            for (const auto& hit : sparse_hits) {
                GpuSparseHit gpu_hit;
                gpu_hit.position = start + hit.position;  // Convert to absolute position
                gpu_hit.distance = hit.distance;
                results[result_idx].hits.push_back(gpu_hit);
            }
        }
    }
    
    return results;
}


// ═══════════════════════════════════════════════════════════════════════════
// GPU Shift-Add Algorithm (Hamming Distance - Substitution-Only)
// ═══════════════════════════════════════════════════════════════════════════

// Shift-add kernel arguments (same structure as Myers, but simpler state)
struct ShiftAddKernelArgs {
    const uint64_t* bv_A;
    const uint64_t* bv_C;
    const uint64_t* bv_G;
    const uint64_t* bv_T;
    uint8_t*        distances;

    size_t genome_start;
    size_t genome_len;
    size_t bv_word_offset;

    uint64_t PM[4];
    int      m;

    size_t stride;
    size_t overlap;
};

// Shift-add multi-pattern kernel arguments
struct ShiftAddMultiPatternArgs {
    const uint64_t* bv_A;
    const uint64_t* bv_C;
    const uint64_t* bv_G;
    const uint64_t* bv_T;

    size_t genome_start;
    size_t genome_len;
    size_t bv_word_offset;

    uint64_t PM[MAX_PATTERNS_PER_KERNEL * 4];
    int      pattern_lengths[MAX_PATTERNS_PER_KERNEL];
    int      num_patterns;
    int      threshold;

    SparseHit* hits;
    uint32_t*  hit_counts;
    size_t     max_hits_per_pattern;

    size_t stride;
    size_t overlap;
};

// Shift-add single-pattern kernel
// Simpler than Myers: no Pv/Mv state, just accumulate Hamming distance.
// Algorithm:
//   D = word with all bits set to 1 (initially all positions mismatched)
//   At each text position:
//     D = (D << 1) | ~Eq   // Shift previous distances, set bit if mismatch
//     distance = popcount(D & mask)  // Count total mismatches
// ───────────────────────────────────────────────────────────────────────────

__global__ void shift_add_kernel(ShiftAddKernelArgs args) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    size_t num_chunks = (args.genome_len + args.stride - 1) / args.stride;
    if (tid >= num_chunks) return;

    // Output range (relative to genome_start)
    size_t out_start = tid * args.stride;
    size_t out_end   = out_start + args.stride;
    if (out_end > args.genome_len) out_end = args.genome_len;

    // Processing range includes warm-up overlap
    size_t proc_start = (tid == 0) ? 0 : (out_start - args.overlap);
    size_t proc_end   = out_end;

    int m = args.m;

    // State arrays: R[d] tracks positions matching with exactly d substitutions
    const int max_errors = m;
    uint64_t R[65]; // Support patterns up to 64bp
    for (int d = 0; d <= max_errors && d < 65; ++d) {
        R[d] = 0;
    }

    // Cached word index and loaded bit-vector words
    size_t   cached_word = ~(size_t)0;
    uint64_t wA = 0, wC = 0, wG = 0, wT = 0;

    // Process each genome position
    for (size_t rel = proc_start; rel < proc_end; ++rel) {
        size_t abs_pos = args.genome_start + rel;
        size_t  word   = abs_pos / 64;
        uint64_t bit   = 1ULL << (abs_pos % 64);

        // Load bit-vector words if needed
        if (word != cached_word) {
            size_t dev_idx = word - args.bv_word_offset;
            wA = args.bv_A[dev_idx];
            wC = args.bv_C[dev_idx];
            wG = args.bv_G[dev_idx];
            wT = args.bv_T[dev_idx];
            cached_word = word;
        }

        // Decode nucleotide and get character mask
        bool bA = wA & bit;
        bool bC = wC & bit;
        bool bG = wG & bit;
        bool bT = wT & bit;

        uint64_t c_mask;
        // Masked N in genome: doesn't match anything (c_mask = 0)
        // N in pattern (wildcards) are already encoded in PM arrays
        if (!bA && !bC && !bG && !bT) {
            c_mask = 0;  // No match
        } else if (bA) {
            c_mask = args.PM[0];
        } else if (bC) {
            c_mask = args.PM[1];
        } else if (bG) {
            c_mask = args.PM[2];
        } else {
            c_mask = args.PM[3];  // T
        }

        // Bit-parallel shift-and update
        uint64_t old_R = R[0];
        
        // Update R[0]: exact matches
        R[0] = ((old_R << 1) | 1) & c_mask;

        // Update R[1..max_errors]: matches with substitutions
        for (int d = 1; d <= max_errors && d < 65; ++d) {
            uint64_t temp = R[d];
            R[d] = ((temp << 1) & c_mask) | (old_R << 1) | 1;
            old_R = temp;
        }

        // Find minimum distance at this position
        // For semi-global alignment: check bit positions [0..min(rel, m-1)] 
        // This allows partial matches at the beginning
        size_t check_bit = (rel < (size_t)(m - 1)) ? rel : (m - 1);
        uint64_t check_mask = 1ULL << check_bit;
        
        uint8_t distance = static_cast<uint8_t>(m); // Default: no match
        for (int d = 0; d <= max_errors && d < 65; ++d) {
            if (R[d] & check_mask) {
                distance = static_cast<uint8_t>(d);
                break;
            }
        }

        // Write output only for positions inside this thread's output range
        if (rel >= out_start)
            args.distances[rel] = distance;
    }
}


// Shift-add multi-pattern kernel (sparse output)

template <int MAX_D>
__global__ void shift_add_multi_pattern_kernel_t(ShiftAddMultiPatternArgs args) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    size_t num_chunks = (args.genome_len + args.stride - 1) / args.stride;
    if (tid >= num_chunks) return;

    size_t out_start = tid * args.stride;
    size_t out_end   = out_start + args.stride;
    if (out_end > args.genome_len) out_end = args.genome_len;

    size_t proc_start = (tid == 0) ? 0 : (out_start - args.overlap);
    size_t proc_end   = out_end;

    // Compile-time-sized state: threshold+1 slots instead of fixed 12.
    // At threshold=3 this is 4 slots × 32 patterns = 1 KB per thread (was 3 KB).
    uint64_t R[MAX_PATTERNS_PER_KERNEL][MAX_D];

    for (int p = 0; p < args.num_patterns; ++p) {
        #pragma unroll
        for (int d = 0; d < MAX_D; ++d) R[p][d] = 0;
    }

    size_t   cached_word = ~(size_t)0;
    uint64_t wA = 0, wC = 0, wG = 0, wT = 0;

    for (size_t rel = proc_start; rel < proc_end; ++rel) {
        size_t abs_pos = args.genome_start + rel;
        size_t  word   = abs_pos / 64;
        uint64_t bit   = 1ULL << (abs_pos % 64);

        if (word != cached_word) {
            size_t dev_idx = word - args.bv_word_offset;
            wA = args.bv_A[dev_idx];
            wC = args.bv_C[dev_idx];
            wG = args.bv_G[dev_idx];
            wT = args.bv_T[dev_idx];
            cached_word = word;
        }

        bool bA = wA & bit;
        bool bC = wC & bit;
        bool bG = wG & bit;
        bool bT = wT & bit;

        // Masked N: nuc_idx is unused (c_mask forced to 0 below).
        bool   valid   = bA | bC | bG | bT;
        int    nuc_idx = bA ? 0 : (bC ? 1 : (bG ? 2 : 3));

        for (int p = 0; p < args.num_patterns; ++p) {
            uint64_t c_mask = valid ? args.PM[p * 4 + nuc_idx] : 0ULL;

            // Update R[0..MAX_D-1] with unified recurrence.
            // OR-reduce as we go so the hot-path "is there a hit?" test is O(1).
            uint64_t old_R    = R[p][0];
            uint64_t new_R0   = ((old_R << 1) | 1) & c_mask;
            R[p][0]           = new_R0;
            uint64_t hit_mask = new_R0;

            #pragma unroll
            for (int d = 1; d < MAX_D; ++d) {
                uint64_t temp = R[p][d];
                uint64_t new_Rd = ((temp << 1) & c_mask) | (old_R << 1) | 1;
                R[p][d]   = new_Rd;
                hit_mask |= new_Rd;
                old_R     = temp;
            }

            // Skip emission at N positions (match original kernel & scalar kernel
            // semantics). State is already updated above with c_mask=0.
            if (!valid || rel < out_start) continue;

            int m = args.pattern_lengths[p];
            size_t cbit = (rel < (size_t)(m - 1)) ? rel : (size_t)(m - 1);
            uint64_t check_mask = 1ULL << cbit;

            // Fast early-reject: no R[] slot up to MAX_D-1 holds this position.
            if (!(hit_mask & check_mask)) continue;

            // Find minimum distance only now — rare path (~10^-4 per position).
            int distance = 0;
            #pragma unroll
            for (int d = 0; d < MAX_D; ++d) {
                if (R[p][d] & check_mask) { distance = d; break; }
            }

            uint32_t idx = atomicAdd(&args.hit_counts[p], 1u);
            if (idx < args.max_hits_per_pattern) {
                size_t write_pos = p * args.max_hits_per_pattern + idx;
                args.hits[write_pos].position    = static_cast<uint32_t>(rel);
                args.hits[write_pos].distance    = static_cast<uint8_t>(distance);
                args.hits[write_pos].pattern_idx = static_cast<uint8_t>(p);
            }
        }
    }
}

// Dispatch helper: picks the right template instantiation based on threshold.
// MAX_D = threshold + 1; supported range mirrors the original 12-slot bound.
static inline void launch_shift_add_multi(int threshold,
                                          dim3 grid, dim3 block,
                                          cudaStream_t stream,
                                          const ShiftAddMultiPatternArgs& args) {
    switch (threshold + 1) {
        case 1:  shift_add_multi_pattern_kernel_t<1> <<<grid, block, 0, stream>>>(args); break;
        case 2:  shift_add_multi_pattern_kernel_t<2> <<<grid, block, 0, stream>>>(args); break;
        case 3:  shift_add_multi_pattern_kernel_t<3> <<<grid, block, 0, stream>>>(args); break;
        case 4:  shift_add_multi_pattern_kernel_t<4> <<<grid, block, 0, stream>>>(args); break;
        case 5:  shift_add_multi_pattern_kernel_t<5> <<<grid, block, 0, stream>>>(args); break;
        case 6:  shift_add_multi_pattern_kernel_t<6> <<<grid, block, 0, stream>>>(args); break;
        case 7:  shift_add_multi_pattern_kernel_t<7> <<<grid, block, 0, stream>>>(args); break;
        case 8:  shift_add_multi_pattern_kernel_t<8> <<<grid, block, 0, stream>>>(args); break;
        case 9:  shift_add_multi_pattern_kernel_t<9> <<<grid, block, 0, stream>>>(args); break;
        case 10: shift_add_multi_pattern_kernel_t<10><<<grid, block, 0, stream>>>(args); break;
        case 11: shift_add_multi_pattern_kernel_t<11><<<grid, block, 0, stream>>>(args); break;
        case 12: shift_add_multi_pattern_kernel_t<12><<<grid, block, 0, stream>>>(args); break;
        default:
            // Should be unreachable; host-side validation caps threshold at 11.
            shift_add_multi_pattern_kernel_t<12><<<grid, block, 0, stream>>>(args);
            break;
    }
}

// Public API: Single-pattern shift-add GPU

GpuShiftAddResult shift_add_gpu(const std::string& pattern,
                                GenomeView         view,
                                size_t             start,
                                size_t             length) {
    // ── Input validation ──────────────────────────────────────────────────
    if (pattern.empty())
        throw std::invalid_argument("shift_add_gpu: pattern must not be empty");
    if (pattern.size() > 64)
        throw std::invalid_argument("shift_add_gpu: pattern length must be ≤ 64");
    if (length == 0)
        throw std::invalid_argument("shift_add_gpu: length must be > 0");
    if (start + length > view.total_bases)
        throw std::invalid_argument("shift_add_gpu: slice out of range");

    std::string pat = validate_and_upper(pattern, "pattern");
    size_t m = pat.size();

    // ── Build pattern-match masks ─────────────────────────────────────────
    uint64_t PM[4] = {0, 0, 0, 0};
    for (size_t i = 0; i < m; ++i) {
        uint64_t bit = 1ULL << i;
        switch (pat[i]) {
            case 'A': PM[0] |= bit; break;
            case 'C': PM[1] |= bit; break;
            case 'G': PM[2] |= bit; break;
            case 'T': PM[3] |= bit; break;
            case 'N': PM[0] |= bit; PM[1] |= bit;
                      PM[2] |= bit; PM[3] |= bit; break;
        }
    }

    // ── Determine which bit-vector words are needed ───────────────────────
    size_t first_word = start / 64;
    size_t last_word  = (start + length - 1) / 64 + 1;
    size_t num_words  = last_word - first_word;
    size_t bv_bytes   = num_words * sizeof(uint64_t);

    // ── Allocate device memory ────────────────────────────────────────────
    DevicePtr d_bvA, d_bvC, d_bvG, d_bvT, d_dist;
    check_cuda(cudaMalloc(&d_bvA.ptr, bv_bytes), "cudaMalloc bv_A");
    check_cuda(cudaMalloc(&d_bvC.ptr, bv_bytes), "cudaMalloc bv_C");
    check_cuda(cudaMalloc(&d_bvG.ptr, bv_bytes), "cudaMalloc bv_G");
    check_cuda(cudaMalloc(&d_bvT.ptr, bv_bytes), "cudaMalloc bv_T");
    check_cuda(cudaMalloc(&d_dist.ptr, length), "cudaMalloc distances");

    // ── Copy genome bit-vectors to device ─────────────────────────────────
    check_cuda(cudaMemcpy(d_bvA.ptr, view.bv_A + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_A");
    check_cuda(cudaMemcpy(d_bvC.ptr, view.bv_C + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_C");
    check_cuda(cudaMemcpy(d_bvG.ptr, view.bv_G + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_G");
    check_cuda(cudaMemcpy(d_bvT.ptr, view.bv_T + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_T");

    // ── Prepare kernel arguments ──────────────────────────────────────────
    size_t stride  = GPU_DEFAULT_STRIDE;
    size_t overlap = 2 * m;
    size_t num_chunks = (length + stride - 1) / stride;

    ShiftAddKernelArgs args;
    args.bv_A            = static_cast<const uint64_t*>(d_bvA.ptr);
    args.bv_C            = static_cast<const uint64_t*>(d_bvC.ptr);
    args.bv_G            = static_cast<const uint64_t*>(d_bvG.ptr);
    args.bv_T            = static_cast<const uint64_t*>(d_bvT.ptr);
    args.distances       = static_cast<uint8_t*>(d_dist.ptr);
    args.genome_start    = start;
    args.genome_len      = length;
    args.bv_word_offset  = first_word;
    args.PM[0]           = PM[0];
    args.PM[1]           = PM[1];
    args.PM[2]           = PM[2];
    args.PM[3]           = PM[3];
    args.m               = static_cast<int>(m);
    args.stride          = stride;
    args.overlap         = overlap;

    // ── Launch kernel ─────────────────────────────────────────────────────
    dim3 block(256);
    dim3 grid((num_chunks + 255) / 256);
    shift_add_kernel<<<grid, block>>>(args);
    check_cuda(cudaGetLastError(), "shift_add_kernel launch");
    check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

    // ── Copy results back ─────────────────────────────────────────────────
    std::vector<uint8_t> distances(length);
    check_cuda(cudaMemcpy(distances.data(), d_dist.ptr,
                          length, cudaMemcpyDeviceToHost), "memcpy distances");

    return {distances, m, length};
}

// Public API: Multi-pattern shift-add with threshold (full array output)

std::vector<GpuShiftAddResult> shift_add_gpu_batch_threshold(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    
    // Reuse sparse implementation and expand to full arrays
    auto sparse_results = shift_add_gpu_batch_sparse(patterns, view, start, length, threshold);
    
    // Convert sparse to full distance arrays
    std::vector<GpuShiftAddResult> results(patterns.size());
    for (size_t i = 0; i < patterns.size(); ++i) {
        results[i].pattern_len = sparse_results[i].pattern_len;
        results[i].text_len    = sparse_results[i].text_len;
        results[i].distances.resize(length, 255);  // 255 = no hit
        
        for (const auto& hit : sparse_results[i].hits) {
            size_t rel_pos = hit.position - start;
            if (rel_pos < length) {
                results[i].distances[rel_pos] = hit.distance;
            }
        }
    }
    
    return results;
}

// Public API: Multi-pattern shift-add sparse output (preferred)

std::vector<GpuSparseResult> shift_add_gpu_batch_sparse(
    const std::vector<std::string>& patterns,
    GenomeView                      view,
    size_t                          start,
    size_t                          length,
    int                             threshold) {
    
    if (patterns.empty()) {
        return {};
    }
    
    // ── Validate inputs ───────────────────────────────────────────────────
    if (length == 0)
        throw std::invalid_argument("shift_add_gpu_batch_sparse: length must be > 0");
    if (start + length > view.total_bases)
        throw std::invalid_argument("shift_add_gpu_batch_sparse: slice out of range");
    if (threshold < 0)
        throw std::invalid_argument("shift_add_gpu_batch_sparse: threshold must be >= 0");
    
    // Validate and normalize all patterns
    std::vector<std::string> normalized_patterns;
    normalized_patterns.reserve(patterns.size());
    size_t max_pattern_len = 0;
    
    for (size_t i = 0; i < patterns.size(); ++i) {
        if (patterns[i].empty())
            throw std::invalid_argument("shift_add_gpu_batch_sparse: pattern " + 
                std::to_string(i) + " is empty");
        if (patterns[i].size() > 64)
            throw std::invalid_argument("shift_add_gpu_batch_sparse: pattern " + 
                std::to_string(i) + " exceeds 64 bp");
        
        normalized_patterns.push_back(validate_and_upper(patterns[i], "pattern"));
        max_pattern_len = std::max(max_pattern_len, patterns[i].size());
    }
    
    // ── Allocate and copy genome bit-vectors to GPU (once) ────────────────
    size_t first_word = start / 64;
    size_t last_word  = (start + length - 1) / 64 + 1;
    size_t num_words  = last_word - first_word;
    size_t bv_bytes   = num_words * sizeof(uint64_t);
    
    DevicePtr d_bvA, d_bvC, d_bvG, d_bvT;
    check_cuda(cudaMalloc(&d_bvA.ptr, bv_bytes), "cudaMalloc bv_A");
    check_cuda(cudaMalloc(&d_bvC.ptr, bv_bytes), "cudaMalloc bv_C");
    check_cuda(cudaMalloc(&d_bvG.ptr, bv_bytes), "cudaMalloc bv_G");
    check_cuda(cudaMalloc(&d_bvT.ptr, bv_bytes), "cudaMalloc bv_T");
    
    check_cuda(cudaMemcpy(d_bvA.ptr, view.bv_A + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_A");
    check_cuda(cudaMemcpy(d_bvC.ptr, view.bv_C + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_C");
    check_cuda(cudaMemcpy(d_bvG.ptr, view.bv_G + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_G");
    check_cuda(cudaMemcpy(d_bvT.ptr, view.bv_T + first_word,
                          bv_bytes, cudaMemcpyHostToDevice), "memcpy bv_T");
    
    // ── Prepare sparse results ────────────────────────────────────────────
    std::vector<GpuSparseResult> results(patterns.size());
    for (size_t i = 0; i < patterns.size(); ++i) {
        results[i].pattern_len = normalized_patterns[i].size();
        results[i].text_len    = length;
    }
    
    // ── Process patterns in batches ───────────────────────────────────────
    size_t stride = GPU_DEFAULT_STRIDE;
    size_t overlap = 2 * max_pattern_len;
    size_t num_chunks = (length + stride - 1) / stride;
    
    size_t max_hits = MAX_HITS_PER_PATTERN;
    
    for (size_t batch_start = 0; batch_start < patterns.size(); 
         batch_start += MAX_PATTERNS_PER_KERNEL) {
        
        size_t batch_end = std::min(batch_start + MAX_PATTERNS_PER_KERNEL, 
                                    patterns.size());
        int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);
        
        // Allocate sparse hit buffer and counters
        size_t hits_buffer_size = num_patterns_in_batch * max_hits * sizeof(SparseHit);
        size_t counters_size = num_patterns_in_batch * sizeof(uint32_t);
        
        DevicePtr d_hits, d_counters;
        check_cuda(cudaMalloc(&d_hits.ptr, hits_buffer_size), "cudaMalloc hits");
        check_cuda(cudaMalloc(&d_counters.ptr, counters_size), "cudaMalloc counters");
        
        check_cuda(cudaMemset(d_counters.ptr, 0, counters_size), "memset counters");
        
        // Prepare kernel arguments
        ShiftAddMultiPatternArgs args;
        args.bv_A            = static_cast<const uint64_t*>(d_bvA.ptr);
        args.bv_C            = static_cast<const uint64_t*>(d_bvC.ptr);
        args.bv_G            = static_cast<const uint64_t*>(d_bvG.ptr);
        args.bv_T            = static_cast<const uint64_t*>(d_bvT.ptr);
        args.genome_start    = start;
        args.genome_len      = length;
        args.bv_word_offset  = first_word;
        args.num_patterns    = num_patterns_in_batch;
        args.threshold       = threshold;
        args.hits            = static_cast<SparseHit*>(d_hits.ptr);
        args.hit_counts      = static_cast<uint32_t*>(d_counters.ptr);
        args.max_hits_per_pattern = max_hits;
        args.stride          = stride;
        args.overlap         = overlap;
        
        // Build pattern masks
        std::memset(args.PM, 0, sizeof(args.PM));
        std::memset(args.pattern_lengths, 0, sizeof(args.pattern_lengths));
        
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            const std::string& pat = normalized_patterns[batch_start + p];
            args.pattern_lengths[p] = static_cast<int>(pat.size());
            
            for (size_t j = 0; j < pat.size(); ++j) {
                uint64_t bit = 1ULL << j;
                switch (pat[j]) {
                    case 'A': args.PM[p * 4 + 0] |= bit; break;
                    case 'C': args.PM[p * 4 + 1] |= bit; break;
                    case 'G': args.PM[p * 4 + 2] |= bit; break;
                    case 'T': args.PM[p * 4 + 3] |= bit; break;
                    case 'N': args.PM[p * 4 + 0] |= bit; 
                              args.PM[p * 4 + 1] |= bit;
                              args.PM[p * 4 + 2] |= bit; 
                              args.PM[p * 4 + 3] |= bit; break;
                }
            }
        }
        
        // Launch kernel
        dim3 block(256);
        dim3 grid((num_chunks + 255) / 256);
        launch_shift_add_multi(threshold, grid, block, /*stream=*/0, args);
        check_cuda(cudaGetLastError(), "shift_add multi-pattern kernel launch");
        check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
        
        // Read back hit counts
        std::vector<uint32_t> hit_counts(num_patterns_in_batch);
        check_cuda(cudaMemcpy(hit_counts.data(), d_counters.ptr,
                              counters_size, cudaMemcpyDeviceToHost), "memcpy counters");
        
        // Read back sparse hits
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            size_t result_idx = batch_start + p;
            uint32_t count = std::min(hit_counts[p], static_cast<uint32_t>(max_hits));
            
            if (count == 0) continue;
            
            std::vector<SparseHit> sparse_hits(count);
            size_t read_offset = p * max_hits * sizeof(SparseHit);
            check_cuda(cudaMemcpy(sparse_hits.data(), 
                                  static_cast<char*>(d_hits.ptr) + read_offset,
                                  count * sizeof(SparseHit), 
                                  cudaMemcpyDeviceToHost), "memcpy sparse hits");
            
            // Convert to GpuSparseHit format
            results[result_idx].hits.reserve(count);
            for (const auto& hit : sparse_hits) {
                GpuSparseHit gpu_hit;
                gpu_hit.position = start + hit.position;
                gpu_hit.distance = hit.distance;
                results[result_idx].hits.push_back(gpu_hit);
            }
        }
    }

    return results;
}

