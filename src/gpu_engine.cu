#include <gpu_engine.hpp>
#include <bio_utils.hpp>

#include <cuda_runtime.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <atomic>

// Maximum patterns per multi-pattern kernel launch.
// Sized to fit pattern masks in kernel args: 64 × 4 × 8 = 2 KB (well under
// the 4 KB argument-struct cap). Per-thread state still scales linearly with
// actual batch size chosen via STOMATA_SHIFTADD_BATCH / STOMATA_MYERS_BATCH (see
// dispatchers below) — the default batch is 32; 64 is available for library-
// scale (1k–10k spacer) runs where halving the batch count matters more than
// the doubled per-thread L2 footprint.
static constexpr int MAX_PATTERNS_PER_KERNEL = 64;

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

// Threads per block. Overridable via STOMATA_BLOCK=<n> (must be a multiple of
// 32 and a power of two between 32 and 1024). Defaults to 256, which matches
// the legacy launch and keeps `(num_chunks + 255) / 256` arithmetic valid.
static int choose_gpu_block() {
    if (const char* ov = std::getenv("STOMATA_BLOCK")) {
        int v = std::atoi(ov);
        if (v == 64 || v == 128 || v == 256 || v == 512 || v == 1024) return v;
    }
    return 256;
}

// Pick a stride (genome positions per thread) that scales the launch to the
// current device. The multi-pattern kernel spills ~1 KB/thread of local memory
// (R[32][threshold+1] uint64_t state), so total in-flight local memory is
// bounded by L2 cache size — pushing past it thrashes the cache and slows the
// kernel by 2–3×. An RTX 6000 Ada has 48 MB L2; budgeting ~60% for per-thread
// state at 1 KB/thread gives ~29K threads, i.e. stride ≈ genome/29K. On a
// human genome that's ~105K, about 2× the legacy GPU_DEFAULT_STRIDE.
//
// Overridable via STOMATA_STRIDE=<bases> for benchmarking. Falls back to
// GPU_DEFAULT_STRIDE if the device query fails.
// `kind` hints at the kernel's bottleneck character. Shift-add (Hamming) is
// memory-bound and sensitive to L2 footprint: ~1 block/SM is the sweet spot.
// Myers (Levenshtein) is compute-bound (~10 ops per position per pattern,
// including an ADD for carry propagation): it wants more threads to fill the
// SM scheduler and tolerates some L2 spill. The two strides can differ by
// ~4× on the same hardware; see `benchmark_results/myers_stride_20260422`
// and `opt_uint32_stride_20260422`.
enum class GpuKernelKind { Hamming, Myers };

static size_t choose_gpu_stride(size_t length, size_t max_pattern_len,
                                int threshold = 3, size_t word_bytes = 8,
                                GpuKernelKind kind = GpuKernelKind::Hamming,
                                int actual_max_p = MAX_PATTERNS_PER_KERNEL) {
    if (const char* ov = std::getenv("STOMATA_STRIDE")) {
        long v = std::atol(ov);
        if (v > 0) return static_cast<size_t>(v);
    }
    int device = 0;
    if (cudaGetDevice(&device) != cudaSuccess) return GPU_DEFAULT_STRIDE;
    int l2_bytes = 0;
    int sm_count = 0;
    cudaDeviceGetAttribute(&l2_bytes, cudaDevAttrL2CacheSize,        device);
    cudaDeviceGetAttribute(&sm_count, cudaDevAttrMultiProcessorCount, device);
    if (l2_bytes <= 0 || sm_count <= 0) return GPU_DEFAULT_STRIDE;

    // Per-thread stack frame scales with actual batch size (not the compile-
    // time MAX_PATTERNS_PER_KERNEL cap). Shift-add state is R[MAX_P][MAX_D] ×
    // word_bytes; Myers state is (Pv+Mv+mask+top_bit)[MAX_P] × word_bytes +
    // score[MAX_P] × 4.
    const int max_p = std::max(1, actual_max_p);
    size_t stack_per_thread;
    if (kind == GpuKernelKind::Myers) {
        stack_per_thread = static_cast<size_t>(max_p) * (4 * word_bytes + 4);
    } else {
        stack_per_thread = word_bytes *
                           static_cast<size_t>(max_p) *
                           static_cast<size_t>(std::max(1, threshold + 1));
    }

    size_t target_threads;
    if (kind == GpuKernelKind::Myers) {
        // Myers is compute-bound; target ~2 blocks/SM to land in the broad
        // sweet spot (stride 32K–49K on RTX 6000 Ada) without risking the
        // cliff at very small strides (~16K = L2 thrash).
        target_threads = static_cast<size_t>(sm_count) * 2 * 256;
    } else {
        // Shift-add is memory-bound; ~1 block/SM with L2 as a safety cap.
        size_t sm_target = static_cast<size_t>(sm_count) * 256;
        size_t l2_safe = (static_cast<size_t>(l2_bytes) * 3 / 4) /
                         std::max<size_t>(1, stack_per_thread);
        target_threads = std::min(sm_target, l2_safe);
        if (target_threads == 0) target_threads = sm_target;
    }

    size_t stride = (length + target_threads - 1) / target_threads;
    size_t min_stride = std::max<size_t>(4096, 100 * std::max<size_t>(1, max_pattern_len));
    if (stride < min_stride) stride = min_stride;
    return stride;
}

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

    // Sparse output — coalesced: all patterns share one flat buffer and one
    // atomic counter, with pattern_idx embedded in each SparseHit. Host uses
    // a single D2H transfer per batch.
    SparseHit* hits;
    uint32_t*  global_hit_count;
    size_t     max_total_hits;

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

// Multi-pattern Myers kernel: up to MAX_P patterns per genome pass.
// Templated on word type W (uint32_t for patterns ≤ 32 bp, uint64_t for
// up to 64 bp) and on MAX_P (compile-time pattern-slot count). Smaller
// MAX_P shrinks per-thread stack linearly — at MAX_P=16, uint32_t the
// stack is 320 B vs 640 B at MAX_P=32. Caller chooses MAX_P based on
// the batch size it wants per launch.
template <typename W, int MAX_P>
__global__ void myers_multi_pattern_kernel_t(MultiPatternArgs args) {
    constexpr int W_BITS = sizeof(W) * 8;
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
    W   Pv[MAX_P];
    W   Mv[MAX_P];
    int score[MAX_P];
    W   mask[MAX_P];
    W   top_bit[MAX_P];

    // Initialize all patterns
    for (int p = 0; p < args.num_patterns; ++p) {
        int m = args.pattern_lengths[p];
        mask[p]    = (m == W_BITS) ? ~W{0} : ((W{1} << m) - W{1});
        top_bit[p] = W{1} << (m - 1);
        Pv[p]      = mask[p];
        Mv[p]      = W{0};
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
            // Get Eq mask for this pattern and nucleotide. PM is stored as
            // uint64_t on the host; narrow to W (low bits are all that matter
            // for patterns ≤ W_BITS bp).
            W Eq = (nuc_idx >= 0) ? static_cast<W>(args.PM[p * 4 + nuc_idx]) : W{0};

            // Myers bit-parallel update (width-agnostic — the carry in Xh
            // is contained within W because Pv and Eq both have bits only in
            // positions < pattern_len ≤ W_BITS, and mask strips bit W_BITS at
            // the end of the iteration).
            W Xv = Eq | Mv[p];
            W Xh = (((Eq & Pv[p]) + Pv[p]) ^ Pv[p]) | Eq | Mv[p];
            W Ph = Mv[p] | ~(Xh | Pv[p]);
            W Mh = Pv[p] & Xh;

            // Score update
            if (Ph & top_bit[p]) ++score[p];
            if (Mh & top_bit[p]) --score[p];

            // Shift horizontal deltas
            Ph = (Ph << 1) & mask[p];
            Mh = (Mh << 1) & mask[p];

            // Update vertical deltas
            Pv[p] = (Mh | ~(Xv | Ph)) & mask[p];
            Mv[p] = Ph & Xv;

            // Output sparse hit if within output range and passes threshold.
            // Coalesced write: single batch-wide atomic, pattern_idx embedded.
            if (rel >= out_start && score[p] <= args.threshold) {
                uint32_t idx = atomicAdd(args.global_hit_count, 1);
                if (idx < args.max_total_hits) {
                    args.hits[idx].position    = static_cast<uint32_t>(rel);
                    args.hits[idx].distance    = static_cast<uint8_t>(score[p]);
                    args.hits[idx].pattern_idx = static_cast<uint8_t>(p);
                }
            }
        }
    }
}

// Myers per-kernel-launch pattern batch size. Default 16.
// `STOMATA_MYERS_BATCH=32` — recover the previous behavior.
// `STOMATA_MYERS_BATCH=64` — double-batch for library-scale runs (halves the
//                          kernel launch count at the cost of 4× per-thread
//                          Myers state — compute-bound path, budget carefully)
// Smaller value shrinks per-thread stack (Pv/Mv/mask/top_bit arrays) linearly,
// improving L2 fit at the cost of more kernel launches. A half-batch run on
// 225 Lev+NGG spacers is ~5% faster than 32 (27.9 s vs 29.2 s on RTX 6000 Ada;
// see benchmark_results/myers_halfbatch_20260423). The extra launches are
// amortized — for the 225-spacer workload, 15 launches instead of 8,
// ~milliseconds of overhead.
static int myers_batch_size() {
    if (const char* ov = std::getenv("STOMATA_MYERS_BATCH")) {
        int v = std::atoi(ov);
        if (v == 16 || v == 32 || v == 64) return v;
    }
    return 16;
}

// Dispatch helper for the Myers multi-pattern kernel.
//   `use_32bit`   — uint32_t state variant when every pattern ≤ 32 bp
//   `max_p`       — compile-time pattern-slot count (16 or 32); caller must
//                   batch no more than `max_p` patterns per launch. Smaller
//                   MAX_P halves per-thread stack and improves L2 fit.
static inline void launch_myers_multi(bool use_32bit, int max_p,
                                      dim3 grid, dim3 block,
                                      const MultiPatternArgs& args) {
    if (use_32bit) {
        if      (max_p <= 16) myers_multi_pattern_kernel_t<uint32_t, 16><<<grid, block>>>(args);
        else if (max_p <= 32) myers_multi_pattern_kernel_t<uint32_t, 32><<<grid, block>>>(args);
        else                  myers_multi_pattern_kernel_t<uint32_t, 64><<<grid, block>>>(args);
    } else {
        if      (max_p <= 16) myers_multi_pattern_kernel_t<uint64_t, 16><<<grid, block>>>(args);
        else if (max_p <= 32) myers_multi_pattern_kernel_t<uint64_t, 32><<<grid, block>>>(args);
        else                  myers_multi_pattern_kernel_t<uint64_t, 64><<<grid, block>>>(args);
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
    size_t stride     = choose_gpu_stride(length, static_cast<size_t>(m),
                                           3, 8, GpuKernelKind::Myers, 1);
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
    
    // ── Process patterns in batches of `batch_size` (≤ MAX_PATTERNS_PER_KERNEL)
    // Use the uint32_t Myers variant when every pattern fits in 32 bits.
    const bool use_32bit = (max_pattern_len <= 32);
    const size_t word_bytes = use_32bit ? 4 : 8;
    const int batch_size = myers_batch_size();

    size_t stride = choose_gpu_stride(length, max_pattern_len, threshold,
                                       word_bytes, GpuKernelKind::Myers,
                                       batch_size);
    size_t overlap = 2 * max_pattern_len;
    size_t num_chunks = (length + stride - 1) / stride;

    // Allocate sparse output buffers
    size_t max_hits = MAX_HITS_PER_PATTERN;

    for (size_t batch_start = 0; batch_start < patterns.size();
         batch_start += batch_size) {

        size_t batch_end = std::min(batch_start + batch_size,
                                    patterns.size());
        int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);

        // Coalesced output: one flat hit buffer + one global counter per batch.
        size_t max_total_hits     = static_cast<size_t>(num_patterns_in_batch) * max_hits;
        size_t hits_buffer_size   = max_total_hits * sizeof(SparseHit);

        DevicePtr d_hits, d_counter;
        check_cuda(cudaMalloc(&d_hits.ptr,    hits_buffer_size), "cudaMalloc hits");
        check_cuda(cudaMalloc(&d_counter.ptr, sizeof(uint32_t)), "cudaMalloc global_hit_count");
        check_cuda(cudaMemset(d_counter.ptr, 0, sizeof(uint32_t)), "memset global_hit_count");

        MultiPatternArgs args;
        args.bv_A              = static_cast<const uint64_t*>(d_bvA.ptr);
        args.bv_C              = static_cast<const uint64_t*>(d_bvC.ptr);
        args.bv_G              = static_cast<const uint64_t*>(d_bvG.ptr);
        args.bv_T              = static_cast<const uint64_t*>(d_bvT.ptr);
        args.genome_start      = start;
        args.genome_len        = length;
        args.bv_word_offset    = first_word;
        args.num_patterns      = num_patterns_in_batch;
        args.threshold         = threshold;
        args.hits              = static_cast<SparseHit*>(d_hits.ptr);
        args.global_hit_count  = static_cast<uint32_t*>(d_counter.ptr);
        args.max_total_hits    = max_total_hits;
        args.stride            = stride;
        args.overlap           = overlap;

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

        int block_sz = choose_gpu_block();
        dim3 block(block_sz);
        dim3 grid((num_chunks + block_sz - 1) / block_sz);
        launch_myers_multi(use_32bit, batch_size, grid, block, args);
        check_cuda(cudaGetLastError(), "multi-pattern kernel launch");
        check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

        // Single D2H for total count, then single D2H for packed hits.
        uint32_t total_hits = 0;
        check_cuda(cudaMemcpy(&total_hits, d_counter.ptr,
                              sizeof(uint32_t), cudaMemcpyDeviceToHost),
                   "memcpy global_hit_count");
        if (total_hits > max_total_hits) total_hits = static_cast<uint32_t>(max_total_hits);
        if (total_hits == 0) continue;

        std::vector<SparseHit> all_hits(total_hits);
        check_cuda(cudaMemcpy(all_hits.data(), d_hits.ptr,
                              static_cast<size_t>(total_hits) * sizeof(SparseHit),
                              cudaMemcpyDeviceToHost),
                   "memcpy packed hits");

        // Lazily expand per-pattern distances from the packed hits.
        for (uint32_t i = 0; i < total_hits; ++i) {
            const SparseHit& hit = all_hits[i];
            if (hit.pattern_idx >= num_patterns_in_batch) continue;
            size_t result_idx = batch_start + hit.pattern_idx;
            if (results[result_idx].distances.empty()) {
                results[result_idx].distances.assign(length, 255);
            }
            if (hit.position < length) {
                results[result_idx].distances[hit.position] = hit.distance;
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
    
    // ── Process patterns in batches of `batch_size` ──────────────────────
    const bool use_32bit = (max_pattern_len <= 32);
    const size_t word_bytes = use_32bit ? 4 : 8;
    const int batch_size = myers_batch_size();

    size_t stride = choose_gpu_stride(length, max_pattern_len, threshold,
                                       word_bytes, GpuKernelKind::Myers,
                                       batch_size);
    size_t overlap = 2 * max_pattern_len;
    size_t num_chunks = (length + stride - 1) / stride;

    size_t max_hits = MAX_HITS_PER_PATTERN;

    for (size_t batch_start = 0; batch_start < patterns.size();
         batch_start += batch_size) {

        size_t batch_end = std::min(batch_start + batch_size,
                                    patterns.size());
        int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);

        // Coalesced output (see shift-add path for rationale).
        size_t max_total_hits   = static_cast<size_t>(num_patterns_in_batch) * max_hits;
        size_t hits_buffer_size = max_total_hits * sizeof(SparseHit);

        DevicePtr d_hits, d_counter;
        check_cuda(cudaMalloc(&d_hits.ptr,    hits_buffer_size), "cudaMalloc hits");
        check_cuda(cudaMalloc(&d_counter.ptr, sizeof(uint32_t)), "cudaMalloc global_hit_count");
        check_cuda(cudaMemset(d_counter.ptr, 0, sizeof(uint32_t)), "memset global_hit_count");

        MultiPatternArgs args;
        args.bv_A              = static_cast<const uint64_t*>(d_bvA.ptr);
        args.bv_C              = static_cast<const uint64_t*>(d_bvC.ptr);
        args.bv_G              = static_cast<const uint64_t*>(d_bvG.ptr);
        args.bv_T              = static_cast<const uint64_t*>(d_bvT.ptr);
        args.genome_start      = start;
        args.genome_len        = length;
        args.bv_word_offset    = first_word;
        args.num_patterns      = num_patterns_in_batch;
        args.threshold         = threshold;
        args.hits              = static_cast<SparseHit*>(d_hits.ptr);
        args.global_hit_count  = static_cast<uint32_t*>(d_counter.ptr);
        args.max_total_hits    = max_total_hits;
        args.stride            = stride;
        args.overlap           = overlap;

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

        int block_sz = choose_gpu_block();
        dim3 block(block_sz);
        dim3 grid((num_chunks + block_sz - 1) / block_sz);
        launch_myers_multi(use_32bit, batch_size, grid, block, args);
        check_cuda(cudaGetLastError(), "sparse multi-pattern kernel launch");
        check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

        uint32_t total_hits = 0;
        check_cuda(cudaMemcpy(&total_hits, d_counter.ptr,
                              sizeof(uint32_t), cudaMemcpyDeviceToHost),
                   "memcpy global_hit_count");
        if (total_hits > max_total_hits) total_hits = static_cast<uint32_t>(max_total_hits);
        if (total_hits == 0) continue;

        std::vector<SparseHit> all_hits(total_hits);
        check_cuda(cudaMemcpy(all_hits.data(), d_hits.ptr,
                              static_cast<size_t>(total_hits) * sizeof(SparseHit),
                              cudaMemcpyDeviceToHost),
                   "memcpy packed hits");

        std::vector<uint32_t> hits_per_pattern(num_patterns_in_batch, 0);
        for (uint32_t i = 0; i < total_hits; ++i) {
            uint8_t p = all_hits[i].pattern_idx;
            if (p < num_patterns_in_batch) hits_per_pattern[p]++;
        }
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            results[batch_start + p].hits.reserve(hits_per_pattern[p]);
        }
        for (uint32_t i = 0; i < total_hits; ++i) {
            const SparseHit& hit = all_hits[i];
            if (hit.pattern_idx >= num_patterns_in_batch) continue;
            size_t result_idx = batch_start + hit.pattern_idx;
            GpuSparseHit gh;
            gh.position = start + hit.position;
            gh.distance = hit.distance;
            results[result_idx].hits.push_back(gh);
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

    // Coalesced output: all patterns' hits are packed contiguously into
    // `hits` using a single batch-wide atomic counter. Each SparseHit carries
    // its pattern_idx so the host can demux after a single D2H transfer.
    SparseHit* hits;
    uint32_t*  global_hit_count;  // single counter for the whole batch
    size_t     max_total_hits;    // cap on total hits in `hits`

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


// Shift-add multi-pattern kernel (sparse output).
//
// Templated on word type W: uint64_t for patterns up to 64 bp, uint32_t for
// patterns ≤ 32 bp. The 32-bit variant halves per-thread state (512 B at
// threshold=3 instead of 1024 B), doubling the L2-safe thread count.
// Correct for short patterns because the bit-parallel window is 1 bit per
// pattern position — never indexes above (pattern_len - 1).

template <int MAX_D, typename W, int MAX_P>
__global__ void shift_add_multi_pattern_kernel_t(ShiftAddMultiPatternArgs args) {
    size_t tid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;

    size_t num_chunks = (args.genome_len + args.stride - 1) / args.stride;
    if (tid >= num_chunks) return;

    size_t out_start = tid * args.stride;
    size_t out_end   = out_start + args.stride;
    if (out_end > args.genome_len) out_end = args.genome_len;

    size_t proc_start = (tid == 0) ? 0 : (out_start - args.overlap);
    size_t proc_end   = out_end;

    // Compile-time-sized state: threshold+1 slots × MAX_P patterns × sizeof(W).
    // At threshold=3, W=uint32_t, MAX_P=32 this is 512 B per thread; at
    // MAX_P=16 it's 256 B. Caller batches no more than MAX_P patterns per
    // launch.
    W R[MAX_P][MAX_D];

    for (int p = 0; p < args.num_patterns; ++p) {
        #pragma unroll
        for (int d = 0; d < MAX_D; ++d) R[p][d] = W{0};
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
            // PM is stored as uint64_t on the host; narrowing to W just keeps
            // the low pattern_len bits (which is all a ≤32-nt pattern uses).
            W c_mask = valid ? static_cast<W>(args.PM[p * 4 + nuc_idx]) : W{0};

            W old_R    = R[p][0];
            W new_R0   = ((old_R << 1) | W{1}) & c_mask;
            R[p][0]    = new_R0;
            W hit_mask = new_R0;

            #pragma unroll
            for (int d = 1; d < MAX_D; ++d) {
                W temp   = R[p][d];
                W new_Rd = ((temp << 1) & c_mask) | (old_R << 1) | W{1};
                R[p][d]  = new_Rd;
                hit_mask |= new_Rd;
                old_R     = temp;
            }

            if (!valid || rel < out_start) continue;

            int m = args.pattern_lengths[p];
            size_t cbit = (rel < (size_t)(m - 1)) ? rel : (size_t)(m - 1);
            W check_mask = W{1} << cbit;

            if (!(hit_mask & check_mask)) continue;

            int distance = 0;
            #pragma unroll
            for (int d = 0; d < MAX_D; ++d) {
                if (R[p][d] & check_mask) { distance = d; break; }
            }

            // Coalesced write: single batch-wide atomic, pattern_idx embedded.
            uint32_t idx = atomicAdd(args.global_hit_count, 1u);
            if (idx < args.max_total_hits) {
                args.hits[idx].position    = static_cast<uint32_t>(rel);
                args.hits[idx].distance    = static_cast<uint8_t>(distance);
                args.hits[idx].pattern_idx = static_cast<uint8_t>(p);
            }
        }
    }
}

// Genome-chunk concurrency for the multi-pattern shift-add kernel.
// STOMATA_GENOME_CHUNKS=N splits the genome into N disjoint slices and launches
// N kernels per batch on N CUDA streams, each scanning its own slice. Disjoint
// slices → disjoint L2 working sets, unlike STOMATA_PIPELINE_STREAMS=1 which
// overlapped two kernels on the SAME genome bytes and contended on L2.
// Default 1 (single-stream). Supported: 1, 2, 4.
static int genome_chunks_count() {
    if (const char* ov = std::getenv("STOMATA_GENOME_CHUNKS")) {
        int v = std::atoi(ov);
        if (v == 1 || v == 2 || v == 4) return v;
    }
    return 1;
}

// Shift-add per-launch pattern batch size. Default 32.
// `STOMATA_SHIFTADD_BATCH=16` — half-batch (smaller L2 footprint, no default win)
// `STOMATA_SHIFTADD_BATCH=64` — double-batch (halves the batch count; useful at
//                              library scale where ⌈n/batch⌉ dominates wall time)
// In contrast to Myers (where half-batch is the default, see myers_batch_size),
// shift-add is already operating well under its L2 ceiling at MAX_P=32 — see
// §4.9 in docs/TOOL_BENCHMARKS_REPORT.md for measurements.
static int shift_add_batch_size() {
    if (const char* ov = std::getenv("STOMATA_SHIFTADD_BATCH")) {
        int v = std::atoi(ov);
        if (v == 16 || v == 32 || v == 64) return v;
    }
    return 32;
}

// Dispatch helper: picks the right kernel instantiation.
// MAX_D = threshold + 1; supported threshold range mirrors the 12-slot bound.
// `use_32bit` selects the uint32_t state variant (valid when every pattern in
// the batch is ≤ 32 bp); otherwise uses the uint64_t variant.
// `max_p` selects the per-thread R[MAX_P][MAX_D] state size (16 or 32).
static inline void launch_shift_add_multi(int threshold, bool use_32bit, int max_p,
                                          dim3 grid, dim3 block,
                                          cudaStream_t stream,
                                          const ShiftAddMultiPatternArgs& args) {
    #define STOMATA_DISPATCH_WP(W, P)                                                                      \
        switch (threshold + 1) {                                                                         \
            case 1:  shift_add_multi_pattern_kernel_t<1,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 2:  shift_add_multi_pattern_kernel_t<2,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 3:  shift_add_multi_pattern_kernel_t<3,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 4:  shift_add_multi_pattern_kernel_t<4,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 5:  shift_add_multi_pattern_kernel_t<5,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 6:  shift_add_multi_pattern_kernel_t<6,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 7:  shift_add_multi_pattern_kernel_t<7,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 8:  shift_add_multi_pattern_kernel_t<8,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 9:  shift_add_multi_pattern_kernel_t<9,  W, P><<<grid, block, 0, stream>>>(args); break; \
            case 10: shift_add_multi_pattern_kernel_t<10, W, P><<<grid, block, 0, stream>>>(args); break; \
            case 11: shift_add_multi_pattern_kernel_t<11, W, P><<<grid, block, 0, stream>>>(args); break; \
            case 12: shift_add_multi_pattern_kernel_t<12, W, P><<<grid, block, 0, stream>>>(args); break; \
            default: shift_add_multi_pattern_kernel_t<12, W, P><<<grid, block, 0, stream>>>(args); break; \
        }
    if (use_32bit) {
        if      (max_p <= 16) { STOMATA_DISPATCH_WP(uint32_t, 16) }
        else if (max_p <= 32) { STOMATA_DISPATCH_WP(uint32_t, 32) }
        else                  { STOMATA_DISPATCH_WP(uint32_t, 64) }
    } else {
        if      (max_p <= 16) { STOMATA_DISPATCH_WP(uint64_t, 16) }
        else if (max_p <= 32) { STOMATA_DISPATCH_WP(uint64_t, 32) }
        else                  { STOMATA_DISPATCH_WP(uint64_t, 64) }
    }
    #undef STOMATA_DISPATCH_WP
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
    size_t stride  = choose_gpu_stride(length, static_cast<size_t>(m),
                                        3, 8, GpuKernelKind::Hamming, 1);
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
    // Use the uint32_t kernel variant when every pattern fits in 32 bits of
    // state. For CRISPR spacers (17–24 nt) this is always true.
    const bool use_32bit = (max_pattern_len <= 32);
    const size_t word_bytes = use_32bit ? 4 : 8;
    const int batch_size = shift_add_batch_size();

    // Stride adapts to device L2 size (see choose_gpu_stride) so per-thread
    // local-memory state stays in cache. Pass actual batch_size so the stack
    // estimate scales with how many patterns this kernel instantiation holds.
    size_t stride = choose_gpu_stride(length, max_pattern_len, threshold,
                                       word_bytes, GpuKernelKind::Hamming,
                                       batch_size);
    size_t overlap = 2 * max_pattern_len;
    size_t num_chunks = (length + stride - 1) / stride;

    // Coalesced-output layout: every pattern's hits land in a single flat
    // buffer, counted by one batch-wide atomic. This reduces the per-batch
    // D2H transfer count from O(num_patterns) to 2 (one for the total count,
    // one for the packed hit array), eliminating ~30 cudaMemcpy launch
    // latencies per 225-spacer run.
    //
    // Cap the total buffer at MAX_PATTERNS_PER_KERNEL * MAX_HITS_PER_PATTERN
    // entries, matching the old worst-case footprint (3.84 GB).
    const size_t max_total_hits      = MAX_PATTERNS_PER_KERNEL * MAX_HITS_PER_PATTERN;
    const size_t hits_buffer_size    = max_total_hits * sizeof(SparseHit);
    const size_t counter_size        = sizeof(uint32_t);

    // Pipelined (2-stream) double-buffer path lets the host harvest batch N-1
    // while the GPU runs kernel N. Gated by env var; default remains the
    // single-stream path (the kernel is memory-bandwidth-bound, so stream
    // concurrency tends to regress).
    //   STOMATA_PIPELINE_STREAMS=1 → use double-buffered 2-stream pipeline
    //   STOMATA_GENOME_CHUNKS=N (N>1) → disjoint-slice concurrency (see below);
    //     takes precedence over pipeline and runs N kernels per batch on N
    //     streams. Each slot holds one chunk's hits + counter.
    const int n_chunks = genome_chunks_count();
    const bool chunked = (n_chunks > 1);
    const bool pipeline = !chunked && [](){
        const char* v = std::getenv("STOMATA_PIPELINE_STREAMS");
        return v && *v && *v != '0';
    }();
    const int N_SLOTS = chunked ? n_chunks : (pipeline ? 2 : 1);

    std::vector<DevicePtr> d_hits_slots(N_SLOTS);
    std::vector<DevicePtr> d_counter_slots(N_SLOTS);
    std::vector<cudaStream_t> streams(N_SLOTS, (cudaStream_t)0);
    const bool need_streams = (pipeline || chunked);
    for (int s = 0; s < N_SLOTS; ++s) {
        check_cuda(cudaMalloc(&d_hits_slots[s].ptr,    hits_buffer_size),
                   "cudaMalloc hits");
        check_cuda(cudaMalloc(&d_counter_slots[s].ptr, counter_size),
                   "cudaMalloc global_hit_count");
        if (need_streams) {
            check_cuda(cudaStreamCreate(&streams[s]), "cudaStreamCreate");
        }
    }

    // Pending-harvest state per slot: records which batch was launched on
    // each slot so we can harvest it when its stream finishes.
    //
    // For chunked mode, each slot represents a genome chunk of the SAME
    // batch; kernel_genome_start differs per slot (the chunk's genome
    // offset), and warmup_skip is the prefix of the kernel's output range
    // whose hits must be discarded (states at positions [0, warmup_skip)
    // lack real preceding-genome context in the shift-add state machine).
    struct PendingHarvest {
        bool   valid = false;
        int    num_patterns = 0;
        size_t batch_start = 0;
        size_t kernel_genome_start = 0;  // absolute; offset added to hit.position
        size_t warmup_skip = 0;          // hits with hit.position < this are dropped
    };
    std::vector<PendingHarvest> pending(N_SLOTS);

    auto harvest_slot = [&](int s) {
        if (!pending[s].valid) return;
        const int    num_patterns_in_batch = pending[s].num_patterns;
        const size_t batch_start_idx       = pending[s].batch_start;
        const size_t abs_base              = pending[s].kernel_genome_start;
        const size_t warmup_skip           = pending[s].warmup_skip;

        if (need_streams) {
            check_cuda(cudaStreamSynchronize(streams[s]), "cudaStreamSynchronize harvest");
        } else {
            check_cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize harvest");
        }

        // D2H #1: total hit count for this batch (4 bytes).
        uint32_t total_hits = 0;
        check_cuda(cudaMemcpy(&total_hits, d_counter_slots[s].ptr,
                              counter_size, cudaMemcpyDeviceToHost),
                   "memcpy global_hit_count");
        if (total_hits > max_total_hits) total_hits = static_cast<uint32_t>(max_total_hits);
        if (total_hits == 0) { pending[s].valid = false; return; }

        // D2H #2: the full packed hit array in one transfer.
        std::vector<SparseHit> all_hits(total_hits);
        check_cuda(cudaMemcpy(all_hits.data(), d_hits_slots[s].ptr,
                              static_cast<size_t>(total_hits) * sizeof(SparseHit),
                              cudaMemcpyDeviceToHost),
                   "memcpy packed hits");

        // Host-side demux: count per pattern (post warmup filter), reserve,
        // then push.
        std::vector<uint32_t> hits_per_pattern(num_patterns_in_batch, 0);
        for (uint32_t i = 0; i < total_hits; ++i) {
            if (all_hits[i].position < warmup_skip) continue;
            uint8_t p = all_hits[i].pattern_idx;
            if (p < num_patterns_in_batch) hits_per_pattern[p]++;
        }
        // NOTE: in chunked mode, this reserve races across chunks writing to
        // the same result vector (different chunks can contribute hits to the
        // same pattern). reserve() only grows, so multiple .reserve(n) calls
        // are fine; the .push_back() calls are ordered per harvest and run on
        // the host thread sequentially across chunk harvests → safe.
        for (int p = 0; p < num_patterns_in_batch; ++p) {
            auto& hv = results[batch_start_idx + p].hits;
            hv.reserve(hv.size() + hits_per_pattern[p]);
        }
        for (uint32_t i = 0; i < total_hits; ++i) {
            const SparseHit& hit = all_hits[i];
            if (hit.position < warmup_skip) continue;
            if (hit.pattern_idx >= num_patterns_in_batch) continue;
            size_t result_idx = batch_start_idx + hit.pattern_idx;
            GpuSparseHit gh;
            gh.position = abs_base + hit.position;
            gh.distance = hit.distance;
            results[result_idx].hits.push_back(gh);
        }
        pending[s].valid = false;
    };

    // Helper: fill pattern masks into an args struct for a given pattern slice.
    auto fill_pattern_args = [&](ShiftAddMultiPatternArgs& args,
                                 size_t batch_start, int num_in_batch) {
        std::memset(args.PM, 0, sizeof(args.PM));
        std::memset(args.pattern_lengths, 0, sizeof(args.pattern_lengths));
        for (int p = 0; p < num_in_batch; ++p) {
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
    };

    const int block_sz = choose_gpu_block();
    dim3 block(block_sz);

    if (chunked) {
        // Disjoint-genome-chunk concurrency: each batch fires N_SLOTS kernels
        // on N_SLOTS streams, each scanning its own [out_start, out_end) slice
        // with a 2*max_pattern_len left-warmup (except chunk 0, which starts
        // at the genome's true start and has no preceding context).
        std::vector<size_t> chunk_out_start(N_SLOTS + 1);
        for (int c = 0; c <= N_SLOTS; ++c) {
            chunk_out_start[c] = (length * static_cast<size_t>(c)) /
                                 static_cast<size_t>(N_SLOTS);
        }

        for (size_t batch_start = 0; batch_start < patterns.size();
             batch_start += batch_size) {
            size_t batch_end = std::min(batch_start + (size_t)batch_size,
                                        patterns.size());
            int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);

            for (int c = 0; c < N_SLOTS; ++c) {
                const size_t out_start_c = chunk_out_start[c];
                const size_t out_end_c   = chunk_out_start[c + 1];
                // Bug #3 fix (v0.9.0): for chunk 0 we previously used 0,
                // which left the bitap warm-up region (positions 0..m-2 of
                // the search range) unfiltered. The shift-add state machine
                // reports prefix-match distances at those positions — partial
                // pattern alignments that aren't real hits. The single-
                // pattern path filters these via min_valid_pos in
                // search_pipeline.cpp; the batch path was missing the
                // equivalent. Use max_pattern_len - 1 for chunk 0 so the
                // host-side warmup_skip filter at lines below correctly
                // discards them.
                //
                // For chunks > 0, `overlap` (= 2 * max_pattern_len) covers
                // both warm-up and the inter-chunk seam boundary.
                const size_t warmup_c    = (c == 0)
                                           ? (max_pattern_len > 0 ? max_pattern_len - 1 : 0)
                                           : overlap;
                const size_t kernel_start = (start + out_start_c >= warmup_c)
                                            ? start + out_start_c - warmup_c
                                            : 0;
                const size_t kernel_len   = out_end_c - out_start_c + (start + out_start_c - kernel_start);
                if (kernel_len == 0) continue;

                check_cuda(cudaMemsetAsync(d_counter_slots[c].ptr, 0,
                                           counter_size, streams[c]),
                           "memset global_hit_count");

                ShiftAddMultiPatternArgs args;
                args.bv_A             = static_cast<const uint64_t*>(d_bvA.ptr);
                args.bv_C             = static_cast<const uint64_t*>(d_bvC.ptr);
                args.bv_G             = static_cast<const uint64_t*>(d_bvG.ptr);
                args.bv_T             = static_cast<const uint64_t*>(d_bvT.ptr);
                args.genome_start     = kernel_start;
                args.genome_len       = kernel_len;
                args.bv_word_offset   = first_word;
                args.num_patterns     = num_patterns_in_batch;
                args.threshold        = threshold;
                args.hits             = static_cast<SparseHit*>(d_hits_slots[c].ptr);
                args.global_hit_count = static_cast<uint32_t*>(d_counter_slots[c].ptr);
                args.max_total_hits   = max_total_hits;
                args.stride           = stride;
                args.overlap          = overlap;
                fill_pattern_args(args, batch_start, num_patterns_in_batch);

                size_t num_chunks_c = (kernel_len + stride - 1) / stride;
                dim3 grid_c((num_chunks_c + block_sz - 1) / block_sz);
                launch_shift_add_multi(threshold, use_32bit, batch_size,
                                       grid_c, block, streams[c], args);
                check_cuda(cudaGetLastError(),
                           "shift_add multi-pattern kernel launch (chunked)");

                pending[c] = { true, num_patterns_in_batch, batch_start,
                               kernel_start, warmup_c };
            }

            // Harvest every chunk for this batch before moving on.
            for (int c = 0; c < N_SLOTS; ++c) harvest_slot(c);
        }
    } else {
        // Single-stream or pipeline-overlap path (legacy).
        int cur_slot = 0;

        for (size_t batch_start = 0; batch_start < patterns.size();
             batch_start += batch_size) {

            // If the slot we're about to use still has work in flight from an
            // earlier iteration (only possible when pipelining), harvest it
            // first.
            if (pending[cur_slot].valid) {
                harvest_slot(cur_slot);
            }

            size_t batch_end = std::min(batch_start + (size_t)batch_size,
                                        patterns.size());
            int num_patterns_in_batch = static_cast<int>(batch_end - batch_start);

            // Zero the global counter for this slot.
            if (pipeline) {
                check_cuda(cudaMemsetAsync(d_counter_slots[cur_slot].ptr, 0,
                                            counter_size, streams[cur_slot]),
                            "memset global_hit_count");
            } else {
                check_cuda(cudaMemset(d_counter_slots[cur_slot].ptr, 0, counter_size),
                            "memset global_hit_count");
            }

            ShiftAddMultiPatternArgs args;
            args.bv_A             = static_cast<const uint64_t*>(d_bvA.ptr);
            args.bv_C             = static_cast<const uint64_t*>(d_bvC.ptr);
            args.bv_G             = static_cast<const uint64_t*>(d_bvG.ptr);
            args.bv_T             = static_cast<const uint64_t*>(d_bvT.ptr);
            args.genome_start     = start;
            args.genome_len       = length;
            args.bv_word_offset   = first_word;
            args.num_patterns     = num_patterns_in_batch;
            args.threshold        = threshold;
            args.hits             = static_cast<SparseHit*>(d_hits_slots[cur_slot].ptr);
            args.global_hit_count = static_cast<uint32_t*>(d_counter_slots[cur_slot].ptr);
            args.max_total_hits   = max_total_hits;
            args.stride           = stride;
            args.overlap          = overlap;
            fill_pattern_args(args, batch_start, num_patterns_in_batch);

            dim3 grid((num_chunks + block_sz - 1) / block_sz);
            launch_shift_add_multi(threshold, use_32bit, batch_size, grid, block,
                                   streams[cur_slot], args);
            check_cuda(cudaGetLastError(), "shift_add multi-pattern kernel launch");

            // Bug #3 fix (v0.9.0): warmup_skip = max_pattern_len - 1 to drop
            // bitap prefix-match artifacts from the search-range start.
            // Single-pattern path filters via min_valid_pos in
            // search_pipeline.cpp:797; batch path was missing the equivalent
            // and emitted ~m-1 spurious "partial-prefix" Hamming hits per
            // chromosome at distance ≪ threshold. See VERIFICATION_MATRIX.md.
            const size_t warmup_skip_legacy =
                (max_pattern_len > 0) ? max_pattern_len - 1 : 0;
            pending[cur_slot] = { true, num_patterns_in_batch, batch_start,
                                  start, warmup_skip_legacy };

            if (pipeline) {
                cur_slot = (cur_slot + 1) % N_SLOTS;
            } else {
                harvest_slot(cur_slot);
            }
        }

        // Drain any remaining pending harvests (pipeline mode)
        for (int s = 0; s < N_SLOTS; ++s) {
            if (pending[s].valid) harvest_slot(s);
        }
    }

    // Destroy streams
    if (need_streams) {
        for (int s = 0; s < N_SLOTS; ++s) {
            if (streams[s]) check_cuda(cudaStreamDestroy(streams[s]), "cudaStreamDestroy");
        }
    }

    return results;
}

