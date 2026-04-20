# Cystidia Performance Baseline

**Latest update:** 2026-04-19 (GPU kernel optimization sweep, v0.5.0)
**Test:** 225 spacers from cross-tool validation set
**Configuration:** threshold=3, hg38.p14 canonical genome, `.cy` mmap index
**Hardware:** GPU-enabled system

## GPU Kernel Comparison (v0.5.0, search-only, `--pam-filter none --no-compute-mismatches`)

| Configuration | Runtime | Notes |
|---|---|---|
| Cystidia `--distance-mode hamming` (shift-add, default) | **38.1s** | Fastest Cystidia config. Templated R-array, OR-reduce early hit detection. |
| Cystidia `--distance-mode levenshtein` (Myers) | 47.9s | Substitutions + indels. |
| Cas-OFFinder, NGG PAM prefilter | 21.2s | Fast due to PAM-coupled search. |
| Cas-OFFinder, PAM-agnostic (`NNNNNNNNNNNNNNNNNNNNNNN`) | ❌ `CL_OUT_OF_RESOURCES` | Crashes at this scale. |

**Headline:** On the apples-to-apples PAM-agnostic comparison, Cas-OFFinder cannot complete the run; Cystidia finishes in 38s. Cas-OFFinder's 21s number depends on the PAM prefilter rejecting >90% of genome positions before any mismatch counting.

### Kernel variants evaluated and rejected

During the optimization sweep we implemented and benchmarked two alternatives, kept the data honest, and then stripped them from the CLI:

- **Scalar Hamming with per-position early exit** (Cas-OFFinder-style): 101.8s — 2.7× slower than shift-add. GPU-memory divergence and poor coalescing dominate what the early-exit saves.
- **CUDA stream pipelining (2-stream ping-pong)**: 37.9s vs 38.1s plain — within run-to-run noise. Genome H2D is amortized once across the batch, so there is no meaningful overlap to harvest.

Both paths are removed from v0.5.0. See `src/gpu_engine.cu` for the surviving kernels.

---

## Earlier Phase 1A Baseline (2026-02-13, raw FASTA, before shift-add optimization)

## Benchmark Results

| Mode | Runtime | Hits Output | Per-spacer | Annotation Overhead | Notes |
|------|---------|-------------|------------|---------------------|-------|
| **Full annotation + PAM filter NGG** | 155.9s | 278,581 | 0.693s | 104.8s (67.2%) | Current default mode |
| **Summary mode** | 166.5s | Aggregate counts only | 0.740s | 115.4s (69.3%) | ⚠️ SLOWER than full! |
| **Search-only (no annotations)** | 51.1s | 4,094,195 | 0.227s | 0s (baseline) | Pure search performance |

## Key Findings

### 1. Annotation is the Bottleneck ✅

**Search phase:** 51.1 seconds (32.8% of full annotation runtime)
**Annotation overhead:** 104.8 seconds (67.2% of full annotation runtime)

This confirms the bottleneck identified in Investigation 9 of CROSS_TOOL_VALIDATION.md.

### 2. Summary Mode is NOT Optimized ⚠️

**Expected:** Summary mode should be fastest (no per-hit details)
**Actual:** Summary mode is 6.8% SLOWER than full annotation mode (166.5s vs 155.9s)

**Root cause** (from code inspection):
- Summary mode runs the full search pipeline with complete annotations
- Then aggregates results into count summaries via `summarize_batch_results()`
- It's a post-processing approach, NOT a performance optimization
- Still computes all expensive alignments, PAM extractions, CFD scores

**Example summary output:**
```json
{
  "threshold": 3,
  "total_spacers_processed": 225,
  "spacers_skipped": 0,
  "spacers": [
    {
      "name": "AAVS1_site_1",
      "sequence": "GTCACCAATCCTGTCCCTAG",
      "total_hits": 1027,
      "hits_by_distance": {"0": 1, "2": 34, "3": 992}
    },
    ...
  ]
}
```

### 3. PAM Filtering Effectiveness

**Without PAM filter (--pam-filter none):** 4,094,195 hits
**With NGG filter (--pam-filter ngg):** 278,581 hits
**Reduction:** 93.2% of hits filtered out

This shows PAM filtering is highly selective, which validates the lazy PAM extraction optimization (Priority 2) - if we extract PAM first and filter early, we can skip expensive annotation for 93% of hits.

### 4. Performance vs CasOFFinder

**Cystidia full annotation:** 155.9s (0.693s/spacer)
**CasOFFinder GPU (from Investigation 9):** 42s (0.187s/spacer)
**Gap:** 113.9s (Cystidia is 3.7x slower)

**Cystidia search-only:** 51.1s (0.227s/spacer)
**CasOFFinder GPU:** 42s (0.187s/spacer)
**Gap:** 9.1s (Cystidia is 1.22x slower)

**Analysis:**
- 92% of the performance gap (104.8s of 113.9s) is due to annotation overhead
- Only 8% of the gap (9.1s) is due to search phase differences
- Annotation optimizations can close most of the gap
- Matching CasOFFinder's 42s requires BOTH annotation optimization AND GPU kernel optimization

### 5. Hit Distribution Analysis

Examining search-only output shows:
- **With PAM filter:** Mostly NGG/NAG sites (278K hits)
- **Without PAM filter:** All PAM types including OTHER (4.09M hits)
- Even with `--no-compute-mismatches`, PAM classification still happens (see output format)

This suggests `--no-compute-mismatches` doesn't fully disable annotation - it may still extract PAM sequences for classification.

## Conclusions & Recommendations

### ✅ Proceed to Phase 1B: Implementation

The benchmarking confirms:

1. **Annotation is the bottleneck** (105s of 156s runtime = 67%)
2. **Search-only mode achieves ~51s** (close to predicted 48s)
3. **Summary mode needs optimization** (currently slower than full mode)
4. **Lazy PAM extraction will be highly effective** (93% of hits filtered by PAM)

### Projected Impact of Optimizations

Based on these benchmarks, expected improvements:

**Priority 1: Minimal Output Mode**
- Skip: DP alignment, aligned_seq building, edit type classification
- Keep: PAM extraction only (cheap)
- **Projected runtime:** 51s (search) + 15-20s (PAM only) = **66-71s**
- **Speedup:** 54-58% faster than baseline

**Priority 2: Lazy PAM Extraction**
- Extract PAM first (cheap: ~15s for 4.09M hits)
- Filter early (keeps 278K of 4.09M hits = 6.8%)
- Full annotation only for passing hits (6.8% of work)
- **Projected runtime:** 51s (search) + 15s (all PAM) + 7s (6.8% annotation) = **73s**
- **Speedup:** 53% faster than baseline

**Combined (Minimal + Lazy):**
- **Projected runtime:** 51s (search) + 15s (all PAM) + 1s (6.8% minimal annotation) = **67s**
- **Speedup:** 57% faster than baseline
- **vs CasOFFinder:** 67s vs 42s (1.6x slower, down from 3.7x)

**Priority 4: Parallel Annotation (8 cores)**
- Parallelize the 15s PAM extraction across 8 cores
- **Projected runtime:** 51s (search) + 2s (parallel PAM) = **53s**
- **Speedup:** 66% faster than baseline
- **vs CasOFFinder:** 53s vs 42s (1.26x slower)

### Next Steps

1. ✅ **Phase 1A complete** - Baseline metrics established
2. 🔄 **Begin Phase 1B** - Implement Priority 1 (Minimal Output Mode)
3. 🔄 **Implement Priority 2** - Lazy PAM extraction with early filtering
4. 📊 **Re-benchmark** - Measure actual speedup vs projections
5. 🎯 **Evaluate** - Decide whether to proceed to Phase 2 (parallel annotation)

## Comparison with Investigation 9

Investigation 9 reported:
- Full annotation: 152s
- Search-only: 48s
- Annotation overhead: 104s

Our Phase 1A results:
- Full annotation: 155.9s ✅ (within 2.6% of Investigation 9)
- Search-only: 51.1s ✅ (within 6.5% of Investigation 9)
- Annotation overhead: 104.8s ✅ (within 0.8% of Investigation 9)

**Conclusion:** Results are highly consistent with Investigation 9. Projections are validated.

---

## Appendix: Benchmark Commands

```bash
# Create output directory
mkdir -p ~/tmp/argus_bench

# Full annotation mode (baseline)
time ./build/src/cystidia --spacer-file tests/data/validation/cross_tool_spacers.txt \
  --genome path/to/genomes/GRCh38.p14/hg38.p14.canonical.fa.cy \
  --threshold 3 --pam-filter ngg \
  --output ~/tmp/argus_bench/full_annotation.tsv --quiet

# Summary mode
time ./build/src/cystidia --spacer-file tests/data/validation/cross_tool_spacers.txt \
  --genome path/to/genomes/GRCh38.p14/hg38.p14.canonical.fa.cy \
  --threshold 3 --summary --summary-format json \
  --output ~/tmp/argus_bench/summary.json --quiet

# Search-only mode (no annotations)
time ./build/src/cystidia --spacer-file tests/data/validation/cross_tool_spacers.txt \
  --genome path/to/genomes/GRCh38.p14/hg38.p14.canonical.fa.cy \
  --threshold 3 --pam-filter none --no-compute-mismatches \
  --output ~/tmp/argus_bench/search_only.tsv --quiet
```

## Benchmark Output Files

- Full annotation: 278,581 lines (278K hits with NGG PAM filter)
- Search-only: 4,094,195 lines (4.09M hits without PAM filter)
- Summary: 1,357 lines JSON (aggregate counts only)
