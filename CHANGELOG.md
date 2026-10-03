# Changelog

## 1.0.0 — 2026-10-03

First stable release.

**Search**
- Exhaustive Hamming or Levenshtein (mismatch + indel) search of a reference
  genome for one guide or a batch, up to a chosen distance threshold.
- Any IUPAC PAM at the 3′ or 5′ end, or no PAM restriction; both strands or
  one.
- Memory-mapped `.st` reference index; optional region restriction.
- NVIDIA CUDA acceleration, with a CPU-only build and runtime fallback.

**Output**
- Per-hit TSV, BED or JSON with coordinates, strand, alignment, CIGAR and PAM.
- CFD scores, per-guide summaries with aggregate CFD, and BED overlap counts.

**Validation**
- Unit and CLI integration tests (CTest), run in CI for the CPU build and the
  container image.
- Independent dynamic-programming oracle with fixtures, and randomized fuzz
  comparisons against it on both the CPU and GPU paths.
- Cross-tool comparisons with Cas-OFFinder and SWOffinder, and recovery of
  experimentally observed CHANGE-seq and GUIDE-seq sites
  (see [validation/](validation/README.md)).
