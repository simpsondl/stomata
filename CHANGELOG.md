# Changelog

All notable changes to Cystidia are documented here. Versioning follows semantic
versioning (major.minor.patch). Kernel-correctness fixes may change result
sets; these are called out in Bug fixes.

## [0.6.0] — 2026-04-22

### Added
- **PAM agnosticism.** The search pipeline is no longer hardcoded to SpCas9/NGG.
  `--pam` accepts any IUPAC PAM pattern and `--pam-position` selects 3' (Cas9)
  or 5' (Cas12) placement. `--pam-extract-length` controls annotated PAM width.
- **Biological annotations per hit.** CFD activity score (Doench et al., 2016),
  PAM classification (canonical NGG / NAG / non-canonical), seed-region
  disruption (SDR) site counts, mismatch positions, and CIGAR-format alignment
  strings.
- **`--spacer-summary`** flag emits a per-spacer aggregate TSV with:
  hit counts by edit distance (d0, d1, d2, d3+), aggregate CFD (aCFD)
  promiscuity score, seed-disruption-region-1 counts, and optional BED-file
  overlap counts.
- **`--intersect-bed LABEL:FILE`** (repeatable) — for each spacer, count
  hits that overlap intervals in the provided BED. Each label appears as a
  `n_<label>_overlaps` column in the spacer summary.
- **`--acfd-threshold`** configures the aCFD cutoff used for the
  `is_promiscuous` flag (default 4.8).
- **Step-specific timing** in `--verbose`: prints kernel time, CFD scoring
  time, and TSV formatting time separately.

### Changed
- Search pipeline refactored for PAM agnosticism. Callers that previously
  relied on implicit SpCas9 behavior should pass `--pam NGG` explicitly.
- CFD scorer interface adjusted to accept variable-length PAM contexts.

### Fixed
- **Kernel deduplication bug.** The GPU halo-deduplication logic was
  incorrectly collapsing distinct adjacent hits into one. Results from
  v0.5.0 may under-count hits in high-density regions; v0.6.0 restores
  full hit enumeration. Run `--no-deduplicate` to match the v0.5.0
  behavior if needed for reproducibility.

### Performance
- Added detailed per-stage timing instrumentation; no throughput regression.

## [0.5.0] — 2026-04-20

Initial public release.

### Added
- GPU-accelerated exhaustive edit-distance (Myers bit-parallel) off-target
  search for CRISPR guide RNAs.
- Hamming and Levenshtein distance modes.
- Sparse output architecture: per-hit atomic writes, constant-memory GPU
  output regardless of batch size.
- Single-pass batch search across multiple spacers in one kernel invocation.
- Memory-mapped binary genome index (`.cy`) for fast genome loads.
- Command-line interface with batch-mode, strand selection, and configurable
  edit-distance threshold.
- Cross-tool validation harness against Cas-OFFinder.
