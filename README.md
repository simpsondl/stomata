# Cystidia

GPU-accelerated exhaustive CRISPR off-target search.

**Version:** 0.6.0 · **License:** MIT · **Requires:** CUDA 11.8+, C++17, CMake 3.20+

## What it does

Given a guide RNA (or a batch of them) and a reference genome, Cystidia enumerates every position in the genome within a user-specified edit distance. Output is a TSV with per-hit coordinates, strand, CIGAR, PAM, CFD activity score, and seed/distal mismatch breakdown — or an aggregate per-spacer summary.

- **Exhaustive, not heuristic.** No seeds, no pruning. Every position is evaluated.
- **PAM-agnostic search, PAM-annotated output.** Cystidia finds all sequence matches first, then annotates PAM. You can re-filter without re-searching and support any IUPAC PAM (SpCas9 NGG, Cas12a TTTV, engineered variants).
- **Indel-aware.** Levenshtein search is a first-class path, not post-hoc.
- **GPU-resident genome.** The `.cy` index is memory-mapped; the first search pays the load cost, subsequent searches are instant.

## Current benchmarks (hg38, 225 spacers, threshold 3, search-only)

| Tool                    | Config                                 | Runtime  |
|-------------------------|----------------------------------------|----------|
| Cystidia                | `--distance-mode hamming` (shift-add)  | **38.1s** |
| Cystidia                | `--distance-mode levenshtein` (Myers)  | 47.9s    |
| Cas-OFFinder            | NGG PAM prefilter                      | 21.2s    |
| Cas-OFFinder            | PAM-agnostic (all-N)                   | crashes (`CL_OUT_OF_RESOURCES`) |

Run date: 2026-04-19, hardware: NVIDIA GPU (compute ≥7.0), hg38.p14 canonical, `.cy` mmap index. Full methodology and evaluated-but-rejected kernel variants (scalar Hamming, CUDA stream ping-pong): [PERFORMANCE_BASELINE.md](docs/PERFORMANCE_BASELINE.md).

Cas-OFFinder is faster on the canonical NGG case because the PAM prefilter rejects >90% of positions before any mismatch counting. It cannot run PAM-agnostic at this scale. For non-canonical PAMs, indel-aware search, or rich annotations, Cystidia is the only tool that finishes.

## Getting started

```bash
# 1. Build
conda create -n cystidia python=3.10 -y && conda activate cystidia
conda install -c conda-forge cmake compilers gxx=11 spdlog catch2 -y
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 2. Index a genome (one-time; ~30s for hg38)
./build/src/cystidia --index-genome hg38.fa   # produces hg38.fa.cy

# 3. Search
./build/src/cystidia \
  --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam NGG \
  --output hits.tsv
```

Batch mode: replace `--pattern` with `--spacer-file guides.txt` (one sequence per line, optional `name<TAB>seq` format, `#` comments allowed). Use `--spacer-summary summary.tsv` for per-spacer aggregates (hit counts by distance, aCFD promiscuity score, BED overlap counts).

## Documentation

- [User Guide](docs/USER_GUIDE.md) — installation, full CLI, output formats, troubleshooting
- [Performance Baseline](docs/PERFORMANCE_BASELINE.md) — current kernel benchmarks and the optimization decisions behind them
- [Benchmark Results](docs/BENCHMARK_RESULTS.md) — older end-to-end performance data (v0.2.0-era; some numbers predate the current search pipeline)
- `--help` on the binary prints the full flag reference

## Things that will surprise you

- **Default threshold is 4, not 3.** Pass `--threshold N` explicitly if you care.
- **Default PAM filter is none.** Without `--pam NGG` (or similar) the output includes every sequence match regardless of PAM context. This is intentional — PAM is an annotation, not a search constraint — but it means a first run produces more rows than Cas-OFFinder does.
- **GPU is automatic if present.** Pass `--cpu-only` to force CPU. Pattern length must be 1–64 bp for GPU; longer patterns fall back to multi-word Myers on CPU.
- **Hamming mode searches by Hamming distance but aligns by edit distance.** So `--distance-mode hamming` output may still contain indels in the alignment column — the search found the position via mismatch counting, but the reported CIGAR is the optimal edit-distance alignment at that position. If you need strict mismatch-only output, post-filter `BULGE` rows or use `--no-compute-mismatches` to skip alignment entirely.
- **U→T normalization is on by default.** RNA spacers work out of the box. Disable with `--no-treat-u-as-t` if you want strict validation.
- **`.cy` index files are memory-mapped, not copied.** First search pays the page-fault cost; subsequent searches on the same file are effectively free.
- **PAM-agnostic output is big.** On hg38, threshold 3, no PAM filter produces ~4M hits across 225 spacers vs ~280K with `--pam NGG`. Use `--max-hits` / `--max-total-hits` to cap, or `--summary` for aggregates only.

## Citation

```bibtex
@software{cystidia2026,
  author  = {Simpson, Danny and Sanjana, Neville E. and Lappalainen, Tuuli},
  title   = {Cystidia: GPU-accelerated exhaustive CRISPR off-target search},
  year    = {2026},
  version = {0.6.0},
  url     = {https://github.com/simpsondl/cystidia}
}
```

Issues: https://github.com/simpsondl/cystidia/issues
