# Stomata

GPU-accelerated exhaustive CRISPR off-target search.

**Version:** 0.10.0 · **License:** MIT · **Requires:** CUDA 11.8+, C++17, CMake 3.20+

## What it does

Given a guide RNA (or a batch of them) and a reference genome, Stomata enumerates every position in the genome within a user-specified edit distance. Output is a TSV with per-hit coordinates, strand, CIGAR, PAM, CFD activity score, and seed/distal mismatch breakdown — or an aggregate per-spacer summary.

- **Exhaustive.** Every position is evaluated.
- **PAM-agnostic search.** Stomata finds all sequence matches first, then annotates PAM. You can re-filter without re-searching and support any IUPAC PAM (SpCas9 NGG, Cas12a TTTV, engineered variants).
- **Indel-aware.** Levenshtein search implemented through Myers' bit-parallel algorithm.
- **GPU-resident genome.** The `.st` index is memory-mapped, allowing instant search post-load.

## Getting started

```bash
# 1. Build
conda create -n stomata python=3.10 -y && conda activate stomata
conda install -c conda-forge cmake compilers gxx=11 spdlog catch2 -y
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 2. Index a genome (one-time; typically 30 seconds to 5 minutes, hardware-dependent, for hg38)
./build/src/stomata --index-genome hg38.fa   # produces hg38.fa.st

# 3. Search
./build/src/stomata \
  --genome hg38.fa.st \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam NGG \
  --output hits.tsv
```

Batch mode: replace `--pattern` with `--spacer-file guides.txt` (one sequence per line, optional `name<TAB>seq` format, `#` comments allowed). 

Use `--spacer-summary summary.tsv` for per-spacer aggregates (hit counts by distance, CRISPick aCFD promiscuity score, BED overlap counts).

## Documentation

- [User Guide](docs/USER_GUIDE.md) — installation, full CLI, output formats, troubleshooting
- `--help` on the binary prints the full flag reference

## Things that will surprise you

- **Default threshold is 3.** Pass `--threshold N` explicitly if you care.
- **Default PAM filter is none.** Without `--pam NGG` (or similar) the output includes every sequence match regardless of PAM context. This is intentional but it means a first run produces more rows than Cas-OFFinder does. Use `--max-hits` / `--max-total-hits` to cap, or `--summary` for aggregates only.
- **GPU is automatic if present.** Pass `--cpu-only` to force CPU. Pattern length must be 1–64 bp for GPU; longer patterns fall back to multi-word Myers on CPU.
- **Hamming mode searches by Hamming distance but aligns by edit distance.** So `--distance-mode hamming` output may still contain indels in the alignment column — the search found the position via mismatch counting, but the reported CIGAR is the optimal edit-distance alignment at that position. If you need strict mismatch-only output, post-filter rows whose `cigar` contains `I` or `D`, or use `--no-compute-mismatches` to skip alignment entirely.
- **U→T normalization is on by default.** RNA spacers work out of the box. Disable with `--no-treat-u-as-t` if you want strict validation.

Issues: https://github.com/simpsondl/stomata/issues
