# Stomata User Guide

Complete guide to using Stomata for exhaustive CRISPR off-target search.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Basic Usage](#basic-usage)
4. [Distance Modes](#distance-modes)
5. [PAM Filtering](#pam-filtering)
6. [Strand Control](#strand-control)
7. [Region Restriction](#region-restriction)
8. [Spacer Summaries](#spacer-summaries)
9. [Output Formats](#output-formats)
10. [Performance](#performance)
11. [Troubleshooting](#troubleshooting)

---

## Installation

See [installation](INSTALL.md) for the CPU-only build and container.

### System Requirements

- **Operating System:** Linux (tested on Ubuntu 20.04+)
- **CPU:** x86_64
- **GPU:** NVIDIA GPU with CUDA compute capability ≥7.0 (optional but recommended)
- **Memory:** 8 GB RAM minimum (16 GB+ recommended for human genome)
- **Storage:** 10 GB for software + genome indices

### Dependencies

- CMake ≥3.20
- C++17 compiler (GCC 11+ or Clang 14+)
- CUDA Toolkit 11.8+ for a GPU build; omit it with `-DSTOMATA_ENABLE_CUDA=OFF`
- spdlog, zlib, Catch2 (Catch2 only for the test suite)

### Build (conda recommended)

```bash
git clone https://github.com/simpsondl/stomata.git
cd stomata

conda create -n stomata python=3.10 -y
conda activate stomata
conda install -c conda-forge cmake compilers gxx=11 spdlog catch2 -y

cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

./build/src/stomata --version
```

### Verify the install

```bash
./build/src/stomata --quickstart
```

`--quickstart` runs a self-contained smoke test against a tiny embedded
genome under four configurations (Hamming/Levenshtein × NGG/PAM-agnostic)
and exits non-zero if any hit count diverges from the engineered ground
truth. It is the fastest way to confirm GPU dispatch and CFD scoring work
on your hardware.

---

## Quick Start

### 30-second demo

```bash
./build/src/stomata \
  --genome tests/data/synthetic/test_small.fa \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3
```

### Real-world: human genome

```bash
# 1. Index hg38 once (~30 s – 5 min depending on hardware). Produces hg38.fa.st.
./build/src/stomata --index-genome hg38.fa

# 2. Search a single guide with SpCas9 NGG PAM filter
./build/src/stomata \
  --genome hg38.fa.st \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam NGG \
  --output vegfa_offtargets.tsv

# 3. Batch-search many guides (loads genome once, searches all)
./build/src/stomata \
  --genome hg38.fa.st \
  --spacer-file my_guides.txt \
  --threshold 3 \
  --pam NGG \
  --output batch_results.tsv
```

---

## Basic Usage

```
stomata --genome GENOME (--pattern SEQ | --spacer-file FILE) [OPTIONS]
stomata --index-genome FASTA               # one-time index build
stomata --quickstart                       # self-test
```

### Required arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--genome FILE` | FASTA, FASTA.gz, or `.st` index | `hg38.fa.st` |
| `--pattern SEQ` *or* `--spacer-file FILE` | Guide RNA(s) to search | `GAGTCC...` |

### Common options

| Option | Description | Default |
|--------|-------------|---------|
| `--threshold N` | Maximum edit distance to report (0–255) | `3` |
| `--distance-mode MODE` | `levenshtein` or `hamming` | `levenshtein` |
| `--pam PATTERN` | IUPAC PAM filter (e.g. `NGG`, `TTTV`) | none |
| `--pam-position POS` | `3prime` (Cas9) or `5prime` (Cas12) | `3prime` |
| `--pam-extract-length N` | Bases to extract for the PAM column | length of `--pam` |
| `--strand MODE` | `both`, `plus`, `minus` | `both` |
| `--region CHR[:START-END]` | Restrict search to one window (1-based, inclusive-exclusive) | whole genome |
| `--format FMT` | `tsv`, `bed`, or `json` | `tsv` |
| `--output FILE` | Write results to FILE instead of stdout | stdout |
| `--max-hits N` | Cap hits per spacer (`0` = unlimited) | unlimited |
| `--max-total-hits N` | Batch mode: cap total hits across spacers | unlimited |
| `--no-compute-mismatches` | Skip CIGAR / aligned-seq / PAM annotation | annotate |
| `--no-scores` | Skip CFD activity scoring | score |
| `--no-treat-u-as-t` | Reject `U` instead of normalizing to `T` | normalize |
| `--cpu-only` | Force CPU even if a GPU is available | use GPU if present |
| `--threads N` (`-t`) | Threads for batch processing (`0` = auto) | auto |
| `--summary` | Aggregate counts by edit distance instead of per-hit rows | per-hit |
| `--summary-format FMT` | `json` (default) or `tsv` for `--summary` | `json` |
| `--verbose` / `--quiet` | Stderr verbosity | normal |
| `--version` / `--help` / `--quickstart` | Info / self-test | — |

### Spacer file format

```
# Lines starting with # are comments; blank lines are ignored.
GAGTCCGAGCAGAAGAAGAA              # Sequence only (auto-named spacer_1, spacer_2, ...)
EMX1     GAGTCCGAGCAGAAGAAGAA     # Name<TAB>sequence
FANCF    GGAATCCCTTCTGCAGCACC NGG # Extra columns are ignored
```

Use `-` as the path to read from stdin.

---

## Distance Modes

Stomata supports two distance metrics for the search step, selected with
`--distance-mode`:

| Mode | Algorithm | Penalizes | Notes |
|------|-----------|-----------|-------|
| `levenshtein` (default) | Myers' bit-vector | substitutions + indels | Indel-aware, exhaustive |
| `hamming` | Shift-add | substitutions only | Faster, no indels |

**Hamming mode searches by Hamming distance but reports the alignment as
edit distance.** The position is found by mismatch counting; the CIGAR /
aligned-seq columns are still produced from a Levenshtein traceback at
that position. If the optimal alignment at the position contains an indel,
you'll see a `BULGE`-style CIGAR even in `--distance-mode hamming`.
To force a strictly mismatch-only output, post-filter rows whose CIGAR
contains `I` or `D`, or pass `--no-compute-mismatches` to skip alignment
entirely.

```bash
# Indel-aware (default)
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3

# Faster, mismatch-only search (e.g. for tool comparisons)
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3 \
  --distance-mode hamming
```

---

## PAM Filtering

Stomata does PAM filtering by IUPAC pattern matching, not by named tier.
The PAM is extracted from the genome and matched against `--pam`; the
search itself does not constrain on PAM, so the same search can be
filtered from saved annotated output. Re-running the CLI with another `--pam` performs a new search.

### Pattern syntax

`--pam` accepts any string of A, C, G, T, U and IUPAC ambiguity codes
(R, Y, S, W, K, M, B, D, H, V, N).

| Use case | Flag | What it matches |
|----------|------|-----------------|
| SpCas9 canonical | `--pam NGG` | `AGG`, `CGG`, `GGG`, `TGG` |
| SpCas9 NGG + NAG | `--pam NRG` | NGG ∪ NAG |
| SpCas9 relaxed (NAG only) | `--pam NAG` | NAG |
| Cas12a / Cpf1 | `--pam TTTV --pam-position 5prime` | TTTA, TTTC, TTTG |
| Engineered SpRY | `--pam NRN` | broad |
| No filter (default) | (omit `--pam`) | every position kept |

### PAM position

- `--pam-position 3prime` (default): PAM follows the spacer (Cas9-style).
- `--pam-position 5prime`: PAM precedes the spacer (Cas12a-style).

### Caveat: PAM and `--no-compute-mismatches`

In **Levenshtein** mode, PAM filtering needs the alignment-validity
re-filter that lives inside mismatch annotation. Combining `--pam` with
`--no-compute-mismatches` in Levenshtein mode is rejected with an error.
In **Hamming** mode (no bulges), the early PAM filter is exact and
`--no-compute-mismatches` is allowed.

### PAM as annotation

Even without `--pam`, the `pam_seq` column in TSV output reports the
sequence Stomata extracted at each hit's PAM position. You can post-filter
on that column with awk/pandas if you prefer not to re-run the search.

---

## Strand Control

`--strand` selects which search passes run. Skipping a pass also skips
its work, so non-default modes are roughly twice as fast as `both`.

| Mode | Forward pass | RC pass | Output strands |
|------|--------------|---------|----------------|
| `both` (default) | yes | yes | `+` and `-` |
| `plus` | yes | no | `+` only |
| `minus` | no | yes | `-` only |

```bash
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3 --strand plus
```

---

## Region Restriction

Use `--region` to constrain the search to a single chromosome or window.
Coordinates are 1-based with an exclusive end, matching common genomics
conventions:

```bash
# Whole chromosome
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3 \
  --region chr19

# A specific window
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3 \
  --region chr6:43730000-43740000
```

The region is enforced at the search engine level — positions outside the
window are not scanned.

---

## Spacer Summaries

Two complementary summary modes are available.

### `--summary`: per-distance counts

Replaces per-hit output with aggregated counts by edit distance. Useful
for quick burden assessment when you don't need per-hit details.

```bash
# JSON (default)
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 4 --summary

# TSV
./build/src/stomata --genome hg38.fa.st --spacer-file guides.txt --threshold 4 \
  --summary --summary-format tsv
```

JSON example:

```json
{
  "threshold": 4,
  "total_spacers_processed": 1,
  "spacers_skipped": 0,
  "spacers": [
    {
      "name": "EMX1",
      "sequence": "GAGTCCGAGCAGAAGAAGAA",
      "total_hits": 127,
      "hits_by_distance": {"0": 1, "1": 3, "2": 15, "3": 42, "4": 66}
    }
  ]
}
```

TSV header: `spacer`, `sequence`, `total_hits`, `d0`, `d1`, ..., `d<threshold>`.

### `--spacer-summary`: aggregate CFD (aCFD) report

For a CRISPick-style per-spacer activity summary, write a separate TSV
alongside the main per-hit output:

```bash
./build/src/stomata \
  --genome hg38.fa.st \
  --spacer-file guides.txt \
  --threshold 3 \
  --pam NGG \
  --spacer-summary guide_summary.tsv \
  --output offtargets.tsv
```

The aCFD methodology follows CRISPick 2025 (Doench lab):

- **SDR** = guide positions 4–20 (1-indexed); positions 1–3 are PAM-distal
  and tolerated.
- A hit contributes to **aCFD** only if its SDR mismatches ≤ 1 and the
  alignment has no indels (bulged sites are inactive).
- **aCFD** = sum of per-hit CFD scores over qualifying hits.
- A spacer is flagged **promiscuous** when `aggregate_cfd > --acfd-threshold`
  (default `4.8`).
- aCFD is defined for 20 bp Cas9 spacers only; non-20 bp spacers report
  `.` in the aCFD columns.

#### `--spacer-summary` columns

| Column | Description |
|--------|-------------|
| `spacer` | Spacer name |
| `sequence` | Spacer sequence |
| `n_d0` | Hits at edit distance 0 |
| `n_d1` | Hits at edit distance 1 |
| `n_d2` | Hits at edit distance 2 |
| `n_d3p` | Hits at edit distance ≥ 3 |
| `n_sdr1_sites` | Hits with ≤ 1 mismatch in SDR (aCFD-eligible) |
| `aggregate_cfd` | Sum of CFD scores for SDR-1 hits |
| `is_promiscuous` | `true` if `aggregate_cfd > --acfd-threshold` |
| `n_<label>_overlaps` | (one per `--intersect-bed LABEL:FILE`) hits whose coordinates intersect the named BED |

### BED intersection

`--intersect-bed LABEL:FILE` is repeatable and adds one column per BED
file to the spacer-summary report — handy for counting how many of a
guide's off-targets fall in coding regions, regulatory annotations, etc.

```bash
./build/src/stomata \
  --genome hg38.fa.st --spacer-file guides.txt --threshold 3 --pam NGG \
  --spacer-summary summary.tsv \
  --intersect-bed CDS:annotations/cds.bed \
  --intersect-bed Promoter:annotations/promoters.bed
```

---

## Output Formats

### TSV (default, per-hit)

Single-spacer header (12 columns):

```
chrom  start  end  pattern  distance  strand  aligned_seq  cigar  pam_seq  alignment_ambiguous  n_ambiguous_cells  cfd_score
```

Batch mode prepends a `spacer` column (13 columns total):

```
spacer  chrom  start  end  pattern  distance  strand  aligned_seq  cigar  pam_seq  alignment_ambiguous  n_ambiguous_cells  cfd_score
```

| Column | Type | Description |
|--------|------|-------------|
| `spacer` | string | Spacer name (batch mode only) |
| `chrom` | string | Chromosome name |
| `start` | int | 0-based start position |
| `end` | int | 0-based exclusive end position |
| `pattern` | string | Query spacer sequence |
| `distance` | int | Edit distance (Levenshtein, even in Hamming search mode) |
| `strand` | char | `+` (forward) or `-` (reverse) |
| `aligned_seq` | string | Genomic sequence at the alignment, or `.` |
| `cigar` | string | M/I/D run-length CIGAR, or `.` |
| `pam_seq` | string | Extracted PAM bases (length = `--pam-extract-length`), or `.` |
| `alignment_ambiguous` | bool | `true` if multiple equal-cost alignments exist |
| `n_ambiguous_cells` | int | Number of traceback cells with co-optimal ops |
| `cfd_score` | float | CFD activity score (0.0–1.0); −1 when not computed |

### BED

Single-spacer output is BED6:

```
chrom  start  end  pattern  score  strand
```

Batch output adds the spacer name in column 7:

```
chrom  start  end  pattern  score  strand  spacer
```

The BED `score` is `pattern_length − distance` (higher = better match), so
direct loading into IGV / UCSC orders hits by quality.

### JSON

`--format json` produces a JSON array of hit objects. Single-spacer
fields mirror the TSV columns; batch output adds `spacer` and
`spacer_sequence` per record.

---

## Performance

### Hardware acceleration

Stomata uses the GPU automatically when available (CUDA + compute
capability ≥7.0). Pass `--cpu-only` to force CPU. Patterns longer than
64 bp fall back to multi-word CPU Myers regardless of GPU availability.

### Indices

The `.st` index is memory-mapped on the host. Mapping is cheap, but page-in
and GPU transfer still take time. Use `.st` when reusing a reference.

```bash
./build/src/stomata --index-genome hg38.fa     # produces hg38.fa.st (one-time)
./build/src/stomata --genome hg38.fa.st ...    # all later runs
```

### Reproducible benchmarks

Historical timing tables were removed because their raw measurements are not
available in this checkout. Use the [benchmark runner](../validation/benchmarks/README.md)
to record versions, hardware, hashes, repeated timings and failures.

### Tips

- **Use `.st` indices** for repeated searches.
- **Batch spacers** with `--spacer-file` — the genome is loaded once.
- **Filter PAM early** with `--pam` to reduce output volume.
- **Skip annotations** with `--no-compute-mismatches` if you only need
  hit counts (Hamming mode) or coordinates.
- **Cap output** with `--max-hits` and `--max-total-hits` for promiscuous
  spacers or large screens.
- **Skip CFD** with `--no-scores` when you don't need activity scoring.

### Memory

Historical planning estimates below are not current build measurements.
Peak memory also depends on hit density and batch size.

| Genome | CPU memory | GPU VRAM |
|--------|------------|----------|
| E. coli (4.6 Mb) | <100 MB | <50 MB |
| Human (3.0 Gb) | ~5 GB | ~2 GB |
| Wheat (16 Gb) | ~25 GB | ~10 GB |

---

## Behavior worth knowing

- **Default threshold is 3.** Pass `--threshold N` explicitly if you care.
- **Default PAM filter is none.** Without `--pam`, every sequence match
  is reported regardless of PAM context. PAM is annotation, not a search
  constraint — filter saved annotated output to avoid another search. Re-running the CLI performs a new search.
- **PAM-agnostic output is large.** On hg38 at threshold 3, 225 spacers
  produces ~4M hits with no PAM filter vs ~280K with `--pam NGG`. Use
  `--summary` to reduce output volume. Hit caps truncate results after search;
  they do not guarantee lower peak memory and cannot be combined with summaries.
- **U → T normalization is on by default.** RNA spacers work as-is.
  Disable with `--no-treat-u-as-t` for strict validation.
- **Pattern N matches any base at zero cost; genome N is masked
  (cost 1 for any base).** This avoids spurious matches at low-quality
  genome regions.
- **Halo deduplication is on by default.** Myers' semi-global alignment
  reports a distance at every end-position, producing clusters of
  near-duplicate hits around each real match. Stomata clusters consecutive
  same-strand/same-chromosome positions and keeps the minimum-distance
  representative. Pass `--no-deduplicate` for cross-tool validation that
  needs every raw end-position.
- **CFD is undefined for indels.** Hits with bulged alignments receive
  `cfd_score = 0` and are excluded from aCFD aggregation.

---

## Troubleshooting

### Build / install

**`STOMATA_VERSION` undefined** or generic CMake errors
Verify CUDA, spdlog, and zlib are installed in the active environment.
Conda's `compilers gxx=11 spdlog catch2` package set is the supported
configuration.

**`stomata --quickstart` fails** with hit-count mismatches
Most often a CUDA SASS/architecture mismatch. The default fat-binary
covers `sm_70;75;80;86;89;90`; for older or newer GPUs rebuild with
`cmake -DCMAKE_CUDA_ARCHITECTURES=<your_sm> ...`. Then retry.

### Runtime errors

**`Cannot open genome file: ...`**
Verify the path; use absolute paths if in doubt. Stomata accepts
`.fa`, `.fa.gz`, and `.st`.

**`Invalid character in pattern: 'X'`**
Spacers may contain only A/C/G/T/U/N (case-insensitive). U is normalized
unless `--no-treat-u-as-t` is set.

**`--pam requires --compute-mismatches in Levenshtein mode`**
In Levenshtein mode the PAM filter relies on the alignment-validity
re-filter inside mismatch annotation. Either keep mismatch computation
on (default) or switch to `--distance-mode hamming`.

**`Warning: CUDA not available, falling back to CPU`**
GPU was not detected. Verify `nvidia-smi` works, that CUDA matches the
toolkit you built against, and that the GPU's compute capability is in
the SASS targets the binary was built with.

**`std::bad_alloc` / out of memory**
Lower `--threshold`, add `--pam`, or split into smaller batches or regions.
Hit caps and summary formatting are not peak-memory limits.

### Performance issues

**Slow first load**
Build the `.st` index once (`--index-genome`); subsequent runs mmap it.

**GPU slower than CPU on small queries**
GPU launch overhead dominates for tiny inputs. Expected; CPU is fine for
single-spacer / small genome runs.

---

## Examples

### 1. Single guide, SpCas9, capture both NGG and NAG

```bash
./build/src/stomata \
  --genome hg38.fa.st \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam NRG \
  --output vegfa_offtargets.tsv
```

`NRG` (R = A or G) covers NAG ∪ NGG. PAM type is preserved in the
`pam_seq` column so you can split downstream.

### 2. Restrict to one locus

```bash
./build/src/stomata \
  --genome hg38.fa.st \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 4 \
  --region chr6:43730000-43740000
```

### 3. Batch + spacer summary + BED intersection

```bash
./build/src/stomata \
  --genome hg38.fa.st \
  --spacer-file my_library.txt \
  --threshold 3 \
  --pam NGG \
  --output library_offtargets.tsv \
  --spacer-summary library_summary.tsv \
  --intersect-bed CDS:annotations/cds.bed
```

`library_summary.tsv` reports aCFD and flags promiscuous
candidates; sort it explicitly to rank guides; the new `n_CDS_overlaps` column counts off-targets that fall
in coding regions per guide.

### 4. Cas12a (5'-PAM, TTTV)

```bash
./build/src/stomata \
  --genome hg38.fa.st \
  --pattern AAACGAGCAGAAGAAGAAGGT \
  --threshold 3 \
  --pam TTTV --pam-position 5prime \
  --pam-extract-length 4
```

### 5. Strand-specific comparison

```bash
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3 \
  --strand plus  --output fwd.tsv
./build/src/stomata --genome hg38.fa.st --pattern SEQ --threshold 3 \
  --strand minus --output rev.tsv
```

### 6. Downstream pandas

```python
import pandas as pd

df = pd.read_csv("offtargets.tsv", sep="\t")

# Hits with no indels (cigar has no I or D)
no_indels = df[~df["cigar"].str.contains("[ID]", na=False)]

# Top SpCas9 NGG hits by CFD
ngg = df[df["pam_seq"].str.match(r".GG$", na=False)]
ngg.sort_values("cfd_score", ascending=False).head(20).to_csv(
    "top_cfd.tsv", sep="\t", index=False)
```

---

## Reproducibility

The search scans the supplied reference within the selected distance threshold.
Levenshtein output is halo-deduplicated by default; use `--no-deduplicate`
when every raw end-position is required. Hit caps deliberately truncate output.
Record the version, reference/guide hashes, flags and CPU/GPU path.
See [testing](INSTALL.md#testing) and [validation](../validation/README.md) for
reproducible checks and their limits. Candidate sequence matches are not
experimental measurements of cleavage.

---

## Citation

```bibtex
@software{Stomata2026,
  author  = {Simpson, Danny and Sanjana, Neville E. and Lappalainen, Tuuli},
  title   = {Stomata: GPU-accelerated exhaustive CRISPR off-target search},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/simpsondl/stomata}
}
```

## License

MIT — see [LICENSE](../LICENSE).

## Getting help

- Issues: <https://github.com/simpsondl/stomata/issues>
- `stomata --help` for the full flag reference.

---

**Version:** 1.0.0
**Authors:** Danny Simpson, Neville E. Sanjana, Tuuli Lappalainen (New York Genome Center)
