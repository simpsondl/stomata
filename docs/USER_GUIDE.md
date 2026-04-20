# Cystidia User Guide

Complete guide to using Cystidia for exhaustive CRISPR off-target search.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Basic Usage](#basic-usage)
4. [Advanced Features](#advanced-features)
5. [Output Formats](#output-formats)
6. [Performance Optimization](#performance-optimization)
7. [Troubleshooting](#troubleshooting)

---

## Installation

### System Requirements

- **Operating System:** Linux (tested on Ubuntu 20.04+)
- **CPU:** x86_64 with AVX2 support
- **GPU:** NVIDIA GPU with CUDA compute capability ≥7.0 (optional but recommended)
- **Memory:** 8GB RAM minimum (16GB+ recommended for human genome)
- **Storage:** 10GB for software + genome indices

### Dependencies

#### Required
- CMake ≥3.20
- C++17 compiler (GCC 11+ or Clang 14+)
- CUDA Toolkit 11.8+ (for GPU acceleration)
- Python 3.10+

#### Optional (for benchmarking/analysis)
- jinja2, plotly (for HTML reports)
- pandas, matplotlib (for analysis)

### Installation Steps

#### Option 1: Conda (Recommended)

```bash
# Clone repository
git clone https://github.com/simpsondl/cystidia.git
cd cystidia

# Create and activate conda environment
conda create -n cystidia python=3.10 -y
conda activate cystidia

# Install build dependencies
conda install -c conda-forge cmake compilers gxx=11 spdlog catch2 -y

# Install Python packages
pip install pandas numpy matplotlib jinja2 plotly

# Build Cystidia
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Verify installation
./build/src/cystidia --version
```

#### Option 2: System Dependencies

```bash
# Install system packages (Ubuntu/Debian)
sudo apt update
sudo apt install -y build-essential cmake g++ libspdlog-dev catch2 python3-pip

# Install CUDA Toolkit (if not present)
# Follow: https://developer.nvidia.com/cuda-downloads

# Build and install
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
sudo cmake --install build
```

---

## Quick Start

### 30-Second Demo

Search for off-targets of a single spacer sequence in a small genome:

```bash
# Index a genome (one-time)
./build/src/cystidia --index-genome tests/data/synthetic/genome.fa

# Search for off-targets
./build/src/cystidia \
  --genome tests/data/synthetic/genome.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --format tsv
```

### Real-World Example: Human Genome

```bash
# Step 1: Download and index hg38 genome
./scripts/prepare_hg38_canonical.sh

# Step 2: Search single spacer (CPU mode)
./build/src/cystidia \
  --genome data/genomes/hg38_canonical.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam-filter both \
  --format tsv \
  --output vegfa_offtargets.tsv

# Note: Using `--pam-filter both` (not just `ngg`) ensures comprehensive
# SpCas9 off-target search including both NGG and NAG PAMs

# Step 3: Search multiple spacers (GPU batch mode, much faster)
./build/src/cystidia \
  --genome data/genomes/hg38_canonical.cy \
  --spacer-file my_guides.txt \
  --threshold 3 \
  --pam-filter both \
  --format tsv \
  --output batch_results.tsv
```

---

## Basic Usage

### Command-Line Interface

```
cystidia [OPTIONS] --genome GENOME --pattern PATTERN --threshold N
cystidia [OPTIONS] --genome GENOME --spacer-file FILE --threshold N
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--genome FILE` | Genome file (.fa, .fa.gz, or .cy) | `hg38.fa.cy` |
| `--threshold N` | Maximum edit distance (0-255) | `3` |

**One of:**
| Argument | Description | Example |
|----------|-------------|---------|
| `--pattern SEQ` | Single spacer sequence (20bp) | `GAGTCCGAGCAGAAGAAGAA` |
| `--spacer-file FILE` | File with multiple spacers | `guides.txt` |

### Common Options

| Option | Description | Default |
|--------|-------------|---------|
| `--format FMT` | Output format: `tsv`, `bed` | `tsv` |
| `--output FILE` | Output file (or stdout) | stdout |
| `--pam-filter FILTER` | PAM filter: `none`, `ngg`, `nag`, `both`, `any` | `none` |
| `--strand MODE` | Strand selection: `both`, `plus`, `minus` | `both` |
| `--max-hits N` | Maximum hits per spacer | unlimited |
| `--no-mismatches` | Skip alignment annotation (faster) | compute |
| `--no-scores` | Disable CFD activity scoring (faster) | compute |
| `--version` | Show version and build info | |
| `--help` | Show full help message | |

### Spacer File Format

Cystidia supports flexible spacer file formats:

```
# Format 1: Name + Sequence (tab-separated)
VEGFA_site1	GAGTCCGAGCAGAAGAAGAA
EMX1_site2	GGAATCCCTTCTGCAGCACC
FANCF_site3	GAAGATGGACCTGATCGACA

# Format 2: Sequence only (one per line)
GAGTCCGAGCAGAAGAAGAA
GGAATCCCTTCTGCAGCACC
GAAGATGGACCTGATCGACA

# Comments and blank lines are ignored
# This is a comment
```

---

## Advanced Features

### Genome Indexing

For repeated searches, create a binary `.cy` index for instant loading:

```bash
# Create index (one-time, ~30 seconds for hg38)
./build/src/cystidia --index-genome hg38.fa

# This creates: hg38.fa.cy

# Subsequent searches use mmap (instant loading)
./build/src/cystidia --genome hg38.fa.cy --pattern ACGT... --threshold 3
```

**Benefits:**
- **100-1000x faster loading** (3s → 0.003s for hg38)
- **Memory efficient** (OS pages in only accessed regions)
- **Persistent** (reusable across runs)

### Strand-Specific Search

By default, Cystidia runs both search passes (forward and reverse-complement) and
reports hits on both strands. `--strand` selects which pass(es) to run; each
non-default mode skips one pass entirely, so it's also faster:

- `both` (default): run both passes; report + and - hits.
- `plus`: run only the forward pass; output contains only + hits.
- `minus`: run only the RC pass; output contains only - hits.

```bash
# Plus-strand hits only (skips the RC search)
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --strand plus

# Minus-strand hits only (skips the forward search)
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --strand minus
```

### Alignment Modes: Mismatch-Only vs Edit Distance

Cystidia supports two alignment modes that affect how genomic sites are scored and aligned:

| Mode | Algorithm | Penalizes | Use Case |
|------|-----------|-----------|----------|
| Default | Myers' bit-vector (edit distance) | Mismatches + indels | **Exhaustive biological search** |
| `--shift-and` | Shift-and (mismatch-only) | Mismatches only | Direct comparison with Cas-OFFinder |

**Why does this matter?**

For the same genomic site, different alignment algorithms can produce:
- **Different edit distances:** 2 edits with bulge vs 3 mismatches
- **Different PAM sequences:** NAG vs NGG due to shifted alignment boundary
- **Different filtering outcomes:** Site kept vs filtered when using PAM filters

**Example:**

```
Spacer:             GCTCGGGGACACAGGATCCC
Genome (w/ indel):  GCTCG-GGACACAGGATCCC    ← Edit distance alignment (1 deletion)
Genome (mismatch):  GCTCGaGGACACAGGATCCC    ← Mismatch-only alignment (1 substitution)
PAM (w/ indel):     TAG (NAG type)
PAM (mismatch):     AGG (NGG type)
```

With `--threshold 2 --pam-filter ngg`:
- Edit distance mode finds this at distance 1, but PAM is TAG → **filtered out**
- Shift-and mode finds this at distance 1, PAM is AGG → **kept**

#### When to Use Each Mode

```bash
# DEFAULT (edit distance): For complete biological search
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam-filter both

# SHIFT-AND: For comparing with Cas-OFFinder results
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam-filter ngg \
  --shift-and
```

**Recommendations:**

- **For comparing with Cas-OFFinder or other mismatch-only tools:** Use `--shift-and`
- **For complete biological off-target search:** Use default edit distance mode
- **For SpCas9 specificity screening:** Use default mode + `--pam-filter both` (NGG or NAG)
- **For maximum exhaustiveness:** Use default mode + `--pam-filter none`

⚠️ **Important PAM Filtering Interaction:**

When using PAM filters, the alignment mode affects which sites pass through:
- **Edit distance mode** may find alignments with bulges that shift the PAM boundary
- Some sites reported by Cas-OFFinder as "3 mismatches + NGG" might appear in Cystidia as "2 edits with bulge + NAG"
- Using `--pam-filter both` (instead of just `ngg`) captures these alignment variants

**See also:** [Cross-Tool Validation](CROSS_TOOL_VALIDATION.md#investigation-5-pam-filtering-and-alignment-mode-interactions) for detailed analysis of alignment differences.

### PAM Filtering

Cystidia can filter results by PAM sequence compatibility:

| Filter | Description | Use Case |
|--------|-------------|----------|
| `none` | No PAM filtering | **Complete exhaustive search** |
| `ngg` | Only NGG PAMs | SpCas9 strict (use with `--shift-and` for Cas-OFFinder compatibility) |
| `nag` | Only NAG PAMs | SpCas9 relaxed |
| `both` | NGG or NAG PAMs | **SpCas9 standard (recommended for most use)** |
| `any` | Any PAM (including incomplete) | Other nucleases |

```bash
# Default: no PAM filtering (all sites reported)
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3

# Recommended for SpCas9: Filter for NGG/NAG PAMs
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam-filter both    # ← Recommended: captures all plausible SpCas9 sites
```

⚠️ **PAM Filtering and Alignment Modes:**

PAM filtering interacts with the alignment mode (see [Alignment Modes](#alignment-modes-mismatch-only-vs-edit-distance)):
- **Edit distance mode** (default) may produce different PAMs than Cas-OFFinder due to indel alignments
- **Shift-and mode** (`--shift-and`) produces the same PAMs as Cas-OFFinder
- Using `--pam-filter both` ensures you don't miss sites where indel alignments shift NGG→NAG

**Examples:**

```bash
# For exhaustive SpCas9 search (recommended)
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam-filter both    # Includes both NGG and NAG

# For Cas-OFFinder-compatible search
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --pam-filter ngg \
  --shift-and          # Mismatch-only alignment like Cas-OFFinder
```

### Output Limiting

For promiscuous spacers, limit output size:

```bash
# Return only top 1000 hits per spacer
./build/src/cystidia --genome hg38.fa.cy \
  --spacer-file guides.txt \
  --threshold 4 \
  --max-hits 1000
```

### CFD Activity Scoring

Cystidia computes CFD (Cutting Frequency Determination) scores for each off-target hit, predicting the likelihood of cleavage activity. CFD scores are based on Doench et al. 2016 and are the industry standard for off-target activity prediction.

**Understanding CFD Scores:**

| CFD Score | Risk Tier | Interpretation |
|-----------|-----------|----------------|
| ≥ 0.1 | HIGH | Significant off-target potential, may cause cleavage |
| 0.01 - 0.1 | MEDIUM | Moderate risk, consider experimental validation |
| < 0.01 | LOW | Unlikely to cause off-target cleavage |

**Key points:**
- Perfect match with NGG PAM = CFD score of 1.0
- Each mismatch multiplies a position-specific penalty (0.0-1.0)
- PAM type affects score: NGG = 1.0, NAG ≈ 0.26, other PAMs ≈ 0
- **Bulges (insertions/deletions) are not scored** - CFD is undefined for indels (score = 0)

**Example:** An off-target with 2 mismatches at positions 7 and 10 with NGG PAM:
```
CFD = penalty(pos7) × penalty(pos10) × PAM_penalty(NGG)
    = 0.57 × 0.88 × 1.0 = 0.50 (HIGH risk)
```

**Disabling scoring for speed:**
```bash
# Skip CFD scoring for faster output
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 4 \
  --no-scores
```

### MIT Specificity Scoring

In batch mode, Cystidia can compute MIT Specificity Scores to assess overall guide RNA specificity by aggregating all off-targets for each spacer. This metric, developed by the Hsu/Zhang lab at MIT, is widely used for guide selection.

**Understanding MIT Specificity:**

| MIT Score | Specificity Tier | Interpretation |
|-----------|-----------------|----------------|
| ≥ 80 | EXCELLENT | Very few predicted off-targets, ideal guide |
| 50 - 80 | GOOD | Acceptable specificity for most applications |
| 20 - 50 | FAIR | Moderate off-target burden, consider alternatives |
| < 20 | POOR | High off-target burden, avoid if possible |

**Formula:**
```
MIT Specificity = 100 / (100 + Σ(CFD_scores))
```

where the sum is over all off-target CFD scores.

**Key points:**
- Score of 100 = perfect specificity (no predicted off-targets)
- Lower scores indicate higher off-target burden
- Commonly used threshold: >50 for acceptable guides
- **Only available in batch mode** (requires `--spacer-file`)
- Requires CFD scoring to be enabled (default)

**Usage:**
```bash
# Compute MIT specificity for multiple guides
./build/src/cystidia --genome hg38.fa.cy \
  --spacer-file my_guides.txt \
  --threshold 4 \
  --mit-score \
  --output results.tsv
```

The MIT score and tier appear in the last two columns of TSV output for batch searches, allowing easy ranking and filtering of guide candidates.

**Example output interpretation:**
```
spacer          mit_specificity  specificity_tier
EMX1_guide1     87.3             EXCELLENT       # Great choice
EMX1_guide2     62.5             GOOD            # Acceptable
EMX1_guide3     42.1             FAIR            # Consider alternatives
EMX1_guide4     15.8             POOR            # Avoid
```

### RNA Input Support (U→T Normalization)

Cystidia accepts spacers containing uracil (`U`) and automatically converts them to thymine (`T`). This is enabled by default:

```bash
# RNA spacer input (U is treated as T)
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGUCCGAGCAGAAGAAGAA \
  --threshold 3

# Equivalent to searching with GAGTCCGAGCAGAAGAAGAA

# Disable normalization if needed
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --no-treat-u-as-t
```

**Notes:**
- Both uppercase `U` and lowercase `u` are normalized
- Works with both `--pattern` and `--spacer-file` input
- Disabled via `--no-treat-u-as-t` (pattern containing U will be rejected)

### Summary Mode (Aggregated Counts)

For quick assessment of off-target burden without per-hit details, use `--summary`:

```bash
# JSON summary (default)
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 4 \
  --summary

# TSV summary
./build/src/cystidia --genome hg38.fa.cy \
  --spacer-file guides.txt \
  --threshold 4 \
  --summary --summary-format tsv
```

**JSON output example:**
```json
{
  "threshold": 4,
  "total_spacers_processed": 1,
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

**TSV output columns:**
| Column | Description |
|--------|-------------|
| `spacer` | Spacer name |
| `sequence` | Spacer sequence |
| `total_hits` | Total hits within threshold |
| `d0, d1, d2, ...` | Hit counts at each distance |

### Performance Tuning

For maximum speed on large batches:

```bash
# GPU batch mode (up to 800x faster for 100+ spacers)
./build/src/cystidia --genome hg38.fa.cy \
  --spacer-file 100_guides.txt \
  --threshold 3 \
  --pam-filter both \
  --no-compute-mismatches  # Skip detailed annotations for speed
```

---

## Output Formats

### TSV Format (Default)

Tab-separated format with 23 biological annotation columns (including CFD and MIT scores):

```tsv
spacer  chrom   start   end     pattern distance        strand  aligned_seq     mismatch_pos    edit_types      pam_seq pam_type        seed_edits      distal_edits    seed_mismatches seed_dna_bulges seed_rna_bulges distal_mismatches       distal_dna_bulges       distal_rna_bulges       alignment_ambiguous     cfd_score       risk_tier       mit_specificity specificity_tier
VEGFA_site3     chr6    43737006        43737026        GAGTCCGAGCAGAAGAAGAA    2       +       GAGTCCAAGCAGAAGAAGAA    7       MISMATCH        AGG     NGG     1       1       1       0       0       1       0       0       false   0.486   HIGH    87.3    EXCELLENT
```

#### Column Descriptions

| Column | Type | Description |
|--------|------|-------------|
| `spacer` | string | Spacer name (or sequence if unnamed) |
| `chrom` | string | Chromosome name |
| `start` | int | 0-based start position |
| `end` | int | 0-based exclusive end position |
| `pattern` | string | Query spacer sequence |
| `distance` | int | Edit distance (Levenshtein) |
| `strand` | char | '+' (forward) or '-' (reverse) |
| `aligned_seq` | string | Aligned genomic sequence |
| `mismatch_pos` | string | 1-based positions with edits (comma-separated) |
| `edit_types` | string | Edit type at each position (MISMATCH, DNA_BULGE, RNA_BULGE) |
| `pam_seq` | string | 3bp PAM sequence downstream |
| `pam_type` | string | NGG, NAG, OTHER, or INCOMPLETE |
| `seed_edits` | int | Total edits in seed region (positions 1-10) |
| `distal_edits` | int | Total edits in distal region (positions 11-20) |
| `seed_mismatches` | int | Mismatches in seed |
| `seed_dna_bulges` | int | DNA bulges in seed |
| `seed_rna_bulges` | int | RNA bulges in seed |
| `distal_mismatches` | int | Mismatches in distal |
| `distal_dna_bulges` | int | DNA bulges in distal |
| `distal_rna_bulges` | int | RNA bulges in distal |
| `alignment_ambiguous` | bool | Multiple optimal alignments exist |
| `cfd_score` | float | CFD activity score (0.0-1.0), higher = more likely to cleave |
| `risk_tier` | string | Risk classification: LOW (<0.01), MEDIUM (0.01-0.1), HIGH (≥0.1) |
| `mit_specificity` | float | MIT specificity score (0-100), only in batch mode with `--mit-score` |
| `specificity_tier` | string | Specificity tier: POOR (<20), FAIR (20-50), GOOD (50-80), EXCELLENT (≥80) |

### BED Format

Minimal 6-column BED for genome browsers:

```bash
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --format bed \
  --output offtargets.bed
```

Output:
```bed
chr6    43737006        43737026        GAGTCCGAGCAGAAGAAGAA    2       +
chr12   25398281        25398301        GAGTCCGAGCAGAAGAAGAA    3       -
```

Load into IGV, UCSC Genome Browser, or other tools.

---

## Performance Optimization

### Hardware Acceleration

Cystidia automatically detects and uses GPUs when available:

| Mode | When Used | Speedup |
|------|-----------|---------|
| CPU | Single spacer or no GPU | 1x baseline |
| GPU | Batch mode (2+ spacers) + GPU detected | **10-800x faster** |

**GPU Requirements:**
- NVIDIA GPU with CUDA compute capability ≥7.0
- CUDA Toolkit 11.8+ installed
- Sufficient VRAM (2GB+ for hg38)

### Performance Guidelines

| Task | Recommended Mode | Expected Time (hg38) |
|------|------------------|----------------------|
| 1 spacer | CPU or GPU | 0.5-70 seconds |
| 10 spacers | GPU batch | 1-2 seconds |
| 100 spacers | GPU batch | 8-10 seconds |
| 1000 spacers | GPU batch | 20-25 seconds |

### Optimization Tips

1. **Use .cy indices** for repeated searches (100-1000x faster loading)
2. **Batch spacers** when searching multiple guides (10-800x faster)
3. **Filter PAMs early** with `--pam-filter both` for SpCas9 (2-5x fewer hits)
4. **Skip annotations** with `--no-mismatches` when only counting hits
5. **Limit output** with `--max-hits` for promiscuous spacers

### Memory Requirements

| Genome | CPU Memory | GPU VRAM |
|--------|------------|----------|
| E. coli (4.6 Mb) | <100 MB | <50 MB |
| Human (3.0 Gb) | ~5 GB | ~2 GB |
| Wheat (16 Gb) | ~25 GB | ~10 GB |

---

## Troubleshooting

### Common Issues

#### "Cannot open genome file"
```
Error: Cannot open genome file: hg38.fa
```
**Solution:** Verify file path and permissions. Use absolute paths if needed.

#### "Invalid character in pattern"
```
Warning: Skipping spacer 'guide_5': Invalid character in pattern: 'S'
```
**Solution:** Spacers must contain only A, C, G, T, or N. Check input file for typos or non-IUPAC characters.

#### "CUDA error" or GPU not detected
```
Warning: CUDA not available, falling back to CPU
```
**Solution:** 
- Verify NVIDIA driver: `nvidia-smi`
- Check CUDA installation: `nvcc --version`
- Ensure GPU compute capability ≥7.0

#### Out of memory
```
terminate called after throwing an instance of 'std::bad_alloc'
```
**Solution:**
- Use smaller threshold (reduces hits)
- Filter by PAM (`--pam-filter both`)
- Limit output (`--max-hits 10000`)
- Use machine with more RAM

### Performance Issues

#### Slow genome loading
**Solution:** Create `.cy` index with `--index-genome` (one-time cost, massive speedup)

#### GPU slower than expected
**Possible causes:**
- Small batch size (< 10 spacers): GPU overhead dominates
- Old GPU (compute capability < 7.5): Limited parallelism
- System bottleneck: Check `nvidia-smi` for GPU utilization

**Solution:** Increase batch size or use CPU for small queries.

### Getting Help

- **Documentation:** [docs/](../docs/)
- **Issues:** [github.com/simpsondl/cystidia/issues](https://github.com/simpsondl/cystidia/issues)
- **Contact:** dsimpson@example.edu

---

## Examples Gallery

### Example 1: Basic Off-Target Search

Find all sites within 3 mismatches for VEGFA spacer:

```bash
./build/src/cystidia \
  --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --output vegfa_offtargets.tsv
```

### Example 2: High-Specificity Screen

Find only NGG PAM sites with ≤2 seed disruptions:

```bash
./build/src/cystidia \
  --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 2 \
  --pam-filter ngg \
  --output high_specificity.tsv
```

Filter output with awk:
```bash
# Only sites with no seed edits
awk -F'\t' '$13 == 0' high_specificity.tsv
```

### Example 3: Batch Processing

Search 100 spacers in parallel:

```bash
# Create spacer file
cat > my_library.txt << EOF
EMX1_site1	GAGTCCGAGCAGAAGAAGAA
VEGFA_site2	GGAATCCCTTCTGCAGCACC
FANCF_site3	GAAGATGGACCTGATCGACA
EOF

# Run batch search (GPU accelerated)
./build/src/cystidia \
  --genome hg38.fa.cy \
  --spacer-file my_library.txt \
  --threshold 3 \
  --pam-filter both \
  --max-hits 5000 \
  --output library_offtargets.tsv
```

### Example 4: Strand-Specific Analysis

Compare forward vs reverse strand targeting:

```bash
# Forward strand only
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 --strand plus \
  --output fwd_hits.tsv

# Reverse strand only
./build/src/cystidia --genome hg38.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 --strand minus \
  --output rev_hits.tsv
```

### Example 5: Downstream Analysis (Python)

```python
import pandas as pd

# Load results
df = pd.read_csv('offtargets.tsv', sep='\t')

# Filter high-risk sites (seed disruptions)
high_risk = df[df['seed_edits'] <= 1]
print(f"High-risk sites: {len(high_risk)}")

# Count by PAM type
pam_dist = df['pam_type'].value_counts()
print(pam_dist)

# Export high-risk sites for validation
high_risk.to_csv('validate_these.csv', index=False)
```

---

## Advanced Topics

### Choosing the Right Settings: Decision Tree

Different use cases require different combinations of alignment mode and PAM filtering. Use this guide to select optimal settings:

```
What's your goal?

├─ Comparing results with Cas-OFFinder or other mismatch-only tools?
│  └─→ Use: --shift-and --pam-filter ngg
│     Example: ./build/src/cystidia --genome hg38.fa.cy --pattern SEQ \
│              --threshold 3 --shift-and --pam-filter ngg
│
├─ Complete off-target search for SpCas9 (NGG/NAG PAMs)?
│  └─→ Use: (default mode) --pam-filter both
│     Example: ./build/src/cystidia --genome hg38.fa.cy --pattern SEQ \
│              --threshold 3 --pam-filter both
│     ✓ Finds all plausible SpCas9 sites (NGG and NAG PAMs)
│     ✓ Includes sites where indel alignments shift PAM boundaries
│     ✓ Recommended for most SpCas9 applications
│
├─ Maximum exhaustiveness: find ALL genomic sites within threshold?
│  └─→ Use: (default mode) --pam-filter none
│     Example: ./build/src/cystidia --genome hg38.fa.cy --pattern SEQ \
│              --threshold 3 --pam-filter none
│     ✓ No PAM restrictions
│     ✓ Good for research or custom PAM requirements
│     ✓ Post-process results by your own PAM requirements
│
├─ Using a non-SpCas9 nuclease (different PAM)?
│  └─→ Use: (default mode) --pam-filter none + post-processing
│     Example: ./build/src/cystidia --genome hg38.fa.cy --pattern SEQ \
│              --threshold 3 --pam-filter none --output results.tsv
│     Then filter results by your specific PAM in downstream analysis
│
└─ Validating Cystidia against other tools?
   └─→ Run both modes and check for alignment differences
       Example: See scripts/run_cross_tool_validation.py --cystidia-mode both
```

**Quick Reference:**

| Use Case | Mode | PAM Filter | Command |
|----------|------|------------|---------|
| SpCas9 screening (recommended) | default | `both` | `--pam-filter both` |
| Cas-OFFinder comparison | `--shift-and` | `ngg` | `--shift-and --pam-filter ngg` |
| Maximum exhaustiveness | default | `none` | `--pam-filter none` |
| Other nucleases | default | `none` | `--pam-filter none` + custom filter |

### N-Handling Semantics

Cystidia uses biologically-correct N-handling:

- **Pattern N (wildcard):** Matches A, C, G, or T with distance 0
- **Genome N (masked):** Matches nothing (distance 1 for all)

This prevents false matches at low-quality genome regions.

### Edit Distance vs Mismatches

Cystidia reports **edit distance** (Levenshtein distance), which includes:
- Mismatches (substitutions)
- DNA bulges (insertions in genome)
- RNA bulges (deletions in genome / insertions in guide)

For traditional "3-mismatch" searches, filter on `total_mismatches` column.

### Alignment Ambiguity

When `alignment_ambiguous` is `true`, multiple optimal alignments exist. This typically occurs with:
- Repetitive sequences
- Homopolymer runs
- Multiple edit patterns yielding same distance

The reported alignment is arbitrary but valid. Consider these sites carefully.

### Reproducibility

Cystidia guarantees:
- **100% sensitivity:** Never misses a site within threshold
- **Bit-identical results:** Same input always produces same output
- **No heuristics:** Exhaustive search, no approximations

---

## Changelog

### v0.4.0 (February 2026) - New Features

**Single-Strand Search (`--strand plus` / `--strand minus`):**
Each non-default `--strand` mode skips one search pass (forward or reverse-complement),
reporting only hits on the selected strand. Useful when only one orientation is
biologically relevant, and reduces runtime by skipping one pass.

**RNA Input Support (`--treat-u-as-t`, default: on):**
Spacer sequences containing uracil (U) are automatically normalized to thymine (T) before
searching. This enables direct use of RNA-centric spacer sequences without manual
preprocessing. Can be disabled with `--no-treat-u-as-t`.

**Summary Mode (`--summary`):**
Outputs aggregated hit counts by edit distance instead of full per-hit details. Supports
both JSON (default) and TSV output (`--summary-format tsv`). Compatible with all other
filters (PAM, strand). Useful for quick specificity assessment and greatly
reduced output size for large screens.

**Bug Fixes:**
- Fixed ToolAdapter TSV parser to handle current 22-column output format
- Fixed ChangeseqValidator to parse current Cystidia TSV format
- Fixed CLI version test to match actual version string

### v0.3.1 (February 2026) - Critical Bug Fixes

**Hit Deduplication (Bug Fix):**
Myers' semi-global alignment reports edit distance at every genome position, creating
"halos" of nearby low-distance positions around each real match. This caused hit counts
to be inflated by approximately 8x compared to other tools. A new deduplication step
clusters consecutive positions on the same strand/chromosome and keeps only the
minimum-distance hit per cluster.

**Minus Strand PAM Extraction (Bug Fix):**
PAM sequences were incorrectly extracted for minus-strand hits. The PAM was being taken
from immediately after the alignment end position on the forward strand, regardless of
strand. For minus-strand hits, the PAM should be extracted from before the alignment
start on the forward strand, then reverse-complemented. This caused approximately 50%
of valid hits to be incorrectly filtered when using `--pam-filter` options.

**N-Handling Consistency (Bug Fix):**
The DP alignment used for mismatch annotation treated genome N as a wildcard (cost 0),
while the Myers search algorithm treats genome N as masked (cost 1). This inconsistency
has been resolved so both stages use the same semantics: pattern N = wildcard (cost 0),
genome N = masked (cost 1).

**Aligned Sequence Extraction (Bug Fix):**
For hits with indels (DNA/RNA bulges), the aligned genomic sequence was extracted using
a fixed-length substring rather than deriving it from the actual DP alignment traceback.
This has been corrected so the `aligned_seq` column accurately reflects the alignment.

**Minus Strand Aligned Sequence (Bug Fix):**
After reverse-complementing the genome sequence for minus-strand comparison, the
alignment region is at the beginning of the string, not the end. The substring
extraction now correctly handles both strands.

### v0.3.0 (February 2026)

- MIT Specificity Score computation
- CHANGE-seq validation framework
- Comprehensive documentation

### v0.2.0 (February 2026)

- Initial public release
- GPU-accelerated Myers' algorithm
- Batch processing, PAM filtering, CFD scoring

---

## Citation

If you use Cystidia in your research, please cite:

```
[Citation to be added upon publication]
```

## License

[License information to be added]

---

**Version:** 0.3.1
**Last Updated:** February 2026
**Authors:** Danny Simpson and contributors
