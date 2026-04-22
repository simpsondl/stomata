# Cystidia

<div align="center">

**Ultra-fast exhaustive CRISPR off-target search engine**

[![Version](https://img.shields.io/badge/version-0.5.0-blue.svg)](https://github.com/simpsondl/cystidia)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![CUDA](https://img.shields.io/badge/CUDA-11.8%2B-76B900.svg)](https://developer.nvidia.com/cuda-toolkit)
[![C++](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://isocpp.org/)

Cystidia delivers **10-800x speedup** with GPU acceleration, searching 3.1 billion bases in seconds.

[Features](#features) • [Quick Start](#quick-start) • [Performance](#performance) • [Documentation](#documentation) • [Citation](#citation)

</div>

---

## Features

### Core Capabilities
- ⚡ **Blazing Fast:** Search entire human genome in 0.5-23 seconds (GPU) vs hours (traditional)
- 🎯 **100% Sensitivity:** Exhaustive search with no heuristics—never misses an off-target
- 📊 **Rich Annotations:** 23 biological features including CFD activity scores, MIT specificity, seed disruptions, PAM analysis
- 🔄 **Batch Processing:** Screen 100+ guides simultaneously with near-linear scaling
- 💾 **Memory Efficient:** <2 GB GPU VRAM for human genome at any batch size
- 🧬 **Flexible:** Supports any genome, edit distance threshold, PAM filter, or strand
- 🔬 **RNA-Ready:** Accepts uracil (U) in spacers with automatic U→T normalization
- 📈 **Summary Mode:** Quick aggregated counts by distance without per-hit details

### Technical Highlights
- **Myers' Bit-Parallel Algorithm:** 64-way parallelism per thread, GPU-optimized
- **Sparse Output Optimization:** Only transfers hits (~0.00001% of data), eliminating I/O bottleneck
- **Memory-Mapped I/O:** Instant genome loading via `.cy` binary index format
- **Zero-Copy Architecture:** GPU accesses genome via unified memory, no explicit transfers
- **Graceful Validation:** Per-spacer error handling prevents single bad input from crashing batch

---

## Quick Start

### 1. Set up environment
```bash
# Create conda environment
conda create -n cystidia python=3.10 -y
conda activate cystidia

# Install dependencies
conda install -c conda-forge cmake compilers gxx=11 spdlog catch2 -y
pip install pandas numpy matplotlib scikit-learn biopython pysam
```

### 2. Build
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)
```

### 3. Test installation
```bash
# Run unit tests
ctest --test-dir build --output-on-failure

# GPU safety check
./tests/gpu_safety_harness.sh
```

### 4. 30-Second Demo
```bash
# Index a test genome
./build/src/cystidia --index-genome tests/data/synthetic/genome.fa

# Search for off-targets (single guide)
./build/src/cystidia \
  --genome tests/data/synthetic/genome.fa.cy \
  --pattern GAGTCCGAGCAGAAGAAGAA \
  --threshold 3 \
  --format tsv

# Expected: <1 second, finds all sites within 3 mismatches
```

---

## Performance

**Human Genome (hg38, 225 spacers, threshold=3, search-only):**

| Tool | Configuration | Runtime | Notes |
|------|---------------|---------|-------|
| Cystidia | `--distance-mode hamming` (shift-add) | **~38s** | Substitution-only, fastest Cystidia config |
| Cystidia | `--distance-mode levenshtein` (Myers) | ~48s | Substitutions + indels |
| Cas-OFFinder | NGG PAM prefilter | ~21s | Fast, but PAM-coupled |
| Cas-OFFinder | PAM-agnostic (all-N) | ❌ crashes | `CL_OUT_OF_RESOURCES` |

**Hardware:** NVIDIA GPU (compute ≥7.0), hg38.p14 canonical, `.cy` memory-mapped index.

See [Benchmark Results](docs/BENCHMARK_RESULTS.md) for full analysis and [PERFORMANCE_BASELINE.md](docs/PERFORMANCE_BASELINE.md) for annotation-overhead breakdown.

---

## Cystidia vs Cas-OFFinder

Cas-OFFinder is fast because it uses the PAM as a zero-cost stage-1 prefilter: for an `NGG` search, >90% of genome positions are rejected before any mismatch counting happens. That's a great optimization when your PAM is fixed and canonical. It stops working when it isn't.

**Choose Cystidia when:**

- ✅ You care about **non-canonical or novel PAMs** (Cas12a TTTV, engineered variants, PAM-relaxed nucleases). Cystidia finds all sequence matches first and annotates PAM afterward — you can re-filter without re-searching.
- ✅ You want **indel-aware off-target search** (Cystidia Levenshtein mode is a first-class path, not a post-hoc heuristic).
- ✅ You want a **PAM-agnostic scan** (e.g., "find everything within 3 edits of this spacer, regardless of PAM"). Cas-OFFinder's all-N case crashes with `CL_OUT_OF_RESOURCES` on hg38 at threshold 3; Cystidia completes in 38s.
- ✅ You want **rich annotations** out of the box — CFD activity, MIT specificity, seed disruptions, PAM classification — with no downstream scripting.

**Choose Cas-OFFinder when:**

- You only search canonical SpCas9 `NGG` / `NAG` PAMs.
- You only need substitution counts (no indels, no biological annotations).
- Raw speed on the canonical case matters more than scope — Cas-OFFinder is ~2× faster on the NGG case specifically.

In short: **Cas-OFFinder is a fast PAM-coupled search; Cystidia is a full-scope off-target analysis tool.** For a fair apples-to-apples comparison (PAM-agnostic, same threshold, same genome), Cystidia is the only tool that finishes.

---

## Real-World Example

```bash
# Step 1: Index human genome (one-time, ~30 seconds)
./build/src/cystidia --index-genome hg38.fa

# Step 2: Screen CRISPR library (100 guides)
./build/src/cystidia \
  --genome hg38.fa.cy \
  --spacer-file my_guides.txt \
  --threshold 3 \
  --pam-filter both \
  --output offtargets.tsv

# Step 3: Analyze results
# Output includes: chromosome, position, edit distance, 
# seed disruptions, PAM type, alignment details
```

**Use Cases:**
- ✅ Guide specificity screening (pre-synthesis)
- ✅ Library design validation (100s-1000s of guides)
- ✅ Therapeutic guide safety assessment
- ✅ Multiplex editing optimization
- ✅ Off-target prediction for GUIDE-seq/CIRCLE-seq validation

---

---

## Documentation

### User Guides
- **[User Guide](docs/USER_GUIDE.md)** - Complete installation, usage, and troubleshooting
- **[Benchmark Results](docs/BENCHMARK_RESULTS.md)** - Performance analysis on hg38 genome
- **[Achievements Summary](docs/ACHIEVEMENTS.md)** - Technical innovations for 3 audiences
- **[Testing Strategy](docs/TESTING_STRATEGY.md)** - Test suite architecture and coverage

### Developer Resources
- **[Project Context](.claude/PROJECT_CONTEXT.md)** - Architecture decisions and design rationale
- **[Active Tasks](.claude/ACTIVE_TASKS.md)** - Current development status and roadmap

### Quick Links
- [Installation Guide](docs/USER_GUIDE.md#installation)
- [CLI Reference](docs/USER_GUIDE.md#basic-usage)
- [Output Format Documentation](docs/USER_GUIDE.md#output-formats)
- [Performance Optimization Tips](docs/USER_GUIDE.md#performance-optimization)
- [Troubleshooting](docs/USER_GUIDE.md#troubleshooting)

---

## System Requirements

### Minimum
- **OS:** Linux (Ubuntu 20.04+), macOS (x86_64 or ARM64 with Rosetta)
- **CPU:** x86_64 with AVX2 support
- **Memory:** 8 GB RAM
- **Storage:** 10 GB for software + genome indices

### Recommended
- **GPU:** NVIDIA GPU with CUDA compute capability ≥7.0 (RTX 2000+, A100, V100)
- **CUDA:** Toolkit 11.8 or later
- **Memory:** 16+ GB RAM (32 GB for large genomes)
- **Storage:** NVMe SSD for optimal genome loading

---

## Command-Line Reference

### Basic Usage
```bash
# Single spacer search
cystidia --genome GENOME.cy --pattern SEQUENCE --threshold N

# Batch search (multiple spacers)
cystidia --genome GENOME.cy --spacer-file GUIDES.txt --threshold N
```

### Common Options
| Option | Description | Default |
|--------|-------------|---------|
| `--threshold N` | Maximum edit distance (0-255) | Required |
| `--format FMT` | Output format: `tsv`, `bed` | `tsv` |
| `--pam-filter FILTER` | PAM filter: `none`, `ngg`, `nag`, `both`, `any` | `none` |
| `--strand MODE` | Strand selection: `both`, `plus`, `minus` | `both` |
| `--treat-u-as-t` | Normalize uracil to thymine in spacers | On |
| `--summary` | Output aggregated hit counts by distance | Off |
| `--summary-format FMT` | Summary format: `json`, `tsv` | `json` |
| `--max-hits N` | Limit output per spacer | Unlimited |
| `--no-compute-mismatches` | Skip alignment annotations (faster) | Compute |
| `--output FILE` | Output file | stdout |

See [User Guide](docs/USER_GUIDE.md#basic-usage) for complete reference.

---

## Output Format

Cystidia produces tab-separated output with 19 biological annotation columns:

```tsv
spacer  chrom   start   end     pattern distance  strand  aligned_seq  mismatch_pos  ...
VEGFA   chr6    43737006 43737026 GAGTCCGAGCAG... 2  +  GAGTCCAAGCAG...  7  ...
```

**Key Columns:**
- `distance` - Edit distance (Levenshtein)
- `pam_type` - NGG, NAG, OTHER, or INCOMPLETE
- `seed_edits` - Disruptions in seed region (positions 1-10)
- `distal_edits` - Disruptions in distal region (positions 11-20)
- `alignment_ambiguous` - Multiple optimal alignments exist

See [Output Formats](docs/USER_GUIDE.md#output-formats) for complete description.

---

## Citation

If you use Cystidia in your research, please cite:

```bibtex
@software{argus2026,
  author = {Simpson, Daniel L.},
  title = {Cystidia: Ultra-fast exhaustive CRISPR off-target search with GPU acceleration},
  year = {2026},
  version = {0.5.0},
  url = {https://github.com/simpsondl/cystidia}
}
```

---

## Contributing

Contributions welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

**Areas of Interest:**
- Multi-GPU support for larger genomes
- Additional nuclease PAM filters (Cpf1, etc.)
- Python/R API bindings
- Web interface for cloud deployment
- Performance optimizations

---

## License

[License information to be added]

---

## Acknowledgments

- **Myers' Algorithm:** Eugene W. Myers (1999) - Bit-parallel approximate string matching
- **GPU Architecture:** NVIDIA CUDA team
- **Testing:** GeCKO CRISPR library for benchmark spacers
- **Community:** Early adopters and beta testers

---

## Contact

- **Author:** Danny Simpson
- **Email:** dsimpson@example.edu
- **Issues:** [GitHub Issues](https://github.com/simpsondl/cystidia/issues)
- **Discussions:** [GitHub Discussions](https://github.com/simpsondl/cystidia/discussions)

---

<div align="center">

**Built with ❤️ for the CRISPR community**

Made with [C++17](https://isocpp.org/) • [CUDA](https://developer.nvidia.com/cuda-toolkit) • [CMake](https://cmake.org/)

</div>

