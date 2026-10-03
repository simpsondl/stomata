# Stomata

GPU-accelerated CRISPR off-target sequence search, with mismatches, indels,
and arbitrary IUPAC PAMs.

**Version:** 0.11.0 - **License:** MIT - **Platform:** Linux, C++17

Stomata scans a reference genome for sequences within a specified Hamming or
Levenshtein distance of one guide or a batch of guides. It can search without
a PAM restriction, or filter for a 3-prime or 5-prime PAM. Output includes coordinates,
strand, alignment and PAM annotations, optional CFD scores, and guide summaries.

The search does not depend on a seed index or a fixed PAM catalog. A memory-mapped
`.st` reference index avoids repeated FASTA parsing. CUDA accelerates searches
on NVIDIA GPUs; a CPU-only build is also available.

## Install and try it

On Ubuntu 24.04, the CPU build needs no CUDA toolkit:

```bash
git clone https://github.com/simpsondl/stomata.git
cd stomata
sudo apt-get update
sudo apt-get install -y cmake ninja-build g++ libspdlog-dev zlib1g-dev catch2
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DSTOMATA_ENABLE_CUDA=OFF
cmake --build build -j2
./build/src/stomata --quickstart --cpu-only
```

For a GPU build, install CUDA Toolkit 11.8+ and omit
`-DSTOMATA_ENABLE_CUDA=OFF`. See [installation](docs/INSTALL.md) for the
container, supported build targets, and Docker/Apptainer use.

```bash
# A small bundled example; no genome download needed.
./build/src/stomata --genome tests/data/synthetic/test_small.fa --pattern ACGTACGT --threshold 1 --output demo.tsv

# Index a reference once, then search one guide or a batch.
./build/src/stomata --index-genome hg38.fa
./build/src/stomata --genome hg38.fa.st --pattern GAGTCCGAGCAGAAGAAGAA --threshold 3 --pam NGG --output hits.tsv
# Replace --pattern with --spacer-file guides.tsv for batch search.
```

Spacer files accept a sequence per line or `name<TAB>sequence`. Use
`--spacer-summary summary.tsv` for guide-level counts and aCFD.

## Search and output conventions

- Default distance is Levenshtein, threshold 3, both strands, **no PAM filter**.
- Levenshtein output is halo-deduplicated by default. `--no-deduplicate` exposes
  raw end-position hits; nearby reports need not be independent biological loci.
- Hit caps deliberately truncate output after search/annotation. They are not
  memory limits and cannot be combined with aggregate summaries.
- Hamming search uses Levenshtein traceback for annotations. Use
  `--no-compute-mismatches` for strict Hamming-distance reporting without traceback.
- RNA `U` is normalized to `T`. Pattern `N` is a wildcard; reference `N` is masked.
- A sequence candidate is not evidence of cleavage. CFD/aCFD have their own
  model assumptions; see the [user guide](docs/USER_GUIDE.md).

## Documentation and reproducibility

- [User guide](docs/USER_GUIDE.md): flags, formats, coordinates and examples
- [Installation and testing](docs/INSTALL.md): builds, CPU-only option, containers, tests
- [Validation](validation/README.md): independent DP oracle and fuzzing, cross-tool and
  biological evaluations
- [Benchmarks and plot](validation/benchmarks/README.md): four search scopes,
  recorded commands, input hashes and repeated measurements
- [Citation](CITATION.cff)

Issues and reproducible bug reports: <https://github.com/simpsondl/stomata/issues>.
