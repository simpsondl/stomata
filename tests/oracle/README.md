# Independent correctness checks

Run these checks from the repository root with Python 3.10+ and a built
Stomata executable. They use small synthetic inputs and need no downloads
or additional Python packages.

```bash
python3 tests/oracle/check_fixture.py --binary build/src/stomata --cpu-only
python3 tests/oracle/fuzz_oracle_vs_stomata.py \
  --stomata-bin build/src/stomata --num-configs 200 --seed-start 1 \
  --out-dir tests/oracle/failures --cpu-only
```

The fixture compares Hamming and Levenshtein search with NGG and unrestricted
PAM. Results are written under `tests/oracle/results/` (ignored by Git).
The fuzz tester varies PAM side, IUPAC pattern, threshold, guide length,
masked bases and repeats, saving failing inputs and outputs for inspection.

On a GPU host with a CUDA build, repeat without `--cpu-only` and confirm
that the CLI reports GPU execution. CPU and GPU search are separate paths.

These checks compare against an independent dynamic-programming distance
calculation. When alignments of different lengths tie at the best distance,
the oracle lists every optimal span and any of them counts as a match. The oracle shares some reporting conventions, including halo
clustering; disagreements can reflect traceback or deduplication as well
as search errors. Preserve and investigate failing cases.

`tests/regression/smoking_guns.tsv` contains historical cases that require
an external reference genome. Its old counts are not release expectations
and are not used by the synthetic fixture check.
