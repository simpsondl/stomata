#!/usr/bin/env python3
"""
Property-based fuzz tester: oracle (slow, obviously correct) vs Stomata
(fast, GPU). Generates random (genome, spacer-pool, threshold, PAM,
distance-mode) tuples, runs both, diffs. Any disagreement is saved as
a complete reproducer for offline investigation.

Designed to find bugs we don't yet know to look for — we've already
found three correctness bugs through targeted experiments; this
script is the "what else is hiding" net.

Output:
  - Real-time progress to stderr.
  - Per-config result line (PASS / FAIL / TIMEOUT / SKIP).
  - On FAIL: full reproducer saved to OUT_DIR/seed_NNNN/
            (genome.fa + spacers.tsv + config.json + oracle.tsv +
             stomata.tsv + diff.txt).
  - End-of-run summary.

Usage:
  python3 tests/oracle/fuzz_oracle_vs_stomata.py \\
      --num-configs 500 \\
      --seed-start 1 \\
      --out-dir tests/oracle/failures/ \\
      --stomata-bin build/src/stomata

Run overnight. Each config takes ~5-30 s (oracle bound). 500 configs
≈ 1-4 hours.
"""

import argparse
import csv
import json
import random
import shutil
import string
import subprocess
import sys
import tempfile
from pathlib import Path

# Configuration space — cartesian product is sampled, not enumerated.
# Each axis has ~uniform sampling weight unless skewed for known-suspect cases.
THRESHOLDS = [0, 1, 2, 3, 4, 5]
DISTANCE_MODES = ["hamming", "levenshtein"]
PAM_PATTERNS = ["", "NGG", "NAG", "NRG", "NNGRRT", "TTTV"]
PAM_POSITIONS = ["3prime", "5prime"]
SPACER_LENGTHS = [18, 19, 20, 21, 22, 23, 24]  # covers uint32 + uint64 branches
N_SPACERS_RANGE = (3, 8)
N_CHROMS_RANGE = (1, 4)
CHROM_LENGTH_RANGE = (200, 4000)
INCLUDE_N_RUNS = True   # 30 % of chromosomes get an N-run somewhere
INCLUDE_REPEATS = True  # 20 % of chromosomes get a homopolymer / dinuc repeat


def random_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choices("ACGT", k=n))


def random_chromosome(rng: random.Random, length: int) -> str:
    """Generate a chromosome with optional N-runs and repetitive regions
    interspersed in random ACGT."""
    seq = list(random_dna(rng, length))
    if INCLUDE_N_RUNS and rng.random() < 0.3 and length > 50:
        run_len = rng.randint(5, min(40, length // 5))
        run_start = rng.randint(10, length - run_len - 10)
        for i in range(run_start, run_start + run_len):
            seq[i] = "N"
    if INCLUDE_REPEATS and rng.random() < 0.2 and length > 60:
        rep_unit = rng.choice(["A", "T", "AT", "AG", "AAG", "GGG"])
        rep_len = rng.randint(20, min(80, length // 3))
        rep_start = rng.randint(10, length - rep_len - 10)
        rep_seq = (rep_unit * (rep_len // len(rep_unit) + 1))[:rep_len]
        for i in range(rep_start, rep_start + rep_len):
            seq[i] = rep_seq[i - rep_start]
    return "".join(seq)


def random_spacer(rng: random.Random, length: int,
                   chroms: dict, embed_match_chance: float = 0.5) -> str:
    """Generate a spacer. With probability embed_match_chance, copy a
    real substring from a chromosome (so we get realistic hit densities);
    otherwise pure random."""
    if embed_match_chance > 0 and chroms and rng.random() < embed_match_chance:
        chrom = rng.choice(list(chroms.values()))
        if len(chrom) >= length:
            start = rng.randint(0, len(chrom) - length)
            substr = chrom[start:start + length]
            if "N" not in substr:
                # Optionally introduce 0-2 substitutions
                n_subs = rng.randint(0, 2)
                substr_l = list(substr)
                for _ in range(n_subs):
                    pos = rng.randint(0, length - 1)
                    substr_l[pos] = rng.choice([b for b in "ACGT"
                                                if b != substr_l[pos]])
                return "".join(substr_l)
    return random_dna(rng, length)


def gen_config(seed: int) -> dict:
    rng = random.Random(seed)
    n_chroms = rng.randint(*N_CHROMS_RANGE)
    chroms = {}
    for i in range(n_chroms):
        length = rng.randint(*CHROM_LENGTH_RANGE)
        chroms[f"chr_{i+1}"] = random_chromosome(rng, length)

    n_spacers = rng.randint(*N_SPACERS_RANGE)
    spacers = []
    for i in range(n_spacers):
        sp_len = rng.choice(SPACER_LENGTHS)
        spacers.append((f"sp_{i+1:02d}", random_spacer(rng, sp_len, chroms)))

    return {
        "seed":          seed,
        "threshold":     rng.choice(THRESHOLDS),
        "distance_mode": rng.choice(DISTANCE_MODES),
        "pam":           rng.choice(PAM_PATTERNS),
        "pam_position":  rng.choice(PAM_POSITIONS),
        "chroms":        chroms,
        "spacers":       spacers,
    }


def write_inputs(cfg: dict, out_dir: Path):
    """Write FASTA + spacer-list to out_dir; return their paths."""
    fasta_path = out_dir / "genome.fa"
    with fasta_path.open("w") as f:
        for name, seq in cfg["chroms"].items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")

    spacer_path = out_dir / "spacers.tsv"
    with spacer_path.open("w") as f:
        for name, seq in cfg["spacers"]:
            f.write(f"{name}\t{seq}\n")

    cfg_path = out_dir / "config.json"
    with cfg_path.open("w") as f:
        json.dump({k: v for k, v in cfg.items() if k != "chroms"}, f,
                  indent=2, default=str)

    return fasta_path, spacer_path, cfg_path


def run_stomata(stomata_bin: Path, fasta: Path, spacer_file: Path, cfg: dict,
              out_tsv: Path, timeout_s: int = 60) -> str:
    """Returns 'OK' or an error string."""
    args = [
        str(stomata_bin),
        "--genome", str(fasta),
        "--spacer-file", str(spacer_file),
        "--threshold", str(cfg["threshold"]),
        "--distance-mode", cfg["distance_mode"],
        "--output", str(out_tsv),
    ]
    if cfg["pam"]:
        args += ["--pam", cfg["pam"], "--pam-position", cfg["pam_position"]]
    if cfg.get("cpu_only"):
        args += ["--cpu-only"]
    try:
        r = subprocess.run(args, capture_output=True, text=True,
                           timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return f"TIMEOUT_after_{timeout_s}s"
    if r.returncode != 0:
        return f"stomata exit={r.returncode}: {r.stderr.strip().splitlines()[-1] if r.stderr else ''}"
    return "OK"


def run_oracle(oracle_py: Path, fasta: Path, spacer_file: Path, cfg: dict,
               out_tsv: Path, timeout_s: int = 300) -> str:
    args = [
        sys.executable, str(oracle_py),
        "--spacers", str(spacer_file),
        "--genome", str(fasta),
        "--threshold", str(cfg["threshold"]),
        "--pam", cfg["pam"],
        "--pam-position", cfg["pam_position"],
        "--distance-mode", cfg["distance_mode"],
        "--out", str(out_tsv),
    ]
    try:
        r = subprocess.run(args, capture_output=True, text=True,
                           timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return f"TIMEOUT_after_{timeout_s}s"
    if r.returncode != 0:
        return f"oracle exit={r.returncode}: {r.stderr.strip().splitlines()[-1] if r.stderr else ''}"
    return "OK"


def canonicalize(tsv_path: Path, source: str) -> list:
    """Read a hit TSV, return a sorted list of (spacer, chrom, start, end,
    strand, distance) tuples. Different sources have slightly different
    columns; pick the right ones."""
    return sorted(row[:6] for row in _read_rows(tsv_path))


def _read_rows(tsv_path: Path) -> set:
    """Rows as (spacer, chrom, start, end, strand, distance, spans), where
    spans holds every optimal (start, end) when the oracle reports them."""
    rows = set()
    with tsv_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"spacer", "chrom", "start", "end", "strand", "distance"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"Missing columns in {tsv_path}")
        for row in reader:
            start, end = int(row["start"]), int(row["end"])
            spans = {(start, end)}
            for span in filter(None, (row.get("optimal_spans") or "").split(",")):
                s, e = span.split("-")
                spans.add((int(s), int(e)))
            rows.add((row["spacer"], row["chrom"], start, end, row["strand"],
                      int(row["distance"]), frozenset(spans)))
    return rows


def compare_files(stomata_tsv: Path, oracle_tsv: Path) -> str:
    """Diff two hit files. A Stomata hit matches an oracle hit with the same
    spacer, chromosome, strand and distance whose optimal spans include the
    Stomata span: co-optimal alignments of different lengths are all correct."""
    stomata_rows = sorted(r[:6] for r in _read_rows(stomata_tsv))
    unmatched_oracle = sorted(_read_rows(oracle_tsv), key=lambda r: r[:6])
    n_oracle = len(unmatched_oracle)
    only_stomata = []
    for a in stomata_rows:
        for i, o in enumerate(unmatched_oracle):
            if a[:2] == o[:2] and a[4:6] == o[4:6] and (a[2], a[3]) in o[6]:
                del unmatched_oracle[i]
                break
        else:
            only_stomata.append(a)
    return _format_diff(len(stomata_rows), n_oracle, only_stomata,
                        [o[:6] for o in unmatched_oracle])


def _format_diff(n_stomata: int, n_oracle: int, only_stomata: list, only_oracle: list) -> str:
    if not only_stomata and not only_oracle:
        return ""
    lines = []
    lines.append(f"stomata has {n_stomata} hits, oracle has {n_oracle} hits")
    lines.append(f"only in stomata ({len(only_stomata)}):")
    for r in only_stomata[:30]:
        lines.append(f"  + {r}")
    if len(only_stomata) > 30:
        lines.append(f"  ... +{len(only_stomata) - 30} more")
    lines.append(f"only in oracle ({len(only_oracle)}):")
    for r in only_oracle[:30]:
        lines.append(f"  - {r}")
    if len(only_oracle) > 30:
        lines.append(f"  ... +{len(only_oracle) - 30} more")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--num-configs", type=int, default=200)
    ap.add_argument("--seed-start", type=int, default=1)
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="Directory to save failure reproducers")
    ap.add_argument("--stomata-bin", required=True, type=Path)
    ap.add_argument("--cpu-only", action="store_true")
    ap.add_argument("--oracle-py", default=Path(__file__).parent / "oracle_lev_search.py", type=Path)
    ap.add_argument("--stomata-timeout", type=int, default=60)
    ap.add_argument("--oracle-timeout", type=int, default=300)
    args = ap.parse_args()
    if args.num_configs < 1:
        ap.error("--num-configs must be positive")
    args.stomata_bin = args.stomata_bin.resolve()

    if not args.stomata_bin.exists():
        sys.exit(f"stomata binary not found: {args.stomata_bin}")
    if not args.oracle_py.exists():
        sys.exit(f"oracle script not found: {args.oracle_py}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Log the build version we're testing
    version = subprocess.run([str(args.stomata_bin), "--version"],
                              capture_output=True, text=True).stdout.strip()
    print(f"# Fuzz run vs {version}", file=sys.stderr)
    print(f"# Configs: seed {args.seed_start} → {args.seed_start + args.num_configs - 1}",
          file=sys.stderr)
    print(f"# Failures saved to: {args.out_dir}", file=sys.stderr)
    print(f"# {'-' * 60}", file=sys.stderr)
    print(f"# {'seed':>5} {'mode':>12} {'pam':>8} {'thr':>3} {'pam_pos':>7} {'spacers':>7} {'stomata_n':>7} {'oracle_n':>8} {'result':>8}",
          file=sys.stderr)

    n_pass = 0
    n_fail = 0
    n_skip = 0
    n_stomata_err = 0
    n_oracle_err = 0
    failures = []

    for i in range(args.num_configs):
        seed = args.seed_start + i
        cfg = gen_config(seed)
        cfg["cpu_only"] = args.cpu_only

        # The oracle and Stomata disagree on the meaning of `--pam ""`
        # passed via the --pam= argument: Stomata accepts it, oracle accepts
        # it. Both should produce no PAM filtering. Skip if either rejects.
        # Skip cfgs where pam_len > pattern_len (no valid alignment possible).
        max_pat_len = max(len(s) for _, s in cfg["spacers"])
        if cfg["pam"] and len(cfg["pam"]) > max_pat_len:
            n_skip += 1
            print(f"  {seed:>5} {cfg['distance_mode']:>12} {cfg['pam']:>8} "
                  f"{cfg['threshold']:>3} {cfg['pam_position']:>7} "
                  f"{len(cfg['spacers']):>7} {'-':>7} {'-':>8} {'SKIP':>8} (pam>pat)",
                  file=sys.stderr)
            continue

        # Stomata rejects --pam + --no-compute-mismatches in Lev mode; we
        # always compute mismatches anyway, so this should never trip.

        with tempfile.TemporaryDirectory(prefix=f"fuzz_seed{seed}_") as td_str:
            td = Path(td_str)
            fasta, spacer_file, _ = write_inputs(cfg, td)
            stomata_out = td / "stomata.tsv"
            oracle_out = td / "oracle.tsv"

            arg_status = run_stomata(args.stomata_bin, fasta, spacer_file, cfg,
                                    stomata_out, args.stomata_timeout)
            ora_status = run_oracle(args.oracle_py, fasta, spacer_file, cfg,
                                     oracle_out, args.oracle_timeout)

            if arg_status != "OK" or ora_status != "OK":
                # Surface as a failure (could be a real bug, e.g. stomata crash)
                if "stomata" in arg_status: n_stomata_err += 1
                if "oracle" in ora_status: n_oracle_err += 1
                # Save the inputs anyway for triage
                fail_dir = args.out_dir / f"seed_{seed:05d}_error"
                fail_dir.mkdir(exist_ok=True)
                shutil.copy(fasta, fail_dir / "genome.fa")
                shutil.copy(spacer_file, fail_dir / "spacers.tsv")
                with (fail_dir / "config.json").open("w") as f:
                    json.dump({k: v for k, v in cfg.items() if k != "chroms"}, f,
                              indent=2, default=str)
                with (fail_dir / "status.txt").open("w") as f:
                    f.write(f"stomata: {arg_status}\noracle: {ora_status}\n")
                n_fail += 1
                print(f"  {seed:>5} {cfg['distance_mode']:>12} {cfg['pam']:>8} "
                      f"{cfg['threshold']:>3} {cfg['pam_position']:>7} "
                      f"{len(cfg['spacers']):>7} ERR     ERR      ERR     "
                      f"stomata:{arg_status[:30]} oracle:{ora_status[:30]}",
                      file=sys.stderr)
                failures.append(seed)
                continue

            stomata_rows = canonicalize(stomata_out, "stomata")
            oracle_rows = canonicalize(oracle_out, "oracle")

            diff_text = compare_files(stomata_out, oracle_out)
            if not diff_text:
                n_pass += 1
                print(f"  {seed:>5} {cfg['distance_mode']:>12} {cfg['pam']:>8} "
                      f"{cfg['threshold']:>3} {cfg['pam_position']:>7} "
                      f"{len(cfg['spacers']):>7} {len(stomata_rows):>7} "
                      f"{len(oracle_rows):>8} {'PASS':>8}", file=sys.stderr)
            else:
                n_fail += 1
                fail_dir = args.out_dir / f"seed_{seed:05d}"
                fail_dir.mkdir(exist_ok=True)
                shutil.copy(fasta, fail_dir / "genome.fa")
                shutil.copy(spacer_file, fail_dir / "spacers.tsv")
                shutil.copy(stomata_out, fail_dir / "stomata.tsv")
                shutil.copy(oracle_out, fail_dir / "oracle.tsv")
                with (fail_dir / "config.json").open("w") as f:
                    json.dump({k: v for k, v in cfg.items() if k != "chroms"}, f,
                              indent=2, default=str)
                with (fail_dir / "diff.txt").open("w") as f:
                    f.write(diff_text + "\n")
                print(f"  {seed:>5} {cfg['distance_mode']:>12} {cfg['pam']:>8} "
                      f"{cfg['threshold']:>3} {cfg['pam_position']:>7} "
                      f"{len(cfg['spacers']):>7} {len(stomata_rows):>7} "
                      f"{len(oracle_rows):>8} {'FAIL':>8}  → {fail_dir}",
                      file=sys.stderr)
                failures.append(seed)

    print(f"# {'-' * 60}", file=sys.stderr)
    print(f"# Result: {n_pass} pass, {n_fail} fail, {n_skip} skip "
          f"({n_stomata_err} stomata errors, {n_oracle_err} oracle errors)",
          file=sys.stderr)
    if failures:
        print(f"# Failure seeds: {failures[:20]}{' ...' if len(failures) > 20 else ''}",
              file=sys.stderr)
        print(f"# Reproducers in: {args.out_dir}", file=sys.stderr)

    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
