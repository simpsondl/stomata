#!/usr/bin/env python3
"""
Clean-room slow oracle implementation of Levenshtein-or-Hamming
off-target search with PAM filtering.

Purpose: an obviously-correct reference that any optimized
implementation (Stomata, SWOffinder, Cas-OFFinder, ...) must agree with
on a small fixture. Designed to catch correctness bugs that even
cross-tool concordance might miss, because Stomata and SWOffinder share
neither codebase nor algorithmic approach but might still happen to
agree on the wrong answer for some specific edge case.

Trade-offs:
  - SLOW. O(genome × spacers × |tl_range| × pattern_len × text_len) DP.
    Suitable for fixtures of ≤100 spacers × ≤1 Mb genome. Not for hg38.
  - SIMPLE. Bog-standard Needleman-Wunsch DP, no bit-parallelism, no
    GPU, no CIGAR construction. Just min edit distance + PAM check.
  - CONSERVATIVE. Reports every position where (a) a valid d ≤ threshold
    alignment exists with (b) a PAM that matches the IUPAC pattern at
    the canonical or any bulge-shifted target length in
    [pattern_len ± threshold]. This matches the v0.8.0 alignment-validity
    semantics — a hit needs both PAM match AND a valid alignment.

Output: TSV with columns
  spacer  chrom  start  end  strand  distance  pam_seq  target_len

Compare to Stomata's BED output by sorting both on (spacer, chrom, start,
strand) and diffing. Bit-equivalent agreement is the success criterion;
any mismatch is a candidate bug in either implementation.

Usage:
  python3 scripts/oracle_lev_search.py \\
      --spacers tests/regression/oracle_fixture.tsv \\
      --genome tests/regression/oracle_genome.fa \\
      --threshold 3 \\
      --pam NGG \\
      --distance-mode levenshtein \\
      --out oracle.tsv
"""

import argparse
import sys
from pathlib import Path

# IUPAC matching for PAM patterns.
# Returns True iff observation `o` is consistent with pattern code `p`.
def iupac_match(p: str, o: str) -> bool:
    p = p.upper()
    o = o.upper()
    table = {
        "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"}, "U": {"T"},
        "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"C", "G"},
        "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
        "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
        "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
        "N": {"A", "C", "G", "T"},
    }
    return o in table.get(p, set())


def pam_matches(pam: str, pattern: str) -> bool:
    if len(pam) < len(pattern):
        return False
    for pc, oc in zip(pattern, pam):
        if not iupac_match(pc, oc):
            return False
    return True


# Reverse complement.
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# Standard global Needleman-Wunsch edit distance with row-min early exit.
# Treats pattern N as wildcard (cost 0 to anything), genome N as masked
# (cost 1 to anything). Same semantics as src/search_pipeline.cpp's
# alignment_exists_at_target_length.
def edit_distance(pattern: str, target: str, threshold: int) -> int:
    m, n = len(pattern), len(target)
    INF = threshold + 1
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        curr = [i] + [0] * n
        row_min = i
        p = pattern[i - 1].upper()
        for j in range(1, n + 1):
            t = target[j - 1].upper()
            if p == "N":
                cost = 0
            elif t == "N":
                cost = 1
            else:
                cost = 0 if p == t else 1
            curr[j] = min(prev[j - 1] + cost, prev[j] + 1, curr[j - 1] + 1)
            if curr[j] < row_min:
                row_min = curr[j]
        if row_min > threshold:
            return INF
        prev = curr
    return prev[n]


def hamming_distance(pattern: str, target: str, threshold: int) -> int:
    """Substitution-only distance, length-fixed (len pattern == len target)."""
    if len(pattern) != len(target):
        return threshold + 1
    d = 0
    for p, t in zip(pattern, target):
        p, t = p.upper(), t.upper()
        if p == "N":
            continue
        if t == "N":
            d += 1
        elif p != t:
            d += 1
        if d > threshold:
            return threshold + 1
    return d


def load_fasta(path: Path) -> dict:
    """Return {chrom_name: sequence_string} (uppercase, no Ns stripped)."""
    chroms = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    chroms[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            elif line.startswith(";") or not line:
                continue
            else:
                chunks.append(line)
    if name is not None:
        chroms[name] = "".join(chunks).upper()
    return chroms


def search_one_strand(spacer_name: str, pattern: str,
                     chrom: str, seq: str, strand: str,
                     threshold: int, pam_pat: str, pam_position: str,
                     distance_mode: str) -> list:
    """
    Naive scan with CONSERVATIVE semantics: a hit at end position E exists
    iff there is at least one target length `tl ∈ [pattern_len ± threshold]`
    (Lev) or `tl == pattern_len` (Hamming) such that BOTH:
      (a) the alignment of `pattern` to `seq[start..end_pos+1]` has
          edit distance ≤ threshold, AND
      (b) the PAM at the position implied by that same tl matches the
          IUPAC pattern.

    These two conditions must hold for the *same* tl. This is the
    biological semantics: Cas9 engages a single physical alignment with
    a flanking PAM; "different alignment for distance vs PAM" is a
    confused row, not a real hit.

    For ranking among multiple valid tls, pick the lowest-d one (and
    among ties, the canonical tl == pattern_len if it matches, else
    the smallest tl).
    """
    m = len(pattern)
    pam_len = len(pam_pat) if pam_pat else 0

    # tl range: Lev allows ±threshold; Hamming locks to m.
    if distance_mode == "levenshtein":
        tl_lo = max(1, m - threshold)
        tl_hi = m + threshold
    else:
        tl_lo = m
        tl_hi = m

    hits = []
    n = len(seq)

    # Iterate over end positions (0-based inclusive).
    # Conservative single-tl filter: for each end_pos, scan all tls,
    # require BOTH d ≤ threshold AND matching PAM at the same tl.
    for end_pos in range(m - 1, n):
        best_d = threshold + 1
        best_tl = None
        best_pam_seq = None

        for tl in range(tl_lo, tl_hi + 1):
            start = end_pos - tl + 1
            if start < 0:
                continue
            target = seq[start:end_pos + 1]

            # Distance check
            if distance_mode == "levenshtein":
                d = edit_distance(pattern, target, threshold)
            else:
                d = hamming_distance(pattern, target, threshold)
            if d > threshold:
                continue

            # PAM check at THIS tl
            if pam_pat:
                if pam_position == "3prime":
                    pam_start_idx = end_pos + 1
                    pam_end_idx = pam_start_idx + pam_len
                    if pam_end_idx > n:
                        continue
                    pam_seq = seq[pam_start_idx:pam_end_idx]
                else:  # 5prime
                    pam_end_idx = start
                    pam_start_idx = pam_end_idx - pam_len
                    if pam_start_idx < 0:
                        continue
                    pam_seq = seq[pam_start_idx:pam_end_idx]
                if not pam_matches(pam_seq, pam_pat):
                    continue
            else:
                pam_seq = ""

            # Both conditions hold at this tl. Track the lowest-d valid (tl, d, PAM).
            # Tie-break: prefer canonical tl == pattern_len, else smallest tl.
            if d < best_d or (
                d == best_d and (
                    (best_tl != m and tl == m) or
                    (best_tl != m and tl < best_tl)
                )
            ):
                best_d = d
                best_tl = tl
                best_pam_seq = pam_seq

        if best_d > threshold:
            continue

        nominal_start = max(0, end_pos - m + 1)
        hits.append({
            "spacer": spacer_name,
            "chrom": chrom,
            "start": nominal_start,
            "end": end_pos + 1,  # half-open
            "strand": strand,
            "distance": best_d,
            "pam_seq": best_pam_seq if best_pam_seq is not None else "",
            "target_len": best_tl,
        })

    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spacers", required=True, type=Path,
                    help="TSV with name<TAB>sequence rows; '#' comment lines OK")
    ap.add_argument("--genome", required=True, type=Path,
                    help="FASTA file (uncompressed)")
    ap.add_argument("--threshold", type=int, default=3)
    ap.add_argument("--pam", default="NGG", help="IUPAC PAM pattern; '' = no filter")
    ap.add_argument("--pam-position", default="3prime", choices=["3prime", "5prime"])
    ap.add_argument("--distance-mode", default="levenshtein",
                    choices=["levenshtein", "hamming"])
    ap.add_argument("--strand", default="both", choices=["both", "plus", "minus"])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    print(f"Oracle Lev/Hamming search", file=sys.stderr)
    print(f"  pattern threshold: {args.threshold}", file=sys.stderr)
    print(f"  pam:               {args.pam or '(none)'}", file=sys.stderr)
    print(f"  distance:          {args.distance_mode}", file=sys.stderr)
    print(f"  strand:            {args.strand}", file=sys.stderr)

    # Load spacers
    spacers = []
    with open(args.spacers) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            name, seq = parts[0], parts[1]
            spacers.append((name, seq.upper()))
    print(f"  spacers: {len(spacers)}", file=sys.stderr)

    # Load genome
    chroms = load_fasta(args.genome)
    total_bp = sum(len(s) for s in chroms.values())
    print(f"  genome: {total_bp:,} bp across {len(chroms)} chrom(s)", file=sys.stderr)

    # Run
    all_hits = []
    for sp_name, sp_seq in spacers:
        print(f"  ... searching {sp_name} ({sp_seq})", file=sys.stderr)
        for chrom_name, chrom_seq in chroms.items():
            if args.strand in ("both", "plus"):
                all_hits.extend(search_one_strand(
                    sp_name, sp_seq, chrom_name, chrom_seq, "+",
                    args.threshold, args.pam, args.pam_position,
                    args.distance_mode))
            if args.strand in ("both", "minus"):
                # For minus strand, scan the RC of the chromosome with the
                # plus-strand pattern. Coordinates need flipping back.
                rc_seq = rc(chrom_seq)
                rc_hits = search_one_strand(
                    sp_name, sp_seq, chrom_name, rc_seq, "-",
                    args.threshold, args.pam, args.pam_position,
                    args.distance_mode)
                # Flip start/end coordinates back to the forward-strand frame
                n_chrom = len(chrom_seq)
                for h in rc_hits:
                    fwd_start = n_chrom - h["end"]
                    fwd_end = n_chrom - h["start"]
                    h["start"] = fwd_start
                    h["end"] = fwd_end
                all_hits.extend(rc_hits)

    # Halo dedup with TRANSITIVE cluster extension (matches Stomata's
    # deduplicate_hits in src/search_pipeline.cpp). For sorted hits in
    # the same (spacer, chrom, strand), extend a cluster while
    # consecutive end-positions are within `threshold + 1` of EACH OTHER
    # (not of the cluster pivot). This means {829, 831, 833} with
    # radius=3 forms ONE cluster (829↔831 within, 831↔833 within), even
    # though 829↔833 are NOT within radius.
    if args.distance_mode == "levenshtein":
        dedup_radius = args.threshold + 1
        all_hits.sort(key=lambda h: (h["spacer"], h["chrom"], h["strand"], h["end"], h["distance"]))
        deduped = []
        i = 0
        while i < len(all_hits):
            cluster_start = i
            best_idx = i
            j = i + 1
            while j < len(all_hits) \
                    and all_hits[j]["spacer"] == all_hits[i]["spacer"] \
                    and all_hits[j]["chrom"] == all_hits[i]["chrom"] \
                    and all_hits[j]["strand"] == all_hits[i]["strand"] \
                    and (all_hits[j]["end"] - all_hits[j-1]["end"]) < dedup_radius:
                # Stomata's tie-break: lower distance > better PAM > rightmost.
                # We don't track PAM priority in the oracle (it's a SpCas9
                # heuristic for halo dedup that has nothing to do with
                # whether the alignment is correct), but we do mirror the
                # rightmost-on-distance-tie. This means with d=2 hits at
                # 829, 831, 833 in one cluster, oracle picks 833 (matches
                # Stomata); without this, oracle picked 829.
                if all_hits[j]["distance"] < all_hits[best_idx]["distance"]:
                    best_idx = j
                elif all_hits[j]["distance"] == all_hits[best_idx]["distance"]:
                    best_idx = j  # rightmost
                j += 1
            deduped.append(all_hits[best_idx])
            i = j
        all_hits = deduped

    # Sort for deterministic comparison
    all_hits.sort(key=lambda h: (h["spacer"], h["chrom"], h["start"], h["strand"]))

    # Emit TSV
    with args.out.open("w") as out:
        out.write("spacer\tchrom\tstart\tend\tstrand\tdistance\tpam_seq\ttarget_len\n")
        for h in all_hits:
            out.write(f"{h['spacer']}\t{h['chrom']}\t{h['start']}\t{h['end']}\t"
                      f"{h['strand']}\t{h['distance']}\t{h['pam_seq']}\t{h['target_len']}\n")

    print(f"  wrote {len(all_hits):,} hits → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
