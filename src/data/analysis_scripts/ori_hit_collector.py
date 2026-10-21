#!/usr/bin/env python3
"""
Collect the top-N local-alignment hits of an ori query against virus genomes.

Layout expected:
    virus_genomes/<subfolder>/<virus_folder>/*.fasta (or .fa/.fna/.fasta)

For every virus folder in the chosen subfolder, the first fasta file (sorted by name)
is used, and within it only the first record. Results are APPENDED to one CSV, one row
per virus:

    matrix_class, ori_1, score_1, position_1, ori_2, score_2, position_2,
    ori_3, score_3, position_3, genome_length

Example:
    uv run src/data/analysis/ori_hit_collector.py --subfolder WDV_iter1_V_R --matrix repensiviricetes
    uv run src/data/analysis/ori_hit_collector.py --subfolder PCV_iter1_V_A --matrix arfiviricetes
"""
import argparse
import csv
import math
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

# ----------------------------- CONFIG (edit me) -----------------------------
GENOMES_ROOT = Path("data/raw/downloads/virus_genomes")
MATRIX_DIR = Path("src/data/scoring_matrices")             # holds files named after the matrix classes
MATRIX_FILES = {                          # matrix class -> file in MATRIX_DIR
    "repensiviricetes": "repensiviricetes",
    "arfiviricetes": "arfiviricetes",
}
OUT_CSV = Path("data/processed/rep_ori_hits.csv")
FASTA_EXTS = {".fasta", ".fa", ".fna", ".fas"}
MASK_CHAR = "N"                           # must be in the matrix alphabet (and not occur in the query)
OPEN_GAP, EXTEND_GAP = -math.inf, -math.inf
QUERY = "HMDWNYTAC"
ORI_LEN = len(QUERY)                      # 9
N_HITS = 3
# ----------------------------------------------------------------------------

FIELDS = (["matrix_class"]
          + [f"{col}_{i}" for i in range(1, N_HITS + 1) for col in ("ori", "score", "position")]
          + ["genome_length"])


def load_matrix(matrix_class):
    return substitution_matrices.read(str(MATRIX_DIR / MATRIX_FILES[matrix_class]))


def make_aligner(matrix):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = OPEN_GAP
    aligner.extend_gap_score = EXTEND_GAP
    return aligner


def top_hits(aligner, genome, query, n=N_HITS):
    """Return [(ori_9mer, score, position), ...] for the n best non-overlapping hits.

    The query ends in the fixed "AC" terminus, so every scoring alignment ends at the
    query's last base. The reported ori is therefore the 9-mer window of the genome
    ending where the alignment ends, which is always exactly ORI_LEN bases long even if
    the best local alignment only covered part of the query.
    position = 0-based start of that 9-mer in the genome.
    """
    masked = list(genome)
    hits = []
    for _ in range(n):
        masked_str = "".join(masked)
        score = aligner.score(masked_str, query)
        if score <= 0:
            break
        # next(iter(...)) avoids enumerating (possibly huge) numbers of tied alignments
        aln = next(iter(aligner.align(masked_str, query)))
        end = int(aln.coordinates[0][-1])
        start = end - ORI_LEN
        if start < 0:  # hit too close to the genome start to form a full 9-mer
            break
        hits.append((genome[start:end], score, start))
        masked[start:end] = [MASK_CHAR] * ORI_LEN  # same length -> coords stay valid
    return hits


def first_fasta(virus_dir):
    files = sorted(p for p in virus_dir.iterdir() if p.suffix.lower() in FASTA_EXTS)
    return files[0] if files else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subfolder", required=True, help="name of subfolder inside virus_genomes")
    ap.add_argument("--matrix", required=True, choices=sorted(MATRIX_FILES))
    ap.add_argument("--out", type=Path, default=OUT_CSV)
    args = ap.parse_args()

    matrix = load_matrix(args.matrix)
    if MASK_CHAR not in matrix.alphabet:
        raise SystemExit(f"Mask character {MASK_CHAR!r} is not in the matrix alphabet "
                         f"({matrix.alphabet}); pick another MASK_CHAR.")
    aligner = make_aligner(matrix)

    sub_dir = GENOMES_ROOT / args.subfolder
    if not sub_dir.is_dir():
        raise SystemExit(f"No such folder: {sub_dir}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_header = not args.out.exists() or args.out.stat().st_size == 0
    rows_written = 0

    with open(args.out, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            writer.writeheader()

        for virus_dir in sorted(p for p in sub_dir.iterdir() if p.is_dir()):
            fasta = first_fasta(virus_dir)
            if fasta is None:
                print(f"WARNING: no fasta in {virus_dir}")
                continue
            record = next(SeqIO.parse(fasta, "fasta"), None)
            if record is None:
                print(f"WARNING: empty fasta {fasta}")
                continue

            genome = str(record.seq).upper()
            try:
                hits = top_hits(aligner, genome, QUERY)
            except Exception as e:  # e.g. characters outside the matrix alphabet
                print(f"WARNING: {virus_dir.name} failed: {e}")
                continue

            row = {"matrix_class": args.matrix, "genome_length": len(genome)}
            for i, (ori, score, pos) in enumerate(hits, start=1):
                row[f"ori_{i}"], row[f"score_{i}"], row[f"position_{i}"] = ori, score, pos
            writer.writerow(row)  # missing hits (fewer than 3) are left blank
            rows_written += 1
            f.flush()

    print(f"Done: {rows_written} rows appended to {args.out}")

# can only be executed as a standalone program
if __name__ == "__main__":
    main()