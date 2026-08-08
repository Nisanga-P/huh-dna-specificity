import Bio
from Bio import SeqIO
from Bio import Align

from pathlib import Path


muscovy_genome = Path("muscovy.txt").read_text(encoding="utf-8")


aligner = Align.PairwiseAligner()

aligner.mode = "local" #configure for local alignment (Smith-Waterman)

ori_given = "TATTATTAC"

score = aligner.score(muscovy_genome, ori_given)
alignments = aligner.align(muscovy_genome, ori_given)

print(f"Optimal Local Alignment Score: {score}\n")

for alignment in alignments:
    print(alignment)
