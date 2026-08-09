import Bio
from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices

from pathlib import Path

import math


muscovy_genome = Path("muscovy.txt").read_text(encoding="utf-8").replace("\n", "")
equine_genome = Path("equine.txt").read_text(encoding="utf-8").replace("\n", "")


aligner = Align.PairwiseAligner()
aligner.open_gap_score = -math.inf
aligner.extend_gap_score = -math.inf


aligner.mode = "local" #configure for local alignment (Smith-Waterman)

ori_given = "TATTATTAC"


matrix = substitution_matrices.load("NUC.4.4")

# ---Edit Substition Matrix---

# +1, +2: Final 2 bases sequences in 9-mer MUST be exactly "AC" 

matrix[:, "A"] = -math.inf
matrix["A", "A"] = 5
matrix[:, "C"] = -math.inf
matrix["C", "C"] = 5

#-2: Likely "T" or "C"

matrix["C", "Y"] = 4
matrix["T", "Y"] = 4
matrix["A", "Y"] = -1
matrix["G", "Y"] = -1


#-5: Likely "A", "G" and possibly "T" 

matrix["A", "D"] = 4
matrix["G", "D"] = 3
matrix["T", "D"] = 2.5
matrix["C", "D"] = 0

#-6: Likely "A", possibly "T" 

matrix["A", "W"] = 4
matrix["T", "W"] = 1.5
matrix["G", "W"] = 0
matrix["C", "W"] = 0


#-7: "T", "A", "C" equally likely

matrix["T", "H"] = 4
matrix["A", "H"] = 3.5
matrix["C", "H"] = 3.5
matrix["G", "H"] = -1


matrix[:, "N"] = 0


aligner.substitution_matrix = matrix


print(matrix)

def local_alignment(genome, ori):
    score = aligner.score(genome, ori)
    alignments = aligner.align(genome, ori)

    print(f"Optimal Local Alignment Score: {score}\n")

    for alignment in alignments:
        print(alignment)
        score_i = 0

        for t_char, q_char in zip(alignment.target, alignment.query):
                score_i += aligner.score(t_char, q_char)
        print(f"Score: {score}")


# print(aligner.substitution_matrix)

# ori_given2 = "HATTNTTAC"
ori_given2 = "HWDTNYTAC"

local_alignment(muscovy_genome, ori_given2)



# from collections import Counter

# def count_characters_fast(text):
#     # Pass the string directly into Counter and convert to a standard dict
#     return dict(Counter(text))

# # Example usage:
# print(count_characters_fast(muscovy_genome.replace("\n", "")))  
# # Output: {'a': 1, 'p': 2, 'l': 1, 'e': 1}