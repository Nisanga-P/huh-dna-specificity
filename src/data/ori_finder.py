import Bio
from Bio import Align
from Bio import SeqIO
from Bio.Align.substitution_matrices import Array

import numpy as np

from pathlib import Path

import math


class ori_finder:
    def __init__(self):

        self.aligner = Align.PairwiseAligner()

        # ---------- Configure settings for local alignment ----------
        self.aligner.mode = "local" # local alignment (Smith-Waterman)

        # prevent aligner from adding gaps during the alignment -> the 9-mer ori query should be aligned exactly!
        self.aligner.open_gap_score = -math.inf 
        self.aligner.extend_gap_score = -math.inf
        

        # average 'ori' sequence -> forms the query sequence for local alignent
        # Deduced from weighted sequence logos in Fig 5. -> results of HUH-seq assay developed by Tompkins et al.
        # see https://www.bioinformatics.org/sms/iupac.html for IUPAC codes for bases
        self.ori_query = "HWDTNYTAC" 


        #---------- Custom Substitution Matrix ----------

        # # Row (Target Base) x Column (Query Base)
        # sub_matrix = substitution_matrices.load("NUC.4.4") -> other option to load existing substitution matrix and editing it


        # Initialise custom matrix (use convention of -4.0 to 5.0 scores as in "NUC.4.4")
        bases = "ATGC"
        alphabet = "ATGCWYHDN"
        default_scores = np.full((9, 9), -4.0)

        
        sub_matrix = Array(alphabet=alphabet, dims=2, data=default_scores)


        # A(+1), C(+2): Final 2 bases sequences in 9-mer ori MUST be exactly "AC" 

        for residue in bases:
            sub_matrix[residue, "A"] = -math.inf
            sub_matrix[residue, "C"] = -math.inf

        sub_matrix["A", "A"] = 30.0 # must be a large value (avoid edge case where terminating sequence gives optimal alignment -> "AC" never matched so -inf penalty not applied)
        sub_matrix["C", "C"] = 30.0

        # T(-4), T(-1): Most likely definitely "T", but in a few gemini-viruses not specific (can be any base)
        sub_matrix[:, "T"] = -1.5
        sub_matrix["T", "T"] = 5.0

        # Y(-2): More likely to be "T" or "C" but not specific in general
        sub_matrix["T", "Y"] = 4.0
        sub_matrix["C", "Y"] = 3.0
        sub_matrix["A", "Y"] = 0.0
        sub_matrix["G", "Y"] = 0.0

        # N(-3): Not specific (entirely random)

        sub_matrix[:, "N"] = 0.0

        # D(-5): For specific viruses, almost always "A" but not specific for others.

        sub_matrix["T", "D"] = 1.0
        sub_matrix["C", "D"] = 0.0
        sub_matrix["A", "D"] = 2.5
        sub_matrix["G", "D"] = 1.0

        # W(-6): More likely to be "A", small chance of "T" but "C" and "G" unlikely
        sub_matrix["T", "W"] = 1.0
        sub_matrix["C", "W"] = -1.0
        sub_matrix["A", "W"] = 3.0
        sub_matrix["G", "W"] = -1.0

        # H(-7): More likely to be "T", equal chance of "A" or "C", "G" quite unlikely
        sub_matrix["T", "H"] = 3.0
        sub_matrix["C", "H"] = 1.5
        sub_matrix["A", "H"] = 1.5
        sub_matrix["G", "H"] = -2.0


        # Target genome will only have bases "A", "T", "G", "C" (make all rows other than first 4 NaN)
        for residue in alphabet:
            if residue not in bases:
                sub_matrix[residue, :] = None
        
        sub_matrix[:, "G"] = None

        self.aligner.substitution_matrix = sub_matrix
        
        print(sub_matrix)


    def local_alignment(self, virus_genome_fasta):

        virus_genome_record = SeqIO.read("equine.txt", "fasta")
        #genome str should be a SINGLE line with no "\n" (newline characters)
        virus_genome_str = str(virus_genome_record.seq).replace("\n", "") 

        # perform local alignemnt
        score = self.aligner.score(virus_genome_str, self.ori_query)
        alignments = self.aligner.align(virus_genome_str, self.ori_query)

        #Optimal alignment score -> confidence of suggested ori sequence (only alignments with optimal score saved in 'alignments')

        for alignment in alignments:
            print(alignment)
            #print(alignment.aligned)
            print(alignment.coordinates)

        print(score)

virus_genome_record = SeqIO.read("equine.txt", "fasta")
virus_genome_str = str(virus_genome_record.seq).replace("\n", "") 
print(virus_genome_str)


ori_finder1 = ori_finder()
# ori_finder2 = ori_finder()


from collections import Counter

def count_characters_fast(text):
    # Pass the string directly into Counter and convert to a standard dict
    return dict(Counter(text))

print(count_characters_fast(virus_genome_str))




ori_finder1.local_alignment("equine.txt")
# ori_finder2.local_alignment("PCV2.txt")
# ori_finder2.local_alignment("WDV.txt")



#matrix = substitution_matrices.load("NUC.4.4")

# ---Edit Substition Matrix---

# +1, +2: Final 2 bases sequences in 9-mer MUST be exactly "AC" 

# matrix[:, "A"] = -math.inf
# matrix["A", "A"] = 5
# matrix[:, "C"] = -math.inf
# matrix["C", "C"] = 5

# #-2: Likely "T" or "C"

# matrix["C", "Y"] = 4
# matrix["T", "Y"] = 4
# matrix["A", "Y"] = -1
# matrix["G", "Y"] = -1


# #-5: Likely "A", "G" and possibly "T" 

# matrix["A", "D"] = 4
# matrix["G", "D"] = 3
# matrix["T", "D"] = 2.5
# matrix["C", "D"] = 0

# #-6: Likely "A", possibly "T" 

# matrix["A", "W"] = 4
# matrix["T", "W"] = 1.5
# matrix["G", "W"] = 0
# matrix["C", "W"] = 0


# #-7: "T", "A", "C" equally likely

# matrix["T", "H"] = 4
# matrix["A", "H"] = 3.5
# matrix["C", "H"] = 3.5
# matrix["G", "H"] = -1


# matrix[:, "N"] = 0


# aligner.substitution_matrix = matrix


# print(matrix)

# def local_alignment(self, genome, ori):
#     score = self.aligner.score(genome, ori)
#     alignments = self.aligner.align(genome, ori)

#     print(f"Optimal Local Alignment Score: {score}\n")

#     for alignment in alignments:
#         print(alignment)
#         score_i = 0

#         for t_char, q_char in zip(alignment.target, alignment.query):
#                 score_i += aligner.score(t_char, q_char)
#         print(f"Score: {score}")


# # print(aligner.substitution_matrix)

# # ori_given2 = "HATTNTTAC"
# ori_given2 = "HWDTNYTAC"

# local_alignment(muscovy_genome, ori_given2)




