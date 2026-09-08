import Bio
from Bio import Align
from Bio import SeqIO
from Bio.Align.substitution_matrices import Array

from fasta_collector import fasta_collector

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
        self.ori_query = "HMDWNYTAC" 



    # HELPER method!
    def _local_alignment(self, virus_genome_str):
        # Returns result (from local alignment), original virus genome: str, optimal score

        # perform local alignemnt
        optimal_score = self.aligner.score(virus_genome_str, self.ori_query)
        alignments = self.aligner.align(virus_genome_str, self.ori_query) 
        #only alignments with optimal score saved in 'alignments' -> adjust substitution matrix to give just one alignment match (i.e. most likely ori)

        #Optimal alignment score -> confidence of suggested ori sequence 
        return alignments, optimal_score


    def obtain_ori_candidates(self, virus_genome_str):
        '''

        Once local alignment has been performed, converts results into list of ori candidates.

        Args:
           alignments: Pairwise alignment object
               -> result of local alignment
           virus_genome_str: str
               -> virus genome as a string (not fasta!)

         Returns:
           Ori candidates: [str]
        
        '''

        alignments, optimal_score = self._local_alignment(virus_genome_str)

        ori_candidates = []

        for alignment in alignments:
            ori_start_index = alignment.coordinates[0][0]
            ori_end_index = alignment.coordinates[0][1]

            ori_candidates.append(virus_genome_str[ori_start_index:ori_end_index])

        return ori_candidates, optimal_score


# The phylum of cressdnaviricota split into two classes. These are treated separately using child classes.

class repensiviricetes_ori_finder(ori_finder):
    # including geminivirus family
    def __init__(self):
        super().__init__()

        bases = "ATGC"
        alphabet = "ACGTRYSWKMHDN" #include "" or "" (edge case where download ref genome has ambiguous base code)


        #---------- Custom Substitution Matrix ----------

        # # Row (Target Base) x Column (Query Base)

        default_scores = np.full((13, 13), -np.inf)
        sub_matrix = Array(alphabet=alphabet, dims=2, data=default_scores)


        # A(+1), C(+2): Final 2 bases sequences in 9-mer ori MUST be exactly "AC" 

        sub_matrix["A", "A"] = 50.0 # must be a large value (avoid edge case where terminating sequence gives optimal alignment -> "AC" never matched so -inf penalty not applied)
        sub_matrix["C", "C"] = 50.0

        #

        # T(-1):
        sub_matrix.update({('A', 'T'): 0.9, ('C', 'T'): 1.0, ('G', 'T'): 0.8, ('T', 'T'): 6.8})

        # Y(-2):
        sub_matrix.update({('A', 'Y'): 1.4, ('C', 'Y'): 2.8, ('G', 'Y'): 0.8, ('T', 'Y'): 4.1})

        # N(-3):
        sub_matrix.update({('A', 'N'): 3.6, ('C', 'N'): 1.8, ('G', 'N'): 1.9, ('T', 'N'): 2.3})

        # W(-4):
        sub_matrix.update({('A', 'W'): 1.4, ('C', 'W'): 0.9, ('G', 'W'): 0.4, ('T', 'W'): 6.6})

        # D(-5):
        sub_matrix.update({('A', 'D'): 6.6, ('C', 'D'): 0.8, ('G', 'D'): 0.7, ('T', 'D'): 1.3})

        # M(-6):
        sub_matrix.update({('A', 'M'): 6.9, ('C', 'M'): 0.9, ('G', 'M'): 0.3, ('T', 'M'): 1.3})

        # H(-7):
        sub_matrix.update({('A', 'H'): 2.0, ('C', 'H'): 2.1, ('G', 'H'): 0.6, ('T', 'H'): 4.7})


        self.aligner.substitution_matrix = sub_matrix




class arfiviricetes_ori_finder(ori_finder):
    # including circovirus and nanovirus families
    def __init__(self):
        super().__init__()

        bases = "ATGC"
        alphabet = "ACGTRYSWKMHDN" #include "" or "" (edge case where download ref genome has ambiguous base code)


        #---------- Custom Substitution Matrix ----------

        # # Row (Target Base) x Column (Query Base)

        default_scores = np.full((13, 13), -np.inf)
        sub_matrix = Array(alphabet=alphabet, dims=2, data=default_scores)


        # A(+1), C(+2): Final 2 bases sequences in 9-mer ori MUST be exactly "AC" 

        sub_matrix["A", "A"] = 50.0 # must be a large value (avoid edge case where terminating sequence gives optimal alignment -> "AC" never matched so -inf penalty not applied)
        sub_matrix["C", "C"] = 50.0

        #

        # T(-1):
        sub_matrix.update({('A', 'T'): 0.6, ('C', 'T'): 0.5, ('G', 'T'): 0.0, ('T', 'T'): 8.6})

        # Y(-2):
        sub_matrix.update({('A', 'Y'): 1.2, ('C', 'Y'): 2.7, ('G', 'Y'): 0.7, ('T', 'Y'): 4.7})

        # N(-3):
        sub_matrix.update({('A', 'N'): 4.3, ('C', 'N'): 1.7, ('G', 'N'): 1.7, ('T', 'N'): 2.2})

        # W(-4):
        sub_matrix.update({('A', 'W'): 0.9, ('C', 'W'): 0.3, ('G', 'W'): 0.1, ('T', 'W'): 8.4})

        # D(-5):
        sub_matrix.update({('A', 'D'): 1.7, ('C', 'D'): 1.6, ('G', 'D'): 3.3, ('T', 'D'): 3.0})

        # M(-6):
        sub_matrix.update({('A', 'M'): 5.9, ('C', 'M'): 1.1, ('G', 'M'): 0.9, ('T', 'M'): 1.6})

        # H(-7):
        sub_matrix.update({('A', 'H'): 2.6, ('C', 'H'): 2.5, ('G', 'H'): 1.1, ('T', 'H'): 3.3})


        self.aligner.substitution_matrix = sub_matrix




# child = repensiviricetes_ori_finder()

# child2 = arfiviricetes_ori_finder()



# parser = fasta_collector("WDV_iter1")
# #parser.parse_custom_fasta_file("virus_protein_sequences.fasta")
# prelim = parser.get_sequence_collection()
# print(len(prelim))


# ori_find = ori_finder()

# for taxid in prelim:

#     # if taxid == 1367126:
#     #     continue

#     name, _, genome_str = prelim.get(taxid)

#     # print(taxid)
#     # print(name)
#     # # print(genome_str)

#     # alignments, score = ori_find.local_alignment(genome_str)


#     # print(ori_find.obtain_ori_candidates(alignments, genome_str))

#     # print(score)


#     # print("-----")


#     print(taxid)
#     print(name)
#     # print(genome_str)

#     alignments2, score2 = child.local_alignment(genome_str)


#     print(child.obtain_ori_candidates(alignments2, genome_str))

#     print(score2)

#     print("----------------------------------------------------------------------")



# _, _, tr = prelim.get(10834)

# print(tr[2743:2752])




# _, _, ws = prelim.get(1367126)

# for char in ws:
#     if not (char == 'A' or char == 'T' or char == 'C' or char == 'G'):
#         print(char)




#         #---------- Custom Substitution Matrix ----------

#         # # Row (Target Base) x Column (Query Base)
#         # sub_matrix = substitution_matrices.load("NUC.4.4") -> other option to load existing substitution matrix and editing it


#         # Initialise custom matrix (use convention of -4.0 to 5.0 scores as in "NUC.4.4")
#         bases = "ATGC"
#         alphabet = "ATGCWYHDN"
#         default_scores = np.full((9, 9), -4.0)

        
#         sub_matrix = Array(alphabet=alphabet, dims=2, data=default_scores)


#         # A(+1), C(+2): Final 2 bases sequences in 9-mer ori MUST be exactly "AC" 

#         for residue in bases:
#             sub_matrix[residue, "A"] = -math.inf
#             sub_matrix[residue, "C"] = -math.inf

#         sub_matrix["A", "A"] = 50.0 # must be a large value (avoid edge case where terminating sequence gives optimal alignment -> "AC" never matched so -inf penalty not applied)
#         sub_matrix["C", "C"] = 50.0

#         # T(-4), T(-1): Most likely definitely "T", but in a few gemini-viruses not specific (can be any base)
#         sub_matrix[:, "T"] = -1.5
#         sub_matrix["T", "T"] = 5.0

#         # Y(-2): More likely to be "T" or "C" but not specific in general
#         sub_matrix["T", "Y"] = 4.0
#         sub_matrix["C", "Y"] = 3.0
#         sub_matrix["A", "Y"] = 0.0
#         sub_matrix["G", "Y"] = 0.0

#         # N(-3): Not specific (entirely random)

#         sub_matrix[:, "N"] = 0.0

#         # D(-5): For specific viruses, almost always "A" but not specific for others.

#         sub_matrix["T", "D"] = 1.0
#         sub_matrix["C", "D"] = 0.0
#         sub_matrix["A", "D"] = 2.5
#         sub_matrix["G", "D"] = 1.0

#         # W(-6): More likely to be "A", small chance of "T" but "C" and "G" unlikely
#         sub_matrix["T", "W"] = 1.0
#         sub_matrix["C", "W"] = -1.0
#         sub_matrix["A", "W"] = 3.0
#         sub_matrix["G", "W"] = -1.0

#         # H(-7): More likely to be "T", equal chance of "A" or "C", "G" quite unlikely
#         sub_matrix["T", "H"] = 3.0
#         sub_matrix["C", "H"] = 1.5
#         sub_matrix["A", "H"] = 1.5
#         sub_matrix["G", "H"] = -2.0


#         # Target genome will only have bases "A", "T", "G", "C" (make all rows other than first 4 NaN)
#         for residue in alphabet:
#             if residue not in bases:
#                 sub_matrix[residue, :] = None
        
#         sub_matrix[:, "G"] = None

#         self.aligner.substitution_matrix = sub_matrix
        
#         print(sub_matrix)
