from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from Bio import SeqIO

class blastp_Searcher:

    def __init__(self, protein_seq_fasta):

     #   self.sequence = SeqIO.read(protein_seq_fasta, "fasta")
     self.sequence = protein_seq_fasta

    def run_search(self, dataset = "swissprot"):

        result_handle = NCBIWWW.qblast("blastp", dataset, self.sequence)

        with open("blast_results.xml", "w") as out_file:
            out_file.write(result_handle.read())

        result_handle.close()


    def parse_and_filter_results(self):
        
        E_VALUE_THRESH = 0.001

        
        with open("blast_results.xml") as result_file:
            
            blast_records = NCBIXML.parse(result_file)
            
            
            for record in blast_records:
                print(f"\nQuery: {record.query[:50]}...")
                
                # Iterate through every matching database sequence (alignment)
                for alignment in record.alignments:
                    # Iterate through individual High-scoring Segment Pairs (HSPs)
                    for hsp in alignment.hsps:
                        # Only show statistically significant matches
                        if hsp.expect < E_VALUE_THRESH:
                            print("\n" + "*" * 40)
                            print(f"Match Title: {alignment.title}")
                            print(f"Sequence Length: {alignment.length}")
                            print(f"E-value: {hsp.expect}")
                            print(f"Identities: {hsp.identities} / {hsp.align_length}")
                            print(f"Query Alignment:  {hsp.query[:50]}")
                            print(f"Match Alignment:  {hsp.sbjct[:50]}")


blastp_search = blastp_Searcher("tomato_protein.fa")
blastp_search.run_search()
blastp_search.parse_and_filter_results()