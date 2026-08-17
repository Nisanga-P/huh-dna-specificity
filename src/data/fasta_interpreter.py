from Bio import SeqIO

class fasta_interpreter():

    # Description in each fasta file should follow specific format:
    # > organism_taxid, organism_name

    # In the pipeline, this applies to 'virus_protein_sequences.fasta' (viral Rep protein sequences)

    def parse_custom_fasta_file(self, fasta_file):
        sequences_str = []
        info_str = []

        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                sequences_str.append(str(record.seq))

                split_description = record.description.split(" ")

                info_str.append((split_description[0], split_description[1], str(record.seq)))

        return sequences_str, info_str


parser = fasta_interpreter()
a, b = parser.parse_custom_fasta_file("virus_protein_sequences.fasta")
print(a)
print(b)