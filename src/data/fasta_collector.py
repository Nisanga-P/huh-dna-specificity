from Bio import SeqIO
from pathlib import Path

class fasta_collector():

    def __init__(self):
        self.virus_rep_details = {}
        self.virus_genome_details = {}

    def virus_genome_folder_search(self):

        root_dir = Path("virus_genomes")

        for path in root_dir.iterdir():
            if path.is_dir():
                split_folder_name = path.name.split("_")

                taxid = split_folder_name[0]
                organism_name = " ".join(split_folder_name[1:])

                files_gen = (p for p in Path(path).iterdir() if p.is_file())
    
                # Return the first file path, or None if the directory has no files
                first_ref_genome_saved = next(files_gen, None)

                # double check if file
                if first_ref_genome_saved.is_file():
                    with open(first_ref_genome_saved, "r") as f:
                        
                        for record in SeqIO.parse(f, "fasta"):

                            ref_genome_str = str(record.seq)

                            self.virus_genome_details[taxid] = (organism_name, ref_genome_str)



    def parse_virus_rep_fasta_file(self, fasta_file):
    # Description in each fasta file should follow specific format:
    # > organism_taxid, organism_name

    # In the pipeline, this applies to 'virus_protein_sequences.fasta' (viral Rep protein sequences)

        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                split_description = record.description.split(" ")

                taxid = split_description[0]
                organism_name = " ".join(split_description[1:]) # join list of str with spaces
                sequence_str = str(record.seq)

                self.virus_rep_details[taxid] = (organism_name, sequence_str)

        return


parser = fasta_collector()
#parser.parse_custom_fasta_file("virus_protein_sequences.fasta")
parser.virus_genome_folder_search()
print(len(parser.virus_genome_details))