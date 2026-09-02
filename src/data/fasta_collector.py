from Bio import SeqIO
from pathlib import Path

class fasta_collector:

    def __init__(self, iteration_name):

        self.rep_sequences_file_path = "data/raw/downloads/rep_protein_sequences/" + iteration_name + ".fasta"
        self.virus_genome_folder_path = "data/raw/downloads/virus_genomes/" + iteration_name

        self.virus_rep_details = {}
        self.virus_genome_details = {}

        self.preliminary_rep_virus_dict = {}

    # PRIVATE METHOD!
    def _virus_genome_folder_search(self):

        '''
        Searches through 'downloads/virus_genomes' for virus genomes searched in the current iteration, and saves them as a dict
        in 'self.virus_genome_details' for the final 'self.preliminary_rep_virus_dict' dict.

        Helper function for 'gather_all_fasta_records'
        
        '''

        # run in try_catch block in case 'iteration_name' is not valid (i.e. no virus details saved under such an iteration)
        try:
            root_dir = Path(self.virus_genome_folder_path)

            for path in root_dir.iterdir():
                if path.is_dir():
                    split_folder_name = path.name.split("_")

                    # extract taxid, organism_name info from folder name
                    taxid = int(split_folder_name[0])
                    organism_name = " ".join(split_folder_name[1:])

                    files_gen = (p for p in Path(path).iterdir() if p.is_file())
        
                    # Return the first file path (i.e. FIRST reference genome if multiple), or None if the directory has no files
                    first_ref_genome_saved = next(files_gen, None)

                    # double check if reference genome has been found
                    if first_ref_genome_saved.is_file():
                        with open(first_ref_genome_saved, "r") as f:
                            
                            for record in SeqIO.parse(f, "fasta"):

                                ref_genome_str = str(record.seq)

                                self.virus_genome_details[taxid] = (organism_name, ref_genome_str)
                    else:
                        # no reference genome has been downloaded 
                        # should be taken care of later
                        self.virus_genome_details[taxid] = (organism_name, None)
        except FileNotFoundError as e:
            print("ERROR: fasta_collector cannot find the virus details for iteration! Double check 'iteration_name' or run methods in 'protein_hit_inquirer.py' first.")



    def _parse_virus_rep_fasta_file(self):
    # Description in each fasta file should follow specific format:
    # > organism_taxid, organism_name

    # In the pipeline, this applies to 'virus_protein_sequences.fasta' (viral Rep protein sequences)

        with open(self.rep_sequences_file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                split_description = record.description.split(" ")

                taxid = int(split_description[0])
                organism_name = " ".join(split_description[1:]) # join list of str with spaces
                sequence_str = str(record.seq)

                self.virus_rep_details[taxid] = (organism_name, sequence_str)


    def get_sequence_collection(self):
        '''
        Searches through the (protein, DNA) fasta files generated upto this point in the pipeline,
        and builds a preliminary dictionary for the dataset in the form:

        {taxid : (virus_name, rep_seq, genome_seq)}
        '''

        # Run helper methods to search through the relevant fasta files, and initialise the required details dicts
        self._virus_genome_folder_search()
        self._parse_virus_rep_fasta_file()

        for taxid in self.virus_rep_details:

            virus_name, rep_seq = self.virus_rep_details.get(taxid)

            _, genome_seq = self.virus_genome_details.get(taxid)

            self.preliminary_rep_virus_dict[taxid] = (virus_name, rep_seq, genome_seq)

        return self.preliminary_rep_virus_dict



parser = fasta_collector("WDV_iter1")
#parser.parse_custom_fasta_file("virus_protein_sequences.fasta")
parser.get_sequence_collection()
print(len(parser.preliminary_rep_virus_dict))
print(parser.preliminary_rep_virus_dict.get(223275))
