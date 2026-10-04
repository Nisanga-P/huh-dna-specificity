from pathlib import Path
from Bio import SeqIO

class msa_prepper:
    def __init__(self):

        self.msa_folder_path = "data/processed/MSA_fasta/"


    def prepare_viral_genome_msa(self, msa_file_name, skip=1):

        counter = 0

        msa_file_path = self.msa_folder_path + msa_file_name + ".fasta"

        virus_genomes_dir = Path("data/raw/downloads/virus_genomes/")

        with open(msa_file_path, "w", encoding="utf-8") as file:

            for subfolders in virus_genomes_dir.iterdir():
                for virus_folder in subfolders.iterdir():

                    if counter % skip != 0:
                        counter += 1
                        continue

                    split_folder_name = virus_folder.name.split("_")

                    header = ">" + split_folder_name[1] + "_" + split_folder_name[2] + "_" + split_folder_name[0]

                    files_gen = (p for p in Path(virus_folder).iterdir() if p.is_file())
    
                    # Return the first file path (i.e. FIRST reference genome if multiple), or None if the directory has no files
                    first_ref_genome_saved = next(files_gen, None)

                    # double check if reference genome has been found
                    if first_ref_genome_saved.is_file():

                        if first_ref_genome_saved is not None:
                            # Read first FASTA record
                            record = next(
                                SeqIO.parse(first_ref_genome_saved, "fasta"),
                                None
                            )

                            if record is not None:

                                file.write(header + "\n")
                                file.write(str(record.seq) + "\n")

                                counter += 1


                        
                                

msa_prepper = msa_prepper()

msa_prepper.prepare_viral_genome_msa("genome_file", 4)
