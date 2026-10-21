# from pathlib import Path
# from Bio import SeqIO

# class msa_prepper:
#     def __init__(self):

#         self.msa_folder_path = "data/processed/MSA_fasta/"


#     def prepare_viral_genome_msa(self, msa_file_name, skip=1):

#         counter = 0

#         msa_file_path = self.msa_folder_path + msa_file_name + ".fasta"

#         virus_genomes_dir = Path("data/raw/downloads/virus_genomes/")

#         with open(msa_file_path, "w", encoding="utf-8") as file:

#             for subfolders in virus_genomes_dir.iterdir():
#                 for virus_folder in subfolders.iterdir():

#                     if counter % skip != 0:
#                         counter += 1
#                         continue

#                     split_folder_name = virus_folder.name.split("_")

#                     header = ">" + split_folder_name[1] + "_" + split_folder_name[2] + "_" + split_folder_name[0]

#                     files_gen = (p for p in Path(virus_folder).iterdir() if p.is_file())
    
#                     # Return the first file path (i.e. FIRST reference genome if multiple), or None if the directory has no files
#                     first_ref_genome_saved = next(files_gen, None)

#                     # double check if reference genome has been found
#                     if first_ref_genome_saved.is_file():

#                         if first_ref_genome_saved is not None:
#                             # Read first FASTA record
#                             record = next(
#                                 SeqIO.parse(first_ref_genome_saved, "fasta"),
#                                 None
#                             )

#                             if record is not None:

#                                 file.write(header + "\n")
#                                 file.write(str(record.seq) + "\n")

#                                 counter += 1

from pathlib import Path
from Bio import SeqIO


class msa_prepper:
    def __init__(self):
        self.msa_folder_path = Path("data/processed/MSA_fasta/")

    def prepare_viral_genome_msa(
        self,
        msa_file_name,
        folder_ending="A",
        skip=1,
    ):
        if folder_ending not in {"A", "R"}:
            raise ValueError("folder_ending must be 'A' or 'R'")
        if skip < 1:
            raise ValueError("skip must be at least 1")

        virus_genomes_dir = Path("data/raw/downloads/virus_genomes/")
        output_path = self.msa_folder_path / f"{msa_file_name}.fasta"
        self.msa_folder_path.mkdir(parents=True, exist_ok=True)

        counter = 0

        with output_path.open("w", encoding="utf-8") as file:
            for parent_folder in virus_genomes_dir.iterdir():
                # Only process parent folders whose names end in A or R
                if not parent_folder.is_dir() or not parent_folder.name.endswith(folder_ending):
                    continue

                for virus_folder in parent_folder.iterdir():
                    if not virus_folder.is_dir():
                        continue

                    # Apply skip to genome folders within the matching parent folders
                    if counter % skip != 0:
                        counter += 1
                        continue
                    counter += 1

                    folder_parts = virus_folder.name.split("_")
                    if len(folder_parts) < 3:
                        continue

                    header = f">{folder_parts[1]}_{folder_parts[2]}_{folder_parts[0]}"

                    fasta_files = sorted(
                        path for path in virus_folder.iterdir()
                        if path.is_file()
                        and path.suffix.lower() in {".fasta", ".fa", ".fna"}
                    )

                    if not fasta_files:
                        continue

                    # Read the first record from the first FASTA file
                    record = next(SeqIO.parse(fasta_files[0], "fasta"), None)
                    if record is not None:
                        file.write(f"{header}\n{record.seq}\n")
                        
                                

msa_prepper = msa_prepper()

msa_prepper.prepare_viral_genome_msa("genome_file", "A", 3)
