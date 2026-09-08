from Bio import Entrez
from Bio import SeqIO

import time

class protein_hit_inquirer:
    def __init__(self, iteration_name):

        """
        Constructor for protein hit inqurier.

        Args:
         - accession_ids_txt (txt file): name of SINGLE file in data/generated/psiblast_hits that follows the naming convention "{ref_rep}_iter{N}.txt"

        """

        ids_file_path = "data/raw/generated/psiblast_hits/" + iteration_name + ".txt"

        try:
            with open(ids_file_path) as f:
                accession_ids = [line.strip() for line in f if line.strip()]
        except FileNotFoundError:
            print("ERROR: No matching file in 'psiblast_hits' folder")

        self.ids_file_name = iteration_name

        self.accession_ids = accession_ids

        # for later use in accumulating organism_name, taxid (as Entrez is used to obtain virus details)
        # dict: {accession_id : (organism_name, taxid)}
        self.organisms_info = {}

        
    #PRIVATE METHOD!
    def get_organism_taxid_and_name(self):


        """
        Helper function for 'download_amino_acid_sequences'.

        Uses Entrez to obtain virus (organism) species records of the rep proteins submitted.

        Saves a file in 'virus_species_info.txt' with every organism_name, taxid in separate line (matching order of accession_ids).
        Uses the same name as the 'psiblast_hits' file.


        """

        # Give Entrez user email (good practice so warned about data usage limits)
        Entrez.email = "nisangapathirana1@gmail.com"

        # .efetch retrieves full data records from NCBI blast using accession number

        # batch_size + time.sleep used to stay within usage limits. Just a good practice -> this will not pose an issue here as we are separating 
        # accession_ids into small files.
        batch_size = 200

        for i in range(0, len(self.accession_ids), batch_size):
            batch_ids = self.accession_ids[i:i + batch_size]

            with Entrez.efetch(
                db="protein", id=",".join(batch_ids), rettype="gb", retmode="text"
            ) as handle:
                for record in SeqIO.parse(handle, "genbank"):

                    # grab virus name
                    organism_name = record.annotations.get("organism")

                    # grab virus taxid
                    for feature in record.features:
                        if feature.type == "source":
                            for xref in feature.qualifiers.get("db_xref", []):
                                if xref.startswith("taxon:"):
                                    taxid = xref.split(":")[1]


                    self.organisms_info[record.id] = (organism_name, taxid)

            time.sleep(0.34) # rate limit between batches


        new_virus_info_file_path = "data/raw/generated/virus_info/" + self.ids_file_name + "_virus_info.txt"

        with open(new_virus_info_file_path, "w", encoding="utf-8") as file:
            for id in self.accession_ids:

                taxid, organism_name = self.organisms_info[id]

                file.write(str(taxid) + "," + str(organism_name) + "\n")
            



    def download_amino_acid_sequences(self):
        """
        Downloads FASTA amino acid sequences of all protein accession ids in the accession_id txt file PROVIDED that the virus was not
        already encountered previously in the dataset (checked using current 'all_viruses_in_dataset.txt').

        Each FASTA sequence follows a custom description header of > organism_name taxid

        If a repeat is found, details are saved in 'repeats.txt' for future reference (indicator to assess if there is good coverage of Rep proteins across phylogenetic tree).
        Updates 'all_viruses_in_dataset.txt' with new additions to the dataset.

        Note: 'get_organism_taxid_and_name' should run prior to this function (if not, it will call it itself and still successfully run
        but may create conflicts with bash script downloading virus genomes as it updates 'all_viruses_in_dataset.txt')
        
        """

        # hashset containing taxid of viruses already in the dataset -> for fast lookup in checking repeats later
        all_viruses_taxids = set()

        repeat_virus_info = []

        new_viruses_taxid = []

        # load existing taxids into the set
        with open("data/raw/generated/all_viruses_in_dataset.txt", "r", encoding="utf-8") as file:
            for line in file:
                all_viruses_taxids.add(int(line.strip()))


        if not self.organisms_info:
            # ensure it was run first to save virus details into 'self.organisms_info'
            # note in the actual pipeline -> ensure that this was run separately prior (amino acids sequences should be downloaded AFTER the virus_genomes)
            self.get_organism_taxid_and_name()

        Entrez.email = "nisangapathirana1@gmail.com"

        destination = "data/raw/downloads/rep_protein_sequences/" + self.ids_file_name + ".fasta"

        with open(destination, "a", encoding="utf-8") as file:

            # as above, just as a good practice
            batch_size = 200

            for i in range(0, len(self.accession_ids), batch_size):
                batch_ids = self.accession_ids[i:i + batch_size]

                with Entrez.efetch(
                    db="protein", id=",".join(batch_ids), rettype="fasta", retmode="text"
                ) as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        organism_name, taxid = self.organisms_info[record.id]

                        if int(taxid) in all_viruses_taxids:
                            print(f"{organism_name} + Rep details already in database! Save details in 'repeats.txt'. Skipping...")
                            repeat_virus_info.append((organism_name, taxid))

                        else:
                            file.write(">" + str(taxid) + " " + str(organism_name) + "\n")
                            file.write(str(record.seq) + "\n")
                            all_viruses_taxids.add(int(taxid))
                            new_viruses_taxid.append(int(taxid))

                time.sleep(0.34) # rate limit between batches


        # if repeats were identified -> record their details in 'repeats.txt' (this is useful to gauge whether dataset has reached a good coverage of the evolutionary tree)
        if repeat_virus_info:
            with open("data/raw/generated/repeats.txt", "a", encoding="utf-8") as file:
                for name, taxid in repeat_virus_info:
                    file.write(str(name) + "," + str(taxid) + "," + str(self.ids_file_name) + "\n")


        # append the new taxid to 'all_viruses_in_dataset' (record of all viruses in our dataset for future references)
        with open("data/raw/generated/all_viruses_in_dataset.txt", "a", encoding="utf-8") as file:
            for id in new_viruses_taxid:
                file.write(str(id) + "\n")



# inq = protein_hit_inquirer("WDV_iter1.txt")

# print(inq.ids_file_name)

# inq.download_amino_acid_sequences()
