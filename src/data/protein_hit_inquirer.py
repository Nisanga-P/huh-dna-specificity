from Bio import Entrez
from Bio import SeqIO

class protein_hit_inquirer:
    def __init__(self, accession_ids_txt):

        with open(accession_ids_txt) as f:
            accession_ids = [line.strip() for line in f if line.strip()]

        self.accession_ids = accession_ids
        
        self.organisms_info = []


    def get_organism_taxid_and_name(self):
        Entrez.email = "nisangapathirana1@gmail.com"

        # # .efetch retrieves full data records from NCBI blast using accession number
        # handle = Entrez.efetch(db = "protein", id=",".join(self.accession_ids), rettype="gb", retmode="text")
        # records = SeqIO.parse(handle, "genbank") #network handle object -> Entrez.read()
        # handle.close()

        for id in self.accession_ids:
            with Entrez.efetch(
                db="protein", id=id, rettype="gb", retmode="text"
                ) as handle:
                    record = SeqIO.read(handle, "genbank")

                    organism_name = record.annotations.get("organism")

                    for feature in record.features:
                        if feature.type == "source":
                            for xref in feature.qualifiers.get("db_xref", []):
                                if xref.startswith("taxon:"):
                                    taxid = xref.split(":")[1]

                    self.organisms_info.append((organism_name, taxid))


        with open("virus_species_info.txt", "w", encoding="utf-8") as file:
            for organism_info in self.organisms_info:
                organism_name, taxid = organism_info
                file.write(str(taxid) + "," + str(organism_name) + "\n")

    def download_amino_acid_sequences(self):
        Entrez.email = "nisangapathirana1@gmail.com"

        cnt = 0 #to match it will info from self.organisms_info
        fasta_sequences = []

        for id in self.accession_ids:
                with Entrez.efetch(
                    db="protein", id=id, rettype="fasta", retmode="text"
                        ) as handle:
                        record = SeqIO.read(handle, "fasta")
                        fasta_sequences.append(str(record.seq))

        with open("virus_protein_sequences.fasta", "w", encoding="utf-8") as file:
            if not self.organisms_info:
                print("'get_organism_taxid_and_name' has not been run. Using primitive numbering system. To " \
                        "to format FASTA file with custom description >taxid,organism_name run this method first.")
                
                for seq in fasta_sequences:
                    file.write("> " + str(cnt) + "\n")
                    file.write(seq + "\n")
                    cnt += 1
            else:
                print("organism_info found. Formatting FASTA file with custom description "
                                ">taxid,organism_name ")

                for seq in fasta_sequences:
                    organism_name, taxid = self.organisms_info[cnt]
                    file.write(">" + str(taxid) + " " + str(organism_name) + "\n")
                    file.write(seq + "\n")
                    cnt += 1
                    



inq = protein_hit_inquirer("text.txt")

inq.get_organism_taxid_and_name()

inq.download_amino_acid_sequences()



# write similar function to obtain protein fasta file using accession id, and write 
# exactly according to fasta file parser format

# configure virus genome downloader to match this