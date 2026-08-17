from Bio import Entrez
from Bio import SeqIO

class protein_hit_inquirer:
    def __init__(self, accession_ids_txt):

        with open(accession_ids_txt) as f:
            accession_ids = [line.strip() for line in f if line.strip()]

        self.accession_ids = accession_ids
        
        self.organisms_info = []


    def get_organism_taxid_and_name(self):
        # uses Entrez to obtain virus species records of protein accession ids in txt file
        # (i.e. virus Rep protein -> organism of origin)
        # saves in 'virus_species_info.txt' with every organism_name, taxid in separate line (matching order of accession_ids)
        Entrez.email = "nisangapathirana1@gmail.com"

        # # .efetch retrieves full data records from NCBI blast using accession number

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
        # downloads FASTA amino acid sequences of all protein accession ids in the txt file.
        # saved in 'virus_proteins_sequences.fasta' following custom description > organism_name taxid


        Entrez.email = "nisangapathirana1@gmail.com"

        cnt = 0 #to match it with correct info tuple from self.organisms_info
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
                    



inq = protein_hit_inquirer("psiblast_hits/WDV_psiblast_hits_iteration1.txt")

inq.get_organism_taxid_and_name()

inq.download_amino_acid_sequences()



# write similar function to obtain protein fasta file using accession id, and write 
# exactly according to fasta file parser format

# configure virus genome downloader to match this