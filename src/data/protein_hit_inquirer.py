from Bio import Entrez
from Bio import SeqIO

import requests

import time

class protein_hit_inquirer:
    def __init__(self, accession_ids_txt):

        with open(accession_ids_txt) as f:
            accession_ids = [line.strip() for line in f if line.strip()]

        self.accession_ids = accession_ids
        
   #     self.organisms_info = []

        # dict: {accession_id : (organism_name, taxid, uniprot_id)}
        self.organisms_info = {}

        
    #private method
    def _get_organism_taxid_and_name(self):
        # uses Entrez to obtain virus species records of protein accession ids in txt file
        # (i.e. virus Rep protein -> organism of origin)
        # saves in 'virus_species_info.txt' with every organism_name, taxid in separate line (matching order of accession_ids)
        Entrez.email = "nisangapathirana1@gmail.com"

        # # .efetch retrieves full data records from NCBI blast using accession number

        batch_size = 200

        for i in range(0, len(self.accession_ids), batch_size):
            batch_ids = self.accession_ids[i:i + batch_size]

            with Entrez.efetch(
                db="protein", id=",".join(batch_ids), rettype="gb", retmode="text"
            ) as handle:
                for record in SeqIO.parse(handle, "genbank"):

                    organism_name = record.annotations.get("organism")

                    for feature in record.features:
                        if feature.type == "source":
                            for xref in feature.qualifiers.get("db_xref", []):
                                if xref.startswith("taxon:"):
                                    taxid = xref.split(":")[1]

                    accession = record.id.split(".")[0]  # e.g. ABC123.1 -> ABC123
                    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

                    response = requests.get(url, timeout=20)

                    if response.ok: 
                        uniprot_id = data.get("uniProtkbId")
                    else:
                        print("failed")
                        uniprot_id = ""

                    self.organisms_info[record.id] = (organism_name, taxid, uniprot_id)

            time.sleep(0.34) # rate limit between batches


        # request = IdMappingClient.submit(
        #     source = "UniProtKB_AC-ID",
        #     dest="UniProtKB",
        #     ids=self.accession_ids
        # )

        # time.sleep(10.0)

        # results = list(request.each_result())

        # for result in results:
        #     name, taxid = self.organisms_info[result['from']]

        #     self.organisms_info[result['from']] = (name, taxid, result['to'])


        # with open("virus_species_info.txt", "w", encoding="utf-8") as file:
        #     # for organism_info in self.organisms_info:
        #     #     organism_name, taxid = organism_info
        #     #     file.write(str(taxid) + "," + str(organism_name) + "\n")

        #     # here id -> accession id of protein
        #     for id in self.accession_ids:

        #         with Entrez.efetch(
        #             db="protein", id=id, rettype="gb", retmode="text"
        #             ) as handle:
        #                 record = SeqIO.read(handle, "genbank")

        #                 organism_name = record.annotations.get("organism")

        #                 for feature in record.features:
        #                     if feature.type == "source":
        #                         for xref in feature.qualifiers.get("db_xref", []):
        #                             if xref.startswith("taxon:"):
        #                                 taxid = xref.split(":")[1]

        #                 # add info into self.organisms only if it is new
        #                 if taxid in self.organisms_info:
        #                     print(f"Info of organism with {taxid} already exists! Skipping...")
        #                     self.repeat_accession_ids.add(id)
        #                     continue
        #                 else:
        #                 #self.organisms_info.append((organism_name, taxid))
        #                     self.organisms_info[id] = (organism_name, taxid)    
        #                     file.write(str(taxid) + "," + str(organism_name) + "\n")

                

        # with open("virus_species_info.txt", "w", encoding="utf-8") as file:
        #     for organism_info in self.organisms_info:
        #         organism_name, taxid = organism_info
        #         file.write(str(taxid) + "," + str(organism_name) + "\n")


    def download_amino_acid_sequences(self, save_file_name):
        # downloads FASTA amino acid sequences of all protein accession ids in the txt file.
        # saved in 'virus_proteins_sequences.fasta' following custom description > organism_name taxid

        # run first to save organism details into 'self.organisms_info'
        self._get_organism_taxid_and_name()

        Entrez.email = "nisangapathirana1@gmail.com"

        destination = save_file_name + ".fasta"

        with open(destination, "w", encoding="utf-8") as file:
            print("Formatting FASTA file with custom description >taxid,organism_name,uniprot_id ")

            batch_size = 200

            for i in range(0, len(self.accession_ids), batch_size):
                batch_ids = self.accession_ids[i:i + batch_size]

                with Entrez.efetch(
                    db="protein", id=",".join(batch_ids), rettype="fasta", retmode="text"
                ) as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        organism_name, taxid, uniprotid = self.organisms_info[record.id]
                        file.write(">" + str(taxid) + " " + str(organism_name) + " " + str(uniprotid) + "\n")
                        file.write(str(record.seq) + "\n")

                time.sleep(0.34) # rate limit between batches



        # cnt = 0 #to match it with correct info tuple from self.organisms_info
        # # fasta_sequences = []

        # destination = save_file_name + ".fasta"

        # with open(destination, "w", encoding="utf-8") as file:
        #     print("Formatting FASTA file with custom description >taxid,organism_name ")

        #     # here id -> accession id of protein
        #     for id in self.accession_ids:
        #         if id in self.repeat_accession_ids:
        #             continue
        #         else:
        #             with Entrez.efetch(
        #                 db="protein", id=id, rettype="fasta", retmode="text"
        #                     ) as handle:
        #                     record = SeqIO.read(handle, "fasta")
        #                     #fasta_sequences.append(str(record.seq))
        #                     organism_name, taxid = self.organisms_info[id]
        #                     file.write(">" + str(taxid) + " " + str(organism_name) + "\n")
        #                     file.write(str(record.seq) + "\n")



        # with open(destination, "w", encoding="utf-8") as file:
        #     print("Formatting FASTA file with custom description >taxid,organism_name ")

        #     for seq in fasta_sequences:
        #         organism_name, taxid = self.organisms_info[id]
        #         print(organism_name)
        #         print(taxid)
        #         print(seq)
        #         file.write(">" + str(taxid) + " " + str(organism_name) + "\n")
        #         file.write(seq + "\n")
        #         cnt += 1
            

        # with open(destination, "w", encoding="utf-8") as file:
        #     if not self.organisms_info:
        #         print("'get_organism_taxid_and_name' has not been run. Using primitive numbering system. To " \
        #                 "to format FASTA file with custom description >taxid,organism_name run this method first.")
                
        #         for seq in fasta_sequences:
        #             file.write("> " + str(cnt) + "\n")
        #             file.write(seq + "\n")
        #             cnt += 1
        #     else:
        #         print("organism_info found. Formatting FASTA file with custom description "
        #                         ">taxid,organism_name ")

        #         for seq in fasta_sequences:
        #             organism_name, taxid = self.organisms_info[cnt]
        #             file.write(">" + str(taxid) + " " + str(organism_name) + "\n")
        #             file.write(seq + "\n")
        #             cnt += 1
                    



inq = protein_hit_inquirer("psiblast_hits/WDV_psiblast_hits_iteration1.txt")

inq.download_amino_acid_sequences("test3")



# write similar function to obtain protein fasta file using accession id, and write 
# exactly according to fasta file parser format

# configure virus genome downloader to match this