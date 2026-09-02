import subprocess
import sys

from src.data.protein_hit_inquirer import protein_hit_inquirer

'''

MAIN for dataset pipeline.


Has the following methods:
 - wipe_existing_dataset
 - remove_records (?)
 - execute_pipeline -> adds to the dataset from a psi_blast hit of Rep proteins (given as a file of NCBI Accession IDs)

'''



# MAIN function!
def execute_pipeline(accession_ids_txt):

    inquirer = protein_hit_inquirer(accession_ids_txt)

    inquirer.get_organism_taxid_and_name()

    # Run 'virus_genome_downloader.sh' -> this should be run BEFORE protein sequences are downloaded! 
    # (the bash scipt uses 'all_viruses_in_dataset.txt' to check if virus is already present in the genome -> this should be updated only once it has been executed)
    process = subprocess.Popen(
        ['bash', 'src/data/virus_genome_downloader.sh', 'data/raw/generated/virus_info/WDV_iter1_virus_info.txt', "WDV_iter2"], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT,  
        text=True                 
    )

    # Prints bash output line-by-line
    for line in process.stdout:
        print(line, end="")        
        sys.stdout.flush()         

    # Wait for the script to finish
    return_code = process.wait()

    inquirer.download_amino_acid_sequences()