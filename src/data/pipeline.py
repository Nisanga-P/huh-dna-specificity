import subprocess
import sys

from pathlib import Path

import csv

import pandas as pd

from protein_hit_inquirer import protein_hit_inquirer
from fasta_collector import fasta_collector
from sDBM_finder import sDBM_finder
from ori_finder import arfiviricetes_ori_finder, repensiviricetes_ori_finder

'''

MAIN for dataset pipeline.


Has the following methods:
 - wipe_existing_dataset
 - remove_records (?)
 - execute_pipeline -> adds to the dataset from a new set of psi_blast hit of Rep proteins (given as a file of NCBI Accession IDs)

'''

def create_new_dataset(name):
    dataset_file = "data/raw/" + name + ".csv"

    columns = ["taxid", "virus_name", "sDBM_pos", "sDBM_sim_gap", "ori_seq", "ori_align_score", "iteration_name", "rep_amino_acid_seq"]

    if not Path(dataset_file).exists():
        df = pd.DataFrame(columns=columns)

        df.to_csv(dataset_file, index=False)
    else:
        print("Dataset already exists!")



def wipe_existing_dataset():
    pass



def add_records(iteration_name, reference_Rep, file_name = "data/raw/positive_rep_virus_dataset.csv"):

    '''
    Edits:
     - positive_rep_virus_dataset.csv (actual raw dataset)
     - Downloads genomes into 'virus_genomes'
     - Possibly edits 'all_viruses_in_dataset.txt' and 'repeats.txt' (depending on whether there records are new/old)
     - Downloads protein sequence in 'rep_protein_sequences'
     - Creates a text file with iteration_name in 'virus_info' folder

    '''


    # #-----------------load existing dataset------------------
    # current_df = pd.read_csv()

    #-----------------analyse psiblast hits-------------------
    inquirer = protein_hit_inquirer(iteration_name)

    # 'virus_info' folder should have a file containing virus organism info
    inquirer.get_organism_taxid_and_name()


    #----------------download virus genomes-------------------

    virus_info_path = "data/raw/generated/virus_info/" + iteration_name + "_virus_info.txt"

    # Run 'virus_genome_downloader.sh' -> this should be run BEFORE protein sequences are downloaded (.download_amino_acid_sequences)! 
    # (the bash scipt uses 'all_viruses_in_dataset.txt' to check if virus is already present in the genome -> this should be updated only after bash script has been executed)
    process = subprocess.Popen(
        # after script name, arguments that the script expects
        ['bash', 'src/data/virus_genome_downloader.sh', virus_info_path, iteration_name], 
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

    #-----------------download amino acid sequnces---------------
    inquirer.download_amino_acid_sequences()

    #-----------------collect saved fasta files and compile into single dict-----------------
    collector = fasta_collector(iteration_name)  

    prelim_virus_rep_dict = collector.get_sequence_collection()


    #------------------------identify sDBM and ori-----------------------------------------

    sDBM_query = sDBM_finder()

    pcv_seq = "MPSKKNGRSGPQPHKRWVFTLNNPSEDERKKIREPPISLFDYFIVGEEGNEEGRTPHLQGFANFVKKQTFNKVKWYLGARCHIEKAKGTDQQNKEYCSKEGNLLIECGAPRSQGQRSDLSTAVSTLLESGSLVTVAEQHPVTFVRNFRGLAELLKVSGKMQKRDWKTNVHVIVGPPGCGKSKWAANFADPETTYWKPPRNKWWDGYHGEEVVVIDDFYGWLPWDDLLRLCDRYPLTVETKGGTVPFLARSILITSNQTPLEWYSSTAVPAVEALYRRITSLVFWKNATEQSTEEGGQFVTLSPPCPEFPYEINY"
    wdv_seq = "MASSSAPRFRVYSKYLFLTYPECTLEPQYALDSLRTLLNKYEPLYIAAVRELHEDGSPHLHVLVQNKLRASITNPNALNLRMDTSPFSIFHPNIQAAKDCNQVRDYITKEVDSDVNTAEWGTFVAVSTPGRKDRDADMKQIIESSSSREEFLSMVCNRFPFEWSIRLKDFEYTARHLFPDPVATYTPEFPTESLICHETIESWKNEHLYSESPGRHKSIYICGPTRTGKTSWARSLGTHNYYNSLVDFTTYDVNAKYNIIDDIPFKFTPNWKCFVGAQRDFTVNPKYGKRKVIRGGIPCIILVNPDEDWLKDMTPEQSDYMYSNTVVHYMYEGETFINYSFASGEDVTASQ"

    if (reference_Rep == "WDV"):
        # 91 -> known starting position of WDV sDBM
        ref_sDBM = sDBM_query.get_reference_sDBM_embedding(wdv_seq, 91)

        ori_query = repensiviricetes_ori_finder()

    elif (reference_Rep == "PCV"):
        # 80 -> known starting position of PCV sDBM
        ref_sDBM = sDBM_query.get_reference_sDBM_embedding(pcv_seq, 80)

        ori_query = arfiviricetes_ori_finder()
    else:
        raise ValueError("Reference Rep protein must be either 'WDV' or 'PCV'")


    # [ [] , [] ,...] -> each list is a row
    new_records = []

    for taxid in prelim_virus_rep_dict:

        virus_name, rep_seq, genome_seq = prelim_virus_rep_dict.get(taxid)

        sDBM_pred, sDBM_pos, sDBM_confidence = sDBM_query.find_sDBM(rep_seq, ref_sDBM)

        ori_candidates, ori_score = ori_query.obtain_ori_candidates(genome_seq)

        if len(ori_candidates) > 1:
            ori_candidate = None
        else:
            ori_candidate = ori_candidates[0]

        new_records.append([taxid, virus_name, sDBM_pos, sDBM_pred, sDBM_confidence, ori_candidate, ori_score, iteration_name, rep_seq])


    if not Path(file_name).exists():
        raise FileNotFoundError("No existing csv file! Create dataset before adding records.")


    with open(file_name, mode='a', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
            
        for record in new_records:
            writer.writerow(record)
        


# MAIN function!
def execute_pipeline():

    #name convention for iteration: "{Abbrev for organism}_iter{N}_{M}_{V/B}_{A/R/None}"

    create_new_dataset()


add_records("CVM_iter1_V_A", "PCV")