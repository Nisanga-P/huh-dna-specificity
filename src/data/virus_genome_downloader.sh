set -e

trap 'rm -rf genome_packages' EXIT

# ensure NCBI datasets in 'tools' is accessible by the bash script (so a script using it can be run from anywhere)
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
export PATH="$REPO_ROOT/.tools:$PATH" 

#no gaps in defining bash variables
INPUT_FILE="$1"

# The name of the 'psi_blast' run (this must be passed in as an argument!)
SUBFOLDER_NAME="$2"

ALL_VIRUS_RECORD_FILE="data/raw/generated/all_viruses_in_dataset.txt"

#ensure that INPUT text file containing taxids exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File $INPUT_FILE not found!"
    exit 1
fi 

#ensure that NCBI 'datasets' tool is installed
if ! command -v datasets &> /dev/null; then
    echo "Error: NCBI 'datasets' tool is not installed."
    exit 1
fi

echo "Found virus info file (containing taxids)"

mkdir genome_packages
mkdir genome_packages/zipped
mkdir genome_packages/extracted

mkdir data/raw/downloads/virus_genomes/"$SUBFOLDER_NAME"

echo "Starting processing..."

while IFS= read -r line; do

    # Skip empty lines or comments 
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    # Strip any hidden carriage returns 
    line=$(echo "$line" | tr -d '\r')

    IFS=',' read -r organism_name taxid <<< "$line"

    echo "----------------------------------------"

    while IFS= read -r second_line; do
        second_line=${second_line//$'\r'/}

        [[ -z "$second_line" || "$second_line" =~ ^# ]] && continue

        if [[ "$second_line" == "$taxid" ]]; then
            echo "'$organism_name' genome already exists in virus_genomes. Skipping..."
            echo "----------------------------------------"
            continue 2
        fi

    done < "$ALL_VIRUS_RECORD_FILE"



    
    # # do check to maker sure that it has not already been downloaded
    # if [ -d "data/raw/downloads/virus_genomes/${taxid}_${organism_name// /_}" ]; then
    #     echo "'$organism_name' genome already exists in virus_genomes. Skipping..."
    #     echo "----------------------------------------"
    #     continue
    # fi

    echo "New virus!"
    echo "Downloading genome package for Organism: $organism_name with taxid $taxid from the NCBI database"
    
    #datasets download virus genome taxon "$taxid" --include genome --filename genome_packages/zipped/$taxid
    datasets download genome taxon "$taxid" --assembly-source RefSeq --filename genome_packages/zipped/$taxid

    echo "Extracting files..."
    unzip -q genome_packages/zipped/$taxid -d genome_packages/extracted/$taxid

    echo "Moving FASTA file to 'Virus_Genomes' directory"

    mkdir data/raw/downloads/virus_genomes/"$SUBFOLDER_NAME"/${taxid}_${organism_name// /_}

   # cp genome_packages/extracted/$taxid/ncbi_dataset/data/genomic.fna virus_genomes/$taxid.fna
    find genome_packages/extracted/$taxid/ncbi_dataset/data -type f -name "*.fna" -exec cp {} data/raw/downloads/virus_genomes/"$SUBFOLDER_NAME"/${taxid}_${organism_name// /_} \;

    echo "----------------------------------------"
   
done < "$INPUT_FILE"

echo "Finished downloading all virus genomes."