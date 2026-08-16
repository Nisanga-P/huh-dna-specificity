set -e

#no gaps in defining bash variables
INPUT_FILE="virus_species_info.txt"

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

echo "Found virus taxids file"

mkdir genome_packages
mkdir genome_packages/zipped
mkdir genome_packages/extracted

echo "Starting processing of virus taxids..."

while IFS= read -r line; do

    # Skip empty lines or comments 
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    # Strip any hidden carriage returns 
    line=$(echo "$line" | tr -d '\r')

    IFS=',' read -r taxid organism_name <<< "$line"

    echo "----------------------------------------"
    echo "Downloading genome package for TaxID: $taxid from the NCBI database"
    
    #datasets download virus genome taxon "$taxid" --include genome --filename genome_packages/zipped/$taxid
    datasets download genome taxon "$taxid" --assembly-source RefSeq --filename genome_packages/zipped/$taxid

    echo "Extracting files..."
    unzip -q genome_packages/zipped/$taxid -d genome_packages/extracted/$taxid

    echo "Moving FASTA file to 'Virus_Genomes' directory"

    mkdir virus_genomes/${organism_name// /_}

   # cp genome_packages/extracted/$taxid/ncbi_dataset/data/genomic.fna virus_genomes/$taxid.fna
    find genome_packages/extracted/$taxid/ncbi_dataset/data -type f -name "*.fna" -exec cp {} virus_genomes/${organism_name// /_} \;
    

done < "$INPUT_FILE"

rm -rf genome_packages

echo "Finished downloading all virus genomes."