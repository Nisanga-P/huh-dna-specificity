set -e

#no gaps in defining bash variables
INPUT_FILE="taxids.txt"

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

while IFS= read -r taxid || [ -n "$taxid" ]; do
    # Skip empty lines or comments 
    [[ -z "$taxid" || "$taxid" =~ ^# ]] && continue

    # Strip any hidden carriage returns 
    taxid=$(echo "$taxid" | tr -d '\r')

    #add extra line to remove spaces

    echo "----------------------------------------"
    echo "Downloading genome package for TaxID: $taxid from the NCBI database"
    
    #datasets download virus genome taxon "$taxid" --include genome --filename genome_packages/zipped/$taxid
    
    #specifically get ref seq
    datasets download genome taxon "$taxid" --assembly-source RefSeq --filename genome_packages/zipped/$taxid

    echo "Extracting files..."
    unzip -q genome_packages/zipped/$taxid -d genome_packages/extracted/$taxid

    echo "Moving FASTA file to 'Virus_Genomes' directory"
   # cp genome_packages/extracted/$taxid/ncbi_dataset/data/genomic.fna virus_genomes/$taxid.fna
    
    # multiple ref sequences can exist -> pulls all fasta files in the the 'data' subdirectory
    find genome_packages/extracted/$taxid/ncbi_dataset/data -type f -name "*.fna" -exec cp {} virus_genomes/ \;
    

done < "$INPUT_FILE"

rm -rf genome_packages

echo "Finished downloading all virus genomes."