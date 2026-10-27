# [Machine Learning for HUH Endonuclease-DNA Specifity Prediction]

HUH-endonucleases are a superfamily of endonuclease enzymes that form stable protein-DNA adducts, culminating in their use as HUH-Tags with many applications in biotechnology. Tompkins et al. identified the main binding motif (named the sDBM) hypothesised
to determine DNA recognition specificity, and facilitates the use of a much simpler machine learning objective in developing a 
predictive model relating to protein-DNA specificity.

This project develops a basic bioinformatics pipeline to obtain a dataset of HUH-endonuclease proteins with sDBM and ori pairs identified, and trains two machine learning models relating to HUH endonuclease-DNA specificity: a binary classifier to predict
DNA binding/cleavage given a ori-sDBM pair, and a ori predictor that predicts the 9-mer ori nic DNA sequence cleaved/bound to a given
sDBM.

## Project Structure

'''
huh-dna-specificity/
├── data/
│   └── embeddings/          # Input tensors (download or generate separately)
│   └── processed/
│       ├── cleaned_datasets/                 
│       ├── final_datasets/           # ready for feature extraction
|       ├── MSA_fasta/                # single fasta files prepped for phylogentic analysis
│       └── negative_datasets/        # negative datasets for each rep and relaxase datasets   
│   └── raw/                # raw protein, reference genome fasta files for viruses, along with initial dataset (raw from pipeline)
├── src/
│   ├── models/                        # neural network architectures
│   ├── data/                          # contains the scripts for the rep-ori pipeline, along those for web-parser/analysis
│   ├── data_setup_classifier.py       # loads data (+ dataloader) for binary classifier model
│   ├── data_setup_predictor.py        # loads data (+ dataloader) for ori predictor model
│   ├── losses.py                      # custom weighted cross-entropy function
│   └── train_classifer.py             # training entry point for classifier models
│   └── train_predictor.py             # training entry point for predictor models
├── notebooks/
│   ├── data_cleaning/                 # scripts used for cleaning the dataset, and obtaining the final (concatenated) csv dataset
│   ├── data_visualisation/            # analysing ori alignment results, PCA, UMAP
│   ├── exploratory_notebooks/         # sketches and initial environment for ML
│   └── feature_extraction.ipynb       # generating input embeddings for models             
├── requirements.txt
├── .gitignore
├── pyproject.toml
└── README.md
'''

## Citations

- Tompkins KJ, Houtti M, Litzau LA, Aird EJ, Everett BA, Nelson AT, Pornschloegl L, Limón-Swanson LK, Evans RL, Evans K, Shi K, Aihara H, Gordon WR. Molecular underpinnings of ssDNA specificity by Rep HUH-endonucleases and implications for HUH-tag multiplexing and engineering. Nucleic Acids Res. 2021 Jan 25;49(2):1046-1064. doi: 10.1093/nar/gkaa1248. PMID: 33410911; PMCID: PMC7826260.