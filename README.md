# Data preparation for DeepGOSE training


### Dependencies

To prepare the datasets from scratch, a GPU will be need to generate ESM embeddings. Please install [DIAMOND](https://github.com/bbuchfink/diamond) and run ```pip install -r requirements.txt``` for other requirements.

However, if you're only interested in GO specificity, these are needed:
```
pandas
numpy
networkx
nxontology
pronto
matplotlib
scipy
seaborn
statsmodels
```


### Prepared data
You can download the datasets (experimental only) that I prepared and used for model training from this [link](https://drive.google.com/file/d/1HYVDtIwJfZSWouZgZeTpAzqGp1XsuJU4/view?usp=sharing) (about 5.3Gb). Please decompress and place this "data" folder in the current repository.

Alternatively, you can download the raw protein data from UniProt/SwissProt (see the Data sources section below), place it in a new folder ```data/```, and run the following commands (or simply run ```./run.sh```):

Extract A.thaliana, O.sativa, and plant proteins from UniProt SwissProt:

```
python extract_sprot_data.py
```

Parse sequences into fasta files to run DIAMOND:
```
python pkl2fasta.py
```

Run DIAMOND (the data is later split into train-valid-test sets based on sequence similarity)
```
diamond makedb --in data/sativa_exp.fa --db data/sativa_exp.dmnd
diamond blastp --very-sensitive -d data/sativa_exp.dmnd -q data/sativa_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/sativa_exp.sim

diamond makedb --in data/thaliana_exp.fa --db data/thaliana_exp.dmnd
diamond blastp --very-sensitive -d data/thaliana_exp.dmnd -q data/thaliana_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/thaliana_exp.sim

diamond makedb --in data/plant_exp.fa --db data/plant_exp.dmnd
diamond blastp --very-sensitive -d data/plant_exp.dmnd -q data/plant_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/plant_exp.sim

diamond makedb --in data/swissprot_exp.fa --db data/swissprot_exp.dmnd
diamond blastp --very-sensitive -d data/swissprot_exp.dmnd -q data/swissprot_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/swissprot_exp.sim
```

Getting PPI network information from STRING
```
python ppi.py 
```

Getting ESM embeddings, this requires a GPU
```
python esm2.py 
```

Split the data into sub-ontologies (BP, MF, and CC) with train, validation, and test sets
```
python split_ontology.py -pre sativa_exp -df data/sativa_exp.pkl -sf data/sativa_exp.sim
python split_ontology.py -pre thaliana_exp -df data/thaliana_exp.pkl -sf data/thaliana_exp.sim
python split_ontology.py -pre plant_exp -df data/plant_exp.pkl -sf data/plant_exp.sim
python split_ontology.py -pre swissprot_exp -df data/swissprot_exp.pkl -sf data/swissprot_exp.sim
```


### Data sources

The datasets were generated from: 

- UniProtKB/SwissProt ([version 2024_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/), in my latest report I referred to this as 2024-07-22, which is the last-modified date, not the version name). The actual file used to generate the dataset is named ```uniprot_sprot.dat.gz```, which is stored in a compressed folder [uniprot_sprot-only2024_03.tar.gz](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2024_03/knowledgebase/uniprot_sprot-only2024_03.tar.gz)

    Please decompress and place the necessary file in the ```data``` folder

- The GO database version 2024-04-24 (this has already been provided in this repository): [go.obo](https://release.geneontology.org/2024-04-24/ontology/go.obo)


### GO specificity
The GO specificity was assessed according to this notebook: [GO_specificity.ipynb](https://github.com/NgocVuMinh/Protein-function-data/blob/main/GO_specificity.ipynb)

