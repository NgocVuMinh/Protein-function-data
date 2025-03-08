# Extract A.thaliana, O.sativa, and plant proteins from UniProt SwissProt
python extract_sprot_data.py

# Parse sequences into fasta files to run DIAMOND
python pkl2fasta.py

# Run DIAMOND (later: DeepGOSE split the data based on sequence similarity)
diamond makedb --in data/sativa_exp.fa --db data/sativa_exp.dmnd
diamond blastp --very-sensitive -d data/sativa_exp.dmnd -q data/sativa_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/sativa_exp.sim

diamond makedb --in data/thaliana_exp.fa --db data/thaliana_exp.dmnd
diamond blastp --very-sensitive -d data/thaliana_exp.dmnd -q data/thaliana_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/thaliana_exp.sim

diamond makedb --in data/plant_exp.fa --db data/plant_exp.dmnd
diamond blastp --very-sensitive -d data/plant_exp.dmnd -q data/plant_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/plant_exp.sim

diamond makedb --in data/swissprot_exp.fa --db data/swissprot_exp.dmnd
diamond blastp --very-sensitive -d data/swissprot_exp.dmnd -q data/swissprot_exp.fa --outfmt 6 qseqid sseqid bitscore pident > data/swissprot_exp.sim

# Getting PPI network information from STRING
python ppi.py 

# Getting ESM embeddings
python esm2.py 

# Split the data into sub-ontologies (BP, MF, and CC) with train, validation, and test sets
python split_ontology.py -pre sativa_exp -df data/sativa_exp.pkl -sf data/sativa_exp.sim
python split_ontology.py -pre thaliana_exp -df data/thaliana_exp.pkl -sf data/thaliana_exp.sim
python split_ontology.py -pre plant_exp -df data/plant_exp.pkl -sf data/plant_exp.sim
python split_ontology.py -pre swissprot_exp -df data/swissprot_exp.pkl -sf data/swissprot_exp.sim
