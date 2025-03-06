# Data preparation for DeepGOSE training

### Prepared data:
You can download all the datasets I prepared and used for model training from here: [Google Drive](https://drive.google.com/drive/folders/1e4xQ8)
Alternatively, you can download the data from the original sources as listed below, put them in ```data/```, and follow the following scripts.

### Data sources

The plant datasets were generated from the UniProtKB/Swiss-Prot database version 2024_03: 

- Previous release of UniProt SwissProt (version 2024_03, in my latest report I mistakenly referred to this as 2024-07-22, which is the last-modified date, not the version name, you can check that via this [link](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/)) is stored in:
https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2024_03/knowledgebase/uniprot_sprot-only2024_03.tar.gz

    Unzip the file, you will see a ```uniprot_sprot.dat.gz``` file, which is used to generate the dataset. 

- The GO database version 2024-04-24: [go.obo](https://release.geneontology.org/2024-04-24/ontology/go.obo)

