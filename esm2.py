import os
import sys
import gzip
import click as ck
import numpy as np
import pandas as pd
from collections import Counter, deque
from deepgo.utils import (
    Ontology, FUNC_DICT, NAMESPACES, MOLECULAR_FUNCTION, BIOLOGICAL_PROCESS,
    CELLULAR_COMPONENT, HAS_FUNCTION,
    is_cafa_target)
from utils import is_exp_code, is_exp_plus_code, load_data
from deepgo.extract_esm import extract_esm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

device = "cuda"

def main(data_file):
    df = pd.read_pickle(data_file)
    print("Data: ", df.shape)

    fasta_file = "dummy.fa"
    with open(fasta_file, "w") as f:
        for row in df.itertuples():
            record = SeqRecord(
                Seq(row.sequences),
                id=row.proteins,
                description=""
            )
            SeqIO.write(record, f, "fasta")
    prots, esm2_data = extract_esm(fasta_file, device=device)
    esm2_data = list(esm2_data)
    esm2_final = []
    for esm in esm2_data:
        esm2_final.append(esm.numpy())
    df["esm2"] = esm2_final
    df.to_pickle(data_file)
    os.remove(fasta_file)
    print(f"Completed adding esm to {data_file}")
    print(df.shape)

if __name__ == '__main__':
    main("data/sativa_exp.pkl") # overwriting the input file
    main("data/thaliana_exp.pkl") 
    main("data/plant_exp.pkl")
    main("data/swissprot_exp.pkl")