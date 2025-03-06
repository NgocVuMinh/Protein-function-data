import os
import sys
import click as ck
import numpy as np
import pandas as pd

def main(data_file, out_file):
    # Load interpro data
    df = pd.read_pickle(data_file)
    print(len(df)) 
    with open(out_file, "w") as f:
        for row in df.itertuples():
            prot_id = row.proteins
            f.write(">" + prot_id + "\n")
            f.write(row.sequences + "\n")

if __name__ == "__main__":
    main("data/thaliana_exp.pkl", "data/thaliana_exp.fa")
    main("data/sativa_exp.pkl", "data/sativa_exp.fa")
    main("data/plant_exp.pkl", "data/plant_exp.fa")
