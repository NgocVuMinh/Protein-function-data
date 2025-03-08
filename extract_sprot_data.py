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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Load ontologies
go = Ontology("go.obo", with_rels=True)

# Load SwissProt
swissprot_file = "data/uniprot_sprot.dat.gz"
proteins, accessions, sequences, annotations, string_ids, orgs, genes, interpros = load_data(swissprot_file)

df = pd.DataFrame({
    "proteins": proteins,
    "accessions": accessions,
    "genes": genes,
    "sequences": sequences,
    "annotations": annotations,
    "string_ids": string_ids,
    "orgs": orgs,
    "interpros": interpros
})

print(f"Number of SwissProt proteins: {len(proteins)}")

# Getting proteins with experimental annotations
index_exp = []
annotations_exp = []
for i, row in enumerate(df.itertuples()):
    annots_exp = []
    for annot in row.annotations:
        go_id, code = annot.split("|")
        if is_exp_code(code):
            annots_exp.append(go_id)
    if len(annots_exp) == 0:
        continue
    index_exp.append(i)
    annotations_exp.append(annots_exp)
df_exp = df.iloc[index_exp]
df_exp = df_exp.reset_index()
print(f"Number of experimentally annotated proteins: {df_exp.shape[0]}")

# Removing obsolete GO annotations
print("Removing obsolete terms...")
valid_go = []
for t in list(go.ont.keys()):
    if go.get_term(t)["is_obsolete"] == False:
        valid_go.append(t)
# print(f"Total number of GO: {len(list(go.ont.keys()))}")       
# print(f"Number of not-obsolete GO: {len(valid_go)}")

exp_cleaned = []
obsolete = []
for annots in annotations_exp:
    annots_cleaned = []
    for annot in annots:
        if annot in valid_go:
            annots_cleaned.append(annot)
        else:
            obsolete.append(annot)
    exp_cleaned.append(annots_cleaned)
print(f"Number of obsolete proteins found in the dataset: {len(list(set(obsolete)))}")

## Propagating annotations
df_exp["exp_annotations"] = exp_cleaned
prop_annotations = []
for i, row in df_exp.iterrows():
    annot_set = set()
    annots = row["exp_annotations"]
    for go_id in annots:
        annot_set |= go.get_ancestors(go_id)
    annots = list(annot_set)
    prop_annotations.append(annots)
df_exp["prop_annotations"] = prop_annotations

print(f"Number of experimentally annotated proteins: {df_exp.shape[0]}")



# Extracting plant proteins:

print("Extracting plant proteins:")

## A.thaliana (Taxonomy ID: 3702)
thaliana_exp = df_exp[df_exp["orgs"] == "3702"]

## O.sativa (Taxonomy ID: 39947 (Japonica group) and 39946 (Indica group))
sativa_exp = pd.concat([df_exp[df_exp["orgs"] == "39947"], df_exp[df_exp["orgs"] == "39946"]])

## Combined (=A.thaliana+O.sativa)
plant_exp = pd.concat([thaliana_exp, sativa_exp])

print(f"A.thaliana (exp): {thaliana_exp.shape[0]}")
print(f"O.sativa (exp): {sativa_exp.shape[0]}")
print(f"Combined (exp): {plant_exp.shape[0]}")

thaliana_exp.to_pickle("data/thaliana_exp.pkl")
sativa_exp.to_pickle("data/sativa_exp.pkl")
plant_exp.to_pickle("data/plant_exp.pkl")
df_exp.to_pickle("data/swissprot_exp.pkl")

