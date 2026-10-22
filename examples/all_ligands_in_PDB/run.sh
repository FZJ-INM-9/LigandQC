#!/bin/bash

# This script will run the analysis for all ligands in the PDB. 
# it is assumed to be run from the root of the repository.

cd examples/all_ligands_in_PDB
ligand_id_file="data/ligand_ids.txt"
pdb_id_file="data/pdb_ids.txt"
quality_csv="data/quality.csv"
quality_png="data/quality.png"

# Download the ligand ids
ligandqc get-all-ids \
    -ol $ligand_id_file \
    -op $pdb_id_file \
    -om data/mapping.txt

# inspect and plot the ligand quality
ligandqc inspect \
    --ligand-id-file $ligand_id_file \
    --output-file $quality_csv \
    --pdb-constraint-file $pdb_id_file \
    --only-SOI

ligandqc plot \
    -i $quality_csv \
    -o $quality_png

ligandqc report \
    -i $quality_csv \
    -o terminal
    
cd -
    