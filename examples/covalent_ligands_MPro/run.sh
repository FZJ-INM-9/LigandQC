#!/bin/bash

# This script is meant to be run from the root of the repository

cd examples/covalent_ligands_MPro

pdb_ids_file="pdb_ids.txt"
quality_csv="quality.csv"
plot_png="quality.png"

ligandqc inspect \
    --output-file "$quality_csv" \
    --pdb-constraint-file "$pdb_ids_file" \
    --only-SOI \
    --only-chain "A" 

ligandqc plot \
    --input-file "$quality_csv" \
    --output-file "$plot_png"

cd -