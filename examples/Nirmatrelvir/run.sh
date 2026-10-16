#!/bin/bash

# This script is meant to be run from the root of the repository
cd examples/Nirmatrelvir

ligands_file="ligands.txt"
quality_csv="quality.csv"
plot_png="quality.png"

cat > "$ligands_file" << EOF
4WI
EOF

ligandqc inspect \
    --only-SOI \
    --only-chain A \
    -i "$ligands_file" \
    -o "$quality_csv"

ligandqc plot \
    -i "$quality_csv" \
    -o "$plot_png"