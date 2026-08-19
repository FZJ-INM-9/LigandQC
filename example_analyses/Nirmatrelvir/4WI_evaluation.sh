cat > ligands.txt << EOF
4WI
EOF

ligand_inspector inspect \
    --only-SOI \
    --only-chain A \
    ligands.txt \
    nirmatrelvir_quality.csv

ligand_inspector plot \
    nirmatrelvir_quality.csv \
    nirmatrelvir_quality.png