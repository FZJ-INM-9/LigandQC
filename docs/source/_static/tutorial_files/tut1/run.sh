cat > ligand_ids.txt << EOF
4WI
EOF

cat > pdb_ids.txt << EOF
7TLL
EOF

ligandqc inspect \
    --ligand-id-file ligand_ids.txt \
    --output-file output.yaml \
    --pdb-constraint-file pdb_ids.txt \
    --only-chain "A" \
    --only-SOI