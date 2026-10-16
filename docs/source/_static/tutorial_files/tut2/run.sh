cat > ligands.txt << EOF
4WI
EOF

ligandqc inspect \
    --ligand-id-file ligands.txt \
    --output-file quality.csv \
    --only-SOI

ligandqc plot \
    --input-file quality.csv \
    --output-file quality.png