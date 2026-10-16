cat > ligands.txt << EOF
4WI
EOF

ligandqc inspect \
    --input-file ligands.txt \
    --output-file ligands_quality.csv \
    --only-SOI

ligandqc plot \
    --input-file ligands_quality.csv \
    --output-file ligands_quality.png