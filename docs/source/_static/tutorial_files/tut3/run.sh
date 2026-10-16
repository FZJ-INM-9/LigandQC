cat > pdb_ids.txt << EOF
5REK
7TLL
5REJ
EOF

ligandqc inspect \
    --pdb-constraint-file pdb_ids.txt \
    --output-file quality.csv \
    --only-SOI

ligandqc plot \
    --input-file quality.csv \
    --output-file quality.png