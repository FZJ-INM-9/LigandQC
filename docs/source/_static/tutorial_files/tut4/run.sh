# ligandqc get-all-ids \
#     -ol ligands.txt \
#     -op pdb-ids.txt \
#     -om mapping.txt 

ligandqc inspect \
    --ligand-id-file ligands.txt \
    --pdb-constraint-file pdb-ids.txt \
    --output-file quality.csv \
    --only-SOI

ligandqc plot \
    --input-file quality.csv \
    --output-file quality.png