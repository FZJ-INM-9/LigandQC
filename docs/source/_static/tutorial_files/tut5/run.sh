echo 4WI > ligands.txt

ligandqc download \
    --ligand-id-file ligands.txt \
    --output-directory .