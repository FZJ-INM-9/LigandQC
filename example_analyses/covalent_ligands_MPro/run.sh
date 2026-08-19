# is meant to be run from the root of the repository

cd analysis/covalent_ligands_MPro

ligand_inspector inspect \
    data/ligand_ids.txt \
    data/ligand_quality_data.csv \
    --pdb-constraint-file data/pdb_ids.txt \
    --only-SOI \
    --only-chain "A"

ligand_inspector plot \
    data/ligand_quality_data.csv \
    data/plot.png

cd -