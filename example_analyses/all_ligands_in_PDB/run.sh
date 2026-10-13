
mkdir -p analysis/all_ligands_in_PDB/data
cd analysis/all_ligands_in_PDB

# Download the ligand ids
python build_dataset.py data/ligand_ids.txt data/pdb_ids.txt data/ligands2pdbs.txt

# inspect and plot the ligand quality
ligand_inspector inspect \
    data/ligand_ids.txt \
    data/ligand_quality_data.csv \
    --pdb-constraint-file data/pdb_ids.txt \
    --only-SOI

ligand_inspector plot \
    data/ligand_quality_data.csv \
    data/plot.png
    
cd -
    