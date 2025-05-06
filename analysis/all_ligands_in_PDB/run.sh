
mkdir -p analysis/data
python analysis/build_dataset.py analysis/data/ligand_ids.txt analysis/data/pdb_ids.txt analysis/data/ligands2pdbs.txt
ligand_inspector inspect \
    analysis/data/ligand_ids.txt \
    analysis/data/ligand_quality_data.csv \
    --pdb-constraint-file analysis/data/pdb_ids.txt \
    --only-SOI