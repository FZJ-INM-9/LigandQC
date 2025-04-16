# single ligand id, no restraints
ligand_inspector inspect \
    tests/input_files/single_ligand.txt \
    tests/output_files/single_ligand.csv

# single ligand id, with pdb restraints
ligand_inspector inspect \
    tests/input_files/single_ligand.txt \
    tests/output_files/single_ligand_single_pdb.csv \
    --pdb-constraint-file tests/input_files/single_pdbid.txt

# single ligand id, with pdb restraints and chain restraints
ligand_inspector inspect \
    tests/input_files/single_ligand.txt \
    tests/output_files/single_ligand_single_pdb.csv \
    --pdb-constraint-file tests/input_files/single_pdbid.txt \
    --only-chain A
