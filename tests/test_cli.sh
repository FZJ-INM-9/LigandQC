# # single ligand id, no restraints
# ligandqc inspect \
#     tests/input_files/single_ligand.txt \
#     tests/output_files/single_ligand.csv

# # single ligand id, with pdb restraints
# ligandqc inspect \
#     tests/input_files/single_ligand.txt \
#     tests/output_files/single_ligand_single_pdb.csv \
#     --pdb-constraint-file tests/input_files/single_pdbid.txt

# single ligand id, with pdb restraints and chain restraints
ligandqc inspect \
    tests/input_files/single_ligand.txt \
    tests/output_files/single_ligand_single_pdb.yaml \
    --pdb-constraint-file tests/input_files/single_pdbid.txt \
    --only-chain A
