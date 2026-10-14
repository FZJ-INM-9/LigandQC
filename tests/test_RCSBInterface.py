import pytest
from ligand_inspector.rcsb import RCSBInterface
import pandas as pd
import json

def test_get_all_pdbids():
    """
    Test the get_all_pdbids method of the RCSBInterface class.
    """
    rcsb_interface = RCSBInterface()
    ligand_id = "4WI"
    pdb_ids = rcsb_interface.get_all_pdbids([ligand_id])
    
    expected_pdb_ids = ['7RFS', '7RFW', '7SI9', '7TE0', '7TLL', '7U28', '7U29', '7UJU', '7UUP', '7VH8', '7VLO', '7VLP', '7VLQ', '7VTC', '7XB4', '7YRZ', '8C9U', '8DCZ', '8DZ2', '8DZ6', '8DZ9', '8DZA', '8E1Y', '8E25', '8E26', '8E5C', '8EYJ', '8H4Y', '8H51', '8H57', '8H5F', '8H5P', '8H6N', '8H7W', '8H82', '8HOZ', '8HVK', '8HVL', '8HVM', '8HVN', '8HVO', '8HZR', '8IGY', '8INU', '8INW', '9AUJ', '9AUK', '9AUL', '9AUM', '9AVQ', '9BO4', '9BO5', '9BO6', '9BO7']

    assert set(expected_pdb_ids).issubset(set(pdb_ids)), f"Expected {expected_pdb_ids} but got {pdb_ids}"
    assert len(pdb_ids) > 0, "No PDB IDs returned"
    assert isinstance(pdb_ids, list), "Returned value is not a list"

def test_get_ligand_quality_data():
    interface = RCSBInterface()
    ligand_id = "4WI"
    data = interface.get_ligand_quality_data([ligand_id], pdb_ids=["7TLL"])
    assert isinstance(data, list), "Returned value is not a list"
    assert len(data) > 0, "No data returned"
    print(data)
    assert isinstance(data[0], dict), "First element of data is not a dictionary"
    print(data[0].keys())
    keys = data[0].keys()
    for i in ["rcsb_id", "nonpolymer_entities", "exptl"]:
        assert i in keys, f"Key {i} not found in data"
    pdb_id = data[0]["rcsb_id"]
    assert pdb_id == "7TLL", f"Expected '7TLL' but got {pdb_id}"
    nonpolymer_entities = data[0]["nonpolymer_entities"]
    assert isinstance(nonpolymer_entities, list), "nonpolymer_entities is not a list"
    assert len(nonpolymer_entities) == 1, "nonpolymer_entities is empty"
    assert isinstance(nonpolymer_entities[0], dict), "First element of nonpolymer_entities is not a dictionary"
    assert nonpolymer_entities[0]["pdbx_entity_nonpoly"].get("comp_id") == "4WI", f"Expected '4WI' but got {nonpolymer_entities[0]['pdbx_entity_nonpoly'].get('comp_id')}"
    nonpolymer_entity_instances = nonpolymer_entities[0]["nonpolymer_entity_instances"]
    assert isinstance(nonpolymer_entity_instances, list), "nonpolymer_entity_instances is not a list"
    assert len(nonpolymer_entity_instances) == 2, "nonpolymer_entity_instances is empty"


def test_data_processing():
    # file = "tests/input_files/raw_rcsb_data.json"
    file = "intermediate_results.json"
    with open(file, 'r') as f:
        data = json.load(f)
    
    file = "tests/input_files/ligand_ids.txt"
    with open(file, 'r') as f:
        ligand_ids = f.read().splitlines()
        ligand_ids = [ligand_id.strip() for ligand_id in ligand_ids]
        while "" in ligand_ids:
            ligand_ids.remove("")
    
    file = "tests/input_files/pdb_ids.txt"
    with open(file, 'r') as f:
        pdb_ids = f.read().splitlines()
        pdb_ids = [pdb_id.strip() for pdb_id in pdb_ids]
        while "" in pdb_ids:
            pdb_ids.remove("")
    
    interface = RCSBInterface()
    processed_data = interface.process_ligand_quality_data(data, ligand_ids, pdb_ids)
    assert isinstance(processed_data, pd.DataFrame), "Returned value is not a DataFrame"
    assert not processed_data.empty, "Processed data is empty"
