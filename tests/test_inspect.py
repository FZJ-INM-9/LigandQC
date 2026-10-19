import pytest
from ligandqc.inspect import *
from ligandqc.const import *
import pandas as pd
import numpy as np


def test_ligandqc_normal():
    pdb_id = "5REJ"
    ligand_id = "T1V"
    li = LigandInspector()
    df = li.inspect(
        ligand_ids=[ligand_id],
        pdb_ids=[pdb_id],
    )
    df = li.evaluate(df)
    print(df.columns)
    print(df.loc[0])
    assert df is not None
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df.loc[0]["PDB ID"] == pdb_id
    assert df.loc[0]["Ligand ID"] == ligand_id
    assert df.loc[0]["Is Subject of Investigation"]
    assert isinstance(df.at[0, "Ranking for goodness of fit"], float)
    assert isinstance(df.at[0, "Ranking for geometry"], float)
    assert isinstance(df.at[0, METRIC_RSR], float)
    assert isinstance(df.at[0, METRIC_RSCC], float)
    assert isinstance(df.at[0, METRIC_RMSZ_BOND_LENGTH], float)
    assert isinstance(df.at[0, METRIC_RMSZ_BOND_ANGLE], float)
    assert isinstance(df.at[0, "Outliers of bond length"], np.int64)
    assert isinstance(df.at[0, "Outliers of bond angle"], np.int64)
    assert isinstance(df.at[0, METRIC_ATOMIC_CLASHES], np.int64)
    assert isinstance(df.at[0, METRIC_STEREO_CHEMICAL_ERRORS], np.int64)
    assert isinstance(df.at[0, "Model completeness"], np.float64)
    assert isinstance(df.at[0, METRIC_AVG_OCCUPANCY], float)
    assert isinstance(df.at[0, "Number of bonds"], np.int64)
    assert isinstance(df.at[0, "Number of angles"], np.int64)
    assert df.loc[0][f"{METRIC_RSCC}{EVAL_SUFFIX}"] == "bad"
    assert df.loc[0][f"{METRIC_RMSZ_BOND_LENGTH}{EVAL_SUFFIX}"] == "good"
    assert df.loc[0][f"{METRIC_RMSZ_BOND_ANGLE}{EVAL_SUFFIX}"] == "good"
    assert df.loc[0][f"{METRIC_STRUCTURAL_OUTLIERS}{EVAL_SUFFIX}"] == "good"
    assert df.loc[0][f"{METRIC_ATOMIC_CLASHES}{EVAL_SUFFIX}"] == "good"
    assert df.loc[0][f"{METRIC_STEREO_CHEMICAL_ERRORS}{EVAL_SUFFIX}"] == "good"


def test_subject_of_investigation():
    pdb_ids = ["5REJ", "5RTV"]
    ligand_ids = ["DMS"]
    li = LigandInspector()
    df = li.inspect(ligand_ids=ligand_ids, pdb_ids=pdb_ids, only_SOI=True)
    df = li.evaluate(df)
    assert df is not None
    print(df[["Ligand ID", "PDB ID", "Is Subject of Investigation"]])
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert len(df) == 1
    assert df.loc[0]["PDB ID"] == "5RTV"
    assert df.loc[0]["Ligand ID"] == "DMS"
    assert df.loc[0]["Is Subject of Investigation"]


# def test_ligandqc_invalid_pdb_id():
#     pdb_id = "INVALID"
#     ligand_id = "T1V"
#     li = LigandInspector()
#     with pytest.raises(ValueError):
#         li.inspect_single(pdb_id, ligand_id)

# def test_ligandqc_invalid_ligand_id():
#     pdb_id = "5REJ"
#     ligand_id = "INVALID"
#     li = LigandInspector()
#     li.inspect_single(pdb_id, ligand_id)

# def test_ligandqc_invalid_pdb_and_ligand_id():
#     pdb_id = "INVALID"
#     ligand_id = "INVALID"
#     li = LigandInspector()
#     li.inspect_single(pdb_id, ligand_id)


def test_ligandqc_batch_inspection():
    pdb_ids = ["5REJ", "7TLL", "5REK"]
    ligand_ids = ["T1V", "4WI", "T1Y"]
    li = LigandInspector()
    results = li.inspect(ligand_ids=ligand_ids, pdb_ids=pdb_ids)
    assert results is not None
    assert isinstance(results, pd.DataFrame)
    assert not results.empty
    assert len(results) == 4


def test_plot():
    file = "tests/input_files/output_table_without_sdf.csv"
    df = pd.read_csv(file)
    li = LigandInspector()
    li.plot(df, "tests/output_files/plot.png")


# @pytest.mark.parametrize("ligand_id, expected_bond_count", [
#     ("EOF", 30),
#     ("4WI", 37),
# ])
# def test_bond_count(ligand_id, expected_bond_count):
#     """
#     Test the bond count of a ligand.

#     Args:
#         ligand_id (str): The ID of the ligand to test.
#         expected_bond_count (int): The expected number of bonds in the ligand.
#     """
#     inspector = LigandInspector()
#     file = f"tests/input_files/{ligand_id}_ideal.sdf"
#     mol = Chem.SDMolSupplier(file, removeHs=False)[0]
#     bond_count = inspector.get_bond_count(mol)
#     assert bond_count == expected_bond_count, f"Expected {expected_bond_count} bonds, but got {bond_count} for ligand {ligand_id}"

# @pytest.mark.parametrize("ligand_id, expected_angle_count", [
#     ("EOF", 39),
#     ("4WI", 57),
# ])
# def test_angle_count(ligand_id, expected_angle_count):
#     """
#     Test the angle count of a ligand.

#     Args:
#         ligand_id (str): The ID of the ligand to test.
#         expected_angle_count (int): The expected number of angles in the ligand.
#     """
#     inspector = LigandInspector()
#     file = f"tests/input_files/{ligand_id}_ideal.sdf"
#     mol = Chem.SDMolSupplier(file, removeHs=False)[0]
#     angle_count = inspector.get_angle_count(mol)
#     assert angle_count == expected_angle_count, f"Expected {expected_angle_count} angles, but got {angle_count} for ligand {ligand_id}"

# def test_bond_count_visualization():
#     """
#     Test the bond count visualization of a ligand.
#     """
#     inspector = LigandInspector()
#     file = "tests/input_files/EOF_ideal.sdf"
#     mol = Chem.SDMolSupplier(file, removeHs=False)[0]
#     inspector.visualize_bond_count(mol, "tests/output_files/EOF_bond_count.svg")

def produce_artificial_data():
    """
    Produce artificial data for testing purposes. each call of this function yields one data dictionary with the following keys:
    - METRIC_RSR
    - METRIC_RSCC
    - METRIC_RMSZ_BOND_LENGTH
    - METRIC_RMSZ_BOND_ANGLE
    - METRIC_STRUCTURAL_OUTLIERS
    - METRIC_AVG_OCCUPANCY
    - METRIC_ATOMIC_CLASHES
    - METRIC_STEREO_CHEMICAL_ERRORS
    - Number of bonds
    - Number of angles
    - PDB ID
    - Ligand ID
    The values of the metrics are randomly generated within reasonable bounds. The PDB ID and Ligand ID are generated as "1ABC_i" and "LIG_i" where i is a random integer. The number of bonds and angles are also randomly generated within reasonable bounds.
    """

    params = [
		# metric, lb, ub, type
		(METRIC_RSR, 0.0, 1.0, float),
        (METRIC_RSCC, 0.0, 1.0, float),
        (METRIC_RMSZ_BOND_LENGTH, 0.9, 3.0, float),
        (METRIC_RMSZ_BOND_ANGLE, 0.9, 5.0, float),
        (METRIC_STRUCTURAL_OUTLIERS, 0, 20, int),
        (METRIC_AVG_OCCUPANCY, 0.0, 1.0, float),
        (METRIC_ATOMIC_CLASHES, 0, 5, int),
        (METRIC_STEREO_CHEMICAL_ERRORS, 0, 5, int),
    ]
    idx = 0
    while True:
        data = {}
        for metric, lb, ub, typ in params:
            if typ == float:
                data[metric] = lb + (ub - lb) * np.random.rand()
            elif typ == int:
                data[metric] = np.random.randint(lb, ub + 1)
        data["Number of bonds"] = np.random.randint(1, 100, 10)
        data["Number of angles"] = np.random.randint(1, 100, 10)
        data["PDB ID"] = [f"1ABC_{i}" for i in range(10)]
        data["Ligand ID"] = [f"LIG_{i}" for i in range(10)]
        yield data
        idx += 1

def test_produce_artificial_data():
    """
    Test the produce_artificial_data function.
    """
    generator = produce_artificial_data()
    data = [next(generator) for _ in range(3)]
    df = pd.DataFrame(data)
    df.to_csv("tests/output_files/artificial_data.csv", index=False)