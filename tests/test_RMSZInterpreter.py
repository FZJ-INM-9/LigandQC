import pytest

from ligand_inspector.RMSZInterpreter import RMSZInterpreter
from rdkit import Chem

def test_rmsz_threshold_finding():
    """
    Test the find_rmsz_threshold method of the RMSZInterpreter class.
    """
    interpreter = RMSZInterpreter()
    
    # Test with a sample size of 2
    n = 6
    rmsz_threshold = interpreter.find_rmsz_threshold(n)
    
    # Check if the threshold is within expected bounds
    assert rmsz_threshold == 1.48799

def test_rational_plot():
    """
    Test the rational_plot method of the RMSZInterpreter class.
    """
    interpreter = RMSZInterpreter()
    interpreter.plot_rationale()

def test_plot_explanation():
    interpreter = RMSZInterpreter()
    interpreter.plot_explanation()


@pytest.mark.parametrize("ligand_id, expected_bond_count", [
    ("EOF", 30),
    ("4WI", 37),
])
def test_bond_count(ligand_id, expected_bond_count):
    """
    Test the bond count of a ligand.
    
    Args:
        ligand_id (str): The ID of the ligand to test.
        expected_bond_count (int): The expected number of bonds in the ligand.
    """
    interpreter = RMSZInterpreter()
    file = f"tests/input_files/{ligand_id}_ideal.sdf"
    mol = Chem.SDMolSupplier(file, removeHs=False)[0]
    bond_count = interpreter.get_bond_count(mol)
    assert bond_count == expected_bond_count, f"Expected {expected_bond_count} bonds, but got {bond_count} for ligand {ligand_id}"

@pytest.mark.parametrize("ligand_id, expected_angle_count", [
    ("EOF", 39),
    ("4WI", 57),
])
def test_angle_count(ligand_id, expected_angle_count):
    """
    Test the angle count of a ligand.
    
    Args:
        ligand_id (str): The ID of the ligand to test.
        expected_angle_count (int): The expected number of angles in the ligand.
    """
    interpreter = RMSZInterpreter()
    file = f"tests/input_files/{ligand_id}_ideal.sdf"
    mol = Chem.SDMolSupplier(file, removeHs=False)[0]
    angle_count = interpreter.get_angle_count(mol)
    assert angle_count == expected_angle_count, f"Expected {expected_angle_count} angles, but got {angle_count} for ligand {ligand_id}"

def test_bond_count_visualization():
    """
    Test the bond count visualization of a ligand.
    """
    interpreter = RMSZInterpreter()
    file = "tests/input_files/EOF_ideal.sdf"
    mol = Chem.SDMolSupplier(file, removeHs=False)[0]
    interpreter.visualize_bond_count(mol, "tests/output_files/EOF_bond_count.svg")


@pytest.mark.parametrize("rmsz_value, mode, expected_interpretation", [
    (0.2, "bond",  "good"),
    (0.1, "angle", "good"),
    (1.2, "bond",  "good"),
    (1.1, "angle", "good"),
    (1.3, "bond", "ok"),
    (1.7, "angle", "ok"),
    (2.0, "bond",  "ok"),
    (2.0, "angle", "ok"),
    (2.5, "bond",  "bad"),
    (2.5, "angle", "bad"),
    (3.0, "bond",  "bad"),
    (3.0, "angle", "bad"),
    (4.0, "bond",  "bad"),
    (4.0, "angle", "bad"),
])
def test_molecule_interpretation(rmsz_value, mode, expected_interpretation):
    """
    Test the interpretation of a molecule.
    """
    interpreter = RMSZInterpreter()
    file = "tests/input_files/EOF_ideal.sdf"
    mol = Chem.SDMolSupplier(file, removeHs=False)[0]
    
    # Use the find_rmsz_threshold method to get the threshold
    if mode == "bond":
        rmsz_threshold = interpreter.find_rmsz_threshold_bond(mol)
    else:
        rmsz_threshold = interpreter.find_rmsz_threshold_angle(mol)


    interpretation = interpreter.interpret_molecule(mol, rmsz_value, mode=mode)
    assert interpretation == expected_interpretation, f"Expected '{expected_interpretation}', but got '{interpretation}' for RMSZ value {rmsz_value}, with threshold {rmsz_threshold} and mode {mode}"