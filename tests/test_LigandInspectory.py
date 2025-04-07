import pytest
from ligand_inspector.LigandInspector import *
import pandas as pd

def test_ligand_inspector():
    pdb_id = "5REJ"
    ligand_id = "T1V"
    li = LigandInspector()
    df = li.inspect_single(pdb_id, ligand_id)
    assert df is not None
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df.loc[0]["PDB ID"] == pdb_id
    assert df.loc[0]["Ligand ID"] == ligand_id
    assert df.loc[0]["Identifier"] == "5REJ_T1V_A_404"
    assert df.loc[0]["Ranking for goodness of fit"] == 0.15
    assert df.loc[0]["Ranking for geometry"] == 0.76
    assert df.loc[0]["Real space R factor"] == 0.199
    assert df.loc[0]["Real space correlation coefficient"] == 0.755
    assert df.loc[0]["RMSZ-bond-length"] == 0.4
    assert df.loc[0]["RMSZ-bond-angle"] == 0.53 	
    assert df.loc[0]["Outliers of bond length"] == 0
    assert df.loc[0]["Outliers of bond angle"] == 0
    assert df.loc[0]["Atomic clashes"] == 0
    assert df.loc[0]["Stereochemical errors"] == 0
    assert df.loc[0]["Model completeness"] == 1
    assert df.loc[0]["Average occupancy"] == 0.64
