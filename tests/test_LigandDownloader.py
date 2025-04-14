import pytest

from ligand_inspector.LigandDownloader import LigandDownloader



@pytest.mark.parametrize("ligand_id, expected_url", [
    ("4WI", "http://ligand-expo.rcsb.org/reports/4/4WI/4WI_ideal.sdf"),
    ("EOF", "http://ligand-expo.rcsb.org/reports/E/EOF/EOF_ideal.sdf"),
])
def test_build_url(ligand_id, expected_url):
    """
    Test the build_url method of the LigandDownloader class.
    """
    downloader = LigandDownloader()
    url = downloader.build_url(ligand_id)
    assert url == expected_url


@pytest.mark.parametrize("ligand_id, sdf_file", [
    ("EOF", "tests/input_files/EOF_ideal.sdf"),
    ("4WI", "tests/input_files/4WI_ideal.sdf"),
    
])
def test_download_single(ligand_id, sdf_file):
    """
    Test the download_ligand method of the LigandDownloader class.
    """
    downloader = LigandDownloader()
    ligand_data = downloader.download_single_ligand(ligand_id)
    
    # Check if the downloaded data is not empty
    assert ligand_data is not None
    assert len(ligand_data) > 0
    
    # Check if the downloaded data matches the expected SDF file
    with open(sdf_file, "r") as f:
        expected_data = f.read()
    
    assert ligand_data == expected_data

def test_download_multiple():
    """
    Test the download_multiple_ligands method of the LigandDownloader class.
    """
    downloader = LigandDownloader()
    ligand_ids = ["EOF", "4WI"]
    
    # Download multiple ligands
    ligands = downloader.download_multiple_ligands(ligand_ids)
    
    # Check if the downloaded data is not empty
    assert len(ligands) == len(ligand_ids)
    assert isinstance(ligands, list)

    ligands = dict(zip(ligand_ids, ligands))
    
    # Check if the downloaded data matches the expected SDF files
    for ligand_id, sdf_file in zip(ligand_ids, ["tests/input_files/EOF_ideal.sdf", "tests/input_files/4WI_ideal.sdf"]):
        with open(sdf_file, "r") as f:
            expected_data = f.read()
        
        assert ligands[ligand_id] == expected_data