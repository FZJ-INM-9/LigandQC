import pytest

from ligand_inspector.LigandDownloader import LigandDownloader


def compare_sdf_file_contents(content1, content2):
    """
    Compare two SDF files, ignoring the second line.

    Args:
        content1 (str): Content of the first SDF file.
        content2 (str): Content of the second SDF file.
    Returns:
        bool: True if the files are equal (ignoring the second line), False otherwise.
    """
    lines1 = content1.splitlines()
    lines2 = content2.splitlines()

    # Remove the second line from both files
    if len(lines1) > 1:
        lines1.pop(1)
    if len(lines2) > 1:
        lines2.pop(1)

    return lines1 == lines2


@pytest.mark.parametrize(
    "ligand_ids, expected_url",
    [
        (["4WI"], "https://download.rcsb.org/batch/ccd/4WI_ideal.sdf"),
        (["EOF"], "https://download.rcsb.org/batch/ccd/EOF_ideal.sdf"),
        (
            ["EOF", "4WI"],
            "https://download.rcsb.org/batch/ccd/EOF_ideal.sdf:4WI_ideal.sdf",
        ),
    ],
)
def test_build_urls(ligand_ids, expected_url):
    """
    Test the build_url method of the LigandDownloader class.
    """
    downloader = LigandDownloader()
    urls = downloader.build_urls(ligand_ids)
    assert urls[0] == expected_url


@pytest.mark.parametrize(
    "ligand_id, sdf_file",
    [
        ("EOF", "tests/input_files/EOF_ideal.sdf"),
        ("4WI", "tests/input_files/4WI_ideal.sdf"),
    ],
)
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

    assert compare_sdf_file_contents(ligand_data, expected_data)


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
    for ligand_id, sdf_file in zip(
        ligand_ids,
        ["tests/input_files/EOF_ideal.sdf", "tests/input_files/4WI_ideal.sdf"],
    ):
        with open(sdf_file, "r") as f:
            expected_data = f.read()

        assert compare_sdf_file_contents(ligands[ligand_id], expected_data)
