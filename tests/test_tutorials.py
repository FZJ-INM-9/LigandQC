import os

import pandas as pd
import pytest
import yaml

from ligandqc.inspect import LigandInspector
from ligandqc.report import ReportGenerator

os.makedirs("tests/output_files/tutorials", exist_ok=True)

@pytest.mark.skip(reason="These tests are for tutorial examples and should be run manually to verify the tutorial outputs as the data may change over time. Therefore these tests are likely failing in CI/CD pipelines.")
class TestTutorials:
    """
    Test suite for tutorial examples.
    """

    def test_tut1(self):
        """
        Test the first tutorial example.
        """
        ligand_ids = ["4WI"]
        pdb_ids = ["7TLL"]

        inspector = LigandInspector()
        df = inspector.inspect(
            ligand_ids=ligand_ids,
            pdb_ids=pdb_ids,
            only_chain="A",
            only_SOI=True,
        )
        df = inspector.evaluate(df)

        with open("tests/output_files/tutorials/output.yaml", "w") as f:
            yaml.dump(df.to_dict('records'), f)

        # This is where the tutorial ends
        
        with open("tests/output_files/tutorials/output.yaml", "r") as f:
            content1 = f.read()

        with open("docs/source/_static/tutorial_files/tut1/output.yaml", "r") as f:
            content2 = f.read()

        assert content1 == content2


    def test_tut2(self, tmp_path):
        """
        Test the second tutorial example.
        """

        ligand_ids = ["4WI"]

        inspector = LigandInspector()
        df = inspector.inspect(ligand_ids=ligand_ids, only_SOI=True)
        df = inspector.evaluate(df)

        df.to_csv(tmp_path / "quality.csv", index=False)

        inspector.plot(df, tmp_path / "quality.png")


        df = pd.read_csv(tmp_path / "quality.csv")
        report = ReportGenerator(df)
        report.generate_report(tmp_path)

        # This is where the tutorial ends

        assert os.path.exists(tmp_path / "quality.csv")
        assert os.path.exists(tmp_path / "quality.png")

        with open(tmp_path / "quality.csv", "r") as f:
            content1 = f.read()

        with open("docs/source/_static/tutorial_files/tut2/quality.csv", "r") as f:
            content2 = f.read()

        assert content1 == content2

        # assert terminal report content
        
    def test_tut3(self, tmp_path):
                    

        pdb_ids = ["5REK", "5REJ", "7TLL"]

        inspector = LigandInspector()
        df = inspector.inspect(pdb_ids=pdb_ids, only_SOI=True)
        df = inspector.evaluate(df)

        df.to_csv(tmp_path / "quality.csv", index=False)

        inspector.plot(df, tmp_path / "quality.png")


        # This is where the tutorial ends
        assert os.path.exists(tmp_path / "quality.csv")
        assert os.path.exists(tmp_path / "quality.png")

        with open(tmp_path / "quality.csv", "r") as f:
            content1 = f.read()

        with open("docs/source/_static/tutorial_files/tut3/quality.csv", "r") as f:
            content2 = f.read()

        assert content1 == content2

