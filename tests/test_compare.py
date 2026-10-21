import pandas as pd
import pytest
from ligandqc.compare import DatasetComparison
import os


class TestDatasetComparison:

    @pytest.fixture
    def data1(self):
        file = "tests/input_files/example_result_table.csv"
        return pd.read_csv(file)

    @pytest.fixture
    def data2(self):
        file = "tests/input_files/example_result_table2.csv"
        return pd.read_csv(file)

    def test_compare_datasets(self, data1, data2):
        comparison = DatasetComparison(data1, data2)
        comparison.compare_datasets(output_dir="tests/output_files/comparison_results")

        assert os.path.exists("tests/output_files/comparison_results/venn_diagram_ligand_ids.png")
        assert os.path.exists("tests/output_files/comparison_results/venn_diagram_pdb_id_ligand_id_pairs.png")
        assert os.path.exists("tests/output_files/comparison_results/comparative_pareto_plot_good.png")
        assert os.path.exists("tests/output_files/comparison_results/comparative_pareto_plot_ok.png")
        assert os.path.exists("tests/output_files/comparison_results/comparative_pareto_plot_bad.png")
        assert os.path.exists("tests/output_files/comparison_results/comparison_single_metric.csv")
        assert os.path.exists("tests/output_files/comparison_results/comparison_across_metrics.csv")
        assert os.path.exists("tests/output_files/comparison_results/comparison.html")
        

