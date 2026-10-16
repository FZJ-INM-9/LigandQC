import pytest
import pandas as pd
from ligandqc.report import ReportGenerator
import os

class TestReport:

    @pytest.fixture
    def data(self):
        file = "tests/input_files/example_result_table.csv"
        return pd.read_csv(file)

    def test_generate_csv_report(self, data):
        output_dir = "tests/output_files"

        report_generator = ReportGenerator(data)
        report_generator.generate_report(output_path=output_dir)

        assert os.path.exists(os.path.join(output_dir, "across_metrics_df.csv"))
        assert os.path.exists(os.path.join(output_dir, "single_metric_df.csv"))

    def test_generate_html_report(self, data):
        output_dir = "tests/output_files"

        report_generator = ReportGenerator(data)
        report_generator.generate_report(output_path=output_dir)

        assert os.path.exists(os.path.join(output_dir, "across_metrics_df.html"))
        assert os.path.exists(os.path.join(output_dir, "single_metric_df.html"))

    def test_generate_terminal_report(self, data, capsys):
        report_generator = ReportGenerator(data)
        report_generator.evaluate_data()
        report_generator._generate_terminal_report()

        output = capsys.readouterr().out

        assert "Category distribution" in output
        assert "number of evaluations" in output.lower()
        assert "evaluation summary" in output.lower()
        assert "good" in output.lower()
        assert "ok" in output.lower()
        assert "bad" in output.lower()
        assert "good" in report_generator.across_metrics_df.columns
        assert "Real Space R-Factor Evaluation" in report_generator.single_metric_df.index
        assert " (10.00%)" in output
    