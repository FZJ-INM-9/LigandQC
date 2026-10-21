"""
report
======

Report generation module for the ligand inspector. This module provides functionality to generate reports based on the analysis of ligand quality.

"""

import os
from pathlib import Path

import click
import pandas as pd
from rich.console import Console
from rich.table import Table

import ligandqc.const as const


class ReportGenerator:
    """
    A class to generate reports based on ligand quality analysis.

    Attributes:
        data (pd.DataFrame): The input data for report generation.
    """

    def __init__(self, data: pd.DataFrame):
        """
        Initializes the ReportGenerator with the provided data.

        Args:
            data (pd.DataFrame): The input data for report generation.
        """
        self.input_data = data
        self.single_metric_df = None
        self.across_metrics_df = None

    def generate_report(self, output_path: str = "."):
        """
        Generates a report based on the input data and saves it to the specified output path.

        Args:
            output_path (str): Directory for default table exports or a legacy file path ending in .csv / .html.

        Returns:
            None
        """
        output_path = os.fspath(output_path)
        if output_path.lower().endswith(".csv"):
            directory = os.path.dirname(output_path) or "."
            self._generate_csv_report(directory)
        elif output_path.lower().endswith(".html"):
            directory = os.path.dirname(output_path) or "."
            self._generate_html_report(directory)
        else:
            directory = output_path or "."
            self._generate_csv_report(directory)
            self._generate_html_report(directory)

    def _generate_csv_report(self, output_dir: str = "."):
        """
        Saves the formatted summary tables as CSV files in the target directory.

        Args:
            output_dir (str): Directory to write the default report files into.

        Returns:
            None
        """
        self.evaluate_data()

        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)

        self.single_metric_df.to_csv(directory / "single_metric_df.csv", index=True)
        self.across_metrics_df.to_csv(directory / "across_metrics_df.csv", index=True)

    def _generate_html_report(self, output_dir: str = "."):
        """
        Saves the formatted summary tables as HTML files in the target directory.

        Args:
            output_dir (str): Directory to write the default report files into.
        Returns:
            None
        """
        self.evaluate_data()

        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)

        self.single_metric_df.to_html(directory / "single_metric_df.html", index=True)
        self.across_metrics_df.to_html(directory / "across_metrics_df.html", index=True)

    def evaluate_data(self):
        """
        Evaluates the input data and prepares it for report generation.

        This method can include any necessary data processing or analysis steps before generating the report.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: the across-metrics and single-metric summary tables.
        """
        categories = [const.CATEGORY_GOOD, const.CATEGORY_OK, const.CATEGORY_BAD]
        total = len(self.input_data)

        def format_count(value: int) -> str:
            percentage = (value / total * 100) if total else 0.0
            return f"{int(value)} ({percentage:.2f}%)"

        row_labels = list(range(len(const.EVAL_COLUMNS) + 1))
        row_labels += [f">={threshold}" for threshold in range(1, len(const.EVAL_COLUMNS) + 1)]
        row_labels += [f"<={threshold}" for threshold in range(1, len(const.EVAL_COLUMNS) + 1)]

        across_metrics_df = pd.DataFrame(index=row_labels, columns=categories, dtype=object)
        across_metrics_df.index.name = f"Number of {const.EVAL_TERM_PLURAL}"
        across_metrics_df.columns.name = "Category"

        for category in categories:
            is_category = self.input_data[const.EVAL_COLUMNS] == category
            num_category_per_row = is_category.sum(axis=1)

            for count in range(len(const.EVAL_COLUMNS) + 1):
                value = int((num_category_per_row == count).sum())
                across_metrics_df.loc[count, category] = format_count(value)

            for threshold in range(1, len(const.EVAL_COLUMNS) + 1):
                value_at_least = int((num_category_per_row >= threshold).sum())
                value_at_most = int((num_category_per_row <= threshold).sum())
                across_metrics_df.loc[f">={threshold}", category] = format_count(value_at_least)
                across_metrics_df.loc[f"<={threshold}", category] = format_count(value_at_most)

        single_metric_df = pd.DataFrame(index=const.EVAL_COLUMNS, columns=categories, dtype=object)
        for col in const.EVAL_COLUMNS:
            counts = self.input_data[col].value_counts().reindex(categories, fill_value=0).astype(int)
            for category in categories:
                count = int(counts.get(category, 0))
                single_metric_df.loc[col, category] = format_count(count)

        self.single_metric_df = single_metric_df
        self.across_metrics_df = across_metrics_df
        return self.across_metrics_df, self.single_metric_df

    def _generate_terminal_report(self):
        """
        Generates a terminal report by printing the evaluated data to the console.

        Returns:
            None
        """
        self.evaluate_data()

        console = Console()
        console.print()

        column_colors = const.COLORS_RICH

        category_table = Table(
            title=f"Category distribution (Number of {const.EVAL_TERM_PLURAL})",
            show_lines=True,
        )
        category_table.add_column(f"Number of {const.EVAL_TERM_PLURAL}", style="bold")
        for column in self.across_metrics_df.columns:
            category_table.add_column(
                str(column),
                justify="center",
                style=column_colors.get(str(column).lower(), "white"),
            )

        for row_label in self.across_metrics_df.index:
            row = [str(row_label)]
            for value in self.across_metrics_df.loc[row_label].tolist():
                row.append(str(value))
            category_table.add_row(*row)

        console.print(category_table)
        console.print()

        summary_table = Table(title=f"{const.EVAL_TERM} summary (counts and percentages)", show_lines=True)
        summary_table.add_column(const.EVAL_TERM, style="bold")
        for column in self.single_metric_df.columns:
            summary_table.add_column(
                str(column),
                justify="center",
                style=column_colors.get(str(column).lower(), "white"),
            )

        for evaluation in self.single_metric_df.index:
            row = [str(evaluation)]
            for value in self.single_metric_df.loc[evaluation].tolist():
                row.append(str(value))
            summary_table.add_row(*row)

        console.print(summary_table)


@click.command()
@click.option(
    "--output",
    "-o",
    required=True,
    help='Path to save the report (CSV or HTML), alternatively "terminal" for terminal output.',
)
@click.option(
    "--input",
    "-i",
    type=click.Path(exists=True),
    required=True,
    help="Path to the input CSV data file.",
)
def report(output, input):
    """
    generate a report based on the input data and save it to the specified output path.
    """

    data = pd.read_csv(input)
    report_generator = ReportGenerator(data)

    if output == "terminal":
        report_generator._generate_terminal_report()
    else:
        report_generator.generate_report(output_path=output)
