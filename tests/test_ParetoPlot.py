import pytest
import matplotlib.pyplot as plt
from ligand_inspector.plots import ParetoPlot
import pandas as pd
from random import choice

sample_data = []
for label, n in zip("ABCDEF", [10, 20, 5, 15, 25, 30]):
    sample_data += [{"Category": label}] * n


class TestParetoPlot:

    def test_plot_sorting_truncated(self):
        # Sample data
        df = pd.DataFrame(sample_data)
        print(df)

        # Create a ParetoPlot instance
        pareto_plot = ParetoPlot(df)

        # Create a matplotlib axis
        fig, ax = plt.subplots()

        # Call the plot method
        pareto_plot.plot(
            ax, column="Category", truncate=5, sort_by_category=True, ascending=True
        )

        plt.tight_layout()
        plt.savefig("tests/output_files/pareto_plot_test_sorted_truncated.svg")

    def test_plot_sorting(self):
        # Sample data
        df = pd.DataFrame(sample_data)
        print(df)

        # Create a ParetoPlot instance
        pareto_plot = ParetoPlot(df)

        # Create a matplotlib axis
        fig, ax = plt.subplots()

        # Call the plot method
        pareto_plot.plot(ax, column="Category", sort_by_category=True, ascending=True)

        plt.tight_layout()
        plt.savefig("tests/output_files/pareto_plot_test_sorted.svg")

    def test_plot_sorting_descending(self):
        # Sample data
        df = pd.DataFrame(sample_data)
        print(df)

        # Create a ParetoPlot instance
        pareto_plot = ParetoPlot(df)

        # Create a matplotlib axis
        fig, ax = plt.subplots()

        # Call the plot method
        pareto_plot.plot(ax, column="Category", sort_by_category=True, ascending=False)

        plt.tight_layout()
        plt.savefig("tests/output_files/pareto_plot_test_sorted_descending.svg")

    def test_plot_unsorted(self):
        # Sample data
        df = pd.DataFrame(sample_data)

        # Create a ParetoPlot instance
        pareto_plot = ParetoPlot(df)

        # Create a matplotlib axis
        fig, ax = plt.subplots()

        # Call the plot method
        pareto_plot.plot(ax, column="Category", truncate=5)

        plt.tight_layout()
        plt.savefig("tests/output_files/pareto_plot_test_unsorted.svg")
