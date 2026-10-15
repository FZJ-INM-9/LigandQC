"""
compare
=======

This module allows for the comparison of two datasets containing ligand quality metrics. 
It 

* creates a venn diagram to visualize the overlap of ligands between the two datasets.
* creates a comparative pareto plot of the metrics for the two datasets.
* creates a table comparing the metrics for the two datasets.
* export the comparisons to a folder in png, csv, and html formats.

"""


import os
import pandas as pd
import click
import matplotlib.pyplot as plt
from ligand_inspector.const import *
from ligand_inspector.plot import ComparativeParetoPlot, VennPlot
from ligand_inspector.log import *
from ligand_inspector.report import ReportGenerator

class DatasetComparison:
    """
    A class to compare multiple datasets containing ligand quality metrics.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes the DatasetComparison class with the provided datasets.

        Args:
            *args: Variable length argument list for datasets.
            **kwargs: Arbitrary keyword arguments for dataset labels.
        """
        self.datasets = args
        self.labels = kwargs.get("labels", [f"Dataset {i+1}" for i in range(len(args))])

        self.df = pd.concat(self.datasets, keys=self.labels, names=["Dataset", "Index"])
        self.df.reset_index(inplace=True)
        self.df.drop(columns=["Index"], inplace=True)

        for cat in CATEGORY_ORDER:
            self.df[f"Count of '{cat}' {EVAL_TERM_PLURAL}"] = self.df.apply(lambda s: sum([1 if v == cat else 0 for v in s]), axis=1)
        

    def compare_datasets(self, output_dir: str, img_format: str = "png", table_format: str = "csv", normalize: bool = False):
        """
        Compares the datasets and generates visualizations and a comparison table.

        Args:
            output_dir (str): Directory to save the comparison results.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Create a Venn diagram for the overlap of ligands
        self._create_venn_diagrams(output_dir, format=img_format)

        # Create a comparative Pareto plot for the metrics
        self._create_comparative_pareto_plots(output_dir, normalize=normalize, format=img_format)

        # Create a comparison table for the metrics
        self._create_comparison_tables(output_dir, format=table_format)

        # Export the comparison results to HTML
        self._export_to_html(output_dir)

    def _create_venn_diagrams(self, output_dir: str, format: str = "png"):
        """
        creates a Venn diagram to visualize the overlap of ligands between the datasets.

        Args:
            output_dir (str): Directory to save the Venn diagram.
            format (str): Format to save the Venn diagram (default is "png").

        Returns:
            None
        """

        # create a Venn diagram for the overlap of ligands
        venn_sets = [set(dataset['Ligand ID']) for dataset in self.datasets]
        venn_plot = VennPlot(venn_sets, self.labels)

        fig, ax = plt.subplots(figsize=(8, 6))
        venn_plot.plot(ax=ax)
        ax.set_title("Venn Diagram of Ligand ID Overlap")
        file = os.path.join(output_dir, f"venn_diagram_ligand_ids.{format}")
        plt.savefig(file)
        plt.close()
        logger.info(f"Saved Venn diagram to {file}")

        # create a Venn diagram for the overlap of 'PDB ID' and 'Ligand ID' pairs
        venn_sets_pairs = [set([f"{pdb_id}_{ligand_id}" for pdb_id, ligand_id in zip(dataset['PDB ID'], dataset['Ligand ID'])]) for dataset in self.datasets]
        venn_plot_pairs = VennPlot(venn_sets_pairs, self.labels)

        fig, ax = plt.subplots(figsize=(8, 6))
        venn_plot_pairs.plot(ax=ax)
        ax.set_title("Venn Diagram of PDB ID and Ligand ID Overlap")
        file = os.path.join(output_dir, f"venn_diagram_pdb_id_ligand_id_pairs.{format}")
        plt.savefig(file)
        plt.close()
        logger.info(f"Saved Venn diagram of PDB ID and Ligand ID pairs to {file}")



    def _create_comparative_pareto_plots(self, output_dir: str, normalize: bool = False, format: str = "png"):
        """
        Creates a comparative Pareto plot for the metrics in the datasets.

        Args:
            output_dir (str): Directory to save the Pareto plot.
            format (str): Format to save the Pareto plot (default is "png").

        Returns:
            None
        """
        pareto_plot = ComparativeParetoPlot(self.df, dataset_column="Dataset")

        for cat in CATEGORY_ORDER:
            fig, ax = plt.subplots(figsize=(10, 6))
            pareto_plot.plot(ax=ax, column=f"Count of '{cat}' {EVAL_TERM_PLURAL}", rotate_xticks=False, sort_by_category=True, ascending=False, normalize=normalize)
            file = os.path.join(output_dir, f"comparative_pareto_plot_{cat}.{format}")
            plt.savefig(file)
            plt.close()
            logger.info(f"Saved comparative Pareto plot for '{cat}' {EVAL_TERM_PLURAL} to {file}")

    def _create_comparison_tables(self, output_dir: str, format: str = "csv"):
        """
        Creates comparison tables by merging the metric tables from the reports module.

        Args:
            output_dir (str): Directory to save the comparison table.
            format (str): Format to save the comparison table (default is "csv").

        Returns:
            None
        """

        report_generators = [ReportGenerator(dataset) for dataset in self.datasets]
        for report_generator in report_generators:
            report_generator.evaluate_data()

        comparison_table = pd.concat([report_generator.single_metric_df for report_generator in report_generators], axis=1, keys=self.labels)
        comparison_table.index.name = "Metric"
        comparison_table.reset_index(inplace=True)
        file = os.path.join(output_dir, f"comparison_single_metric.{format}")
        comparison_table.to_csv(file, index=False)
        logger.info(f"Saved comparison table to {file}")

        comparison_table2 = pd.concat([report_generator.across_metrics_df for report_generator in report_generators], axis=1, keys=self.labels)
        comparison_table2.reset_index(inplace=True)
        file2 = os.path.join(output_dir, f"comparison_across_metrics.{format}")
        comparison_table2.to_csv(file2, index=False)
        logger.info(f"Saved comparison table across metrics to {file2}")

    def _export_to_html(self, output_dir: str):
        """
        Creates an HTML report for the comparison of the datasets, referencing the generated plots and tables.

        Args:
            output_dir (str): Directory to save the HTML report.
        
        Returns:
            None
        
        """
        html_content = f"""
        <html>
        <head>
            <title>Dataset Comparison Report</title>
        </head>
        <body>
            <h1>Dataset Comparison Report</h1>
            <h2>Venn Diagrams</h2>
            <img src="venn_diagram_ligand_ids.png" alt="Venn Diagram of Ligand ID Overlap" width="600">
            <img src="venn_diagram_pdb_id_ligand_id_pairs.png" alt="Venn Diagram of PDB ID and Ligand ID Overlap" width="600">
            <h2>Comparative Pareto Plots</h2>
            <img src="comparative_pareto_plot_good.png" alt="Comparative Pareto Plot for Good Metrics" width="600">
            <img src="comparative_pareto_plot_ok.png" alt="Comparative Pareto Plot for OK Metrics" width="600">
            <img src="comparative_pareto_plot_bad.png" alt="Comparative Pareto Plot for Bad Metrics" width="600">
            <h2>Comparison Tables</h2>
            <a href="comparison_single_metric.csv">Download Comparison Table (Single Metric)</a><br>
            <a href="comparison_across_metrics.csv">Download Comparison Table (Across Metrics)</a>
        </body>
        </html>
        """

        file = os.path.join(output_dir, "comparison.html")
        with open(file, "w") as f:
            f.write(html_content)
        logger.info(f"Saved HTML report to {file}")


@click.command()
@click.option('--output-dir', 
              default='comparison_results', 
              help='Directory to save the comparison results.',
              show_default=True)
@click.option('--labels',
              default=None,
              help='Comma-separated labels for the datasets. If not provided, default labels will be used.',
              show_default=True)
@click.option('--img-format',
              default='png',
              help='Format to save the plots.',
              show_default=True)
@click.option('--table-format',
              default='csv',
              help='Format to save the comparison tables.',
              show_default=True)
@click.option('--normalize', 
              is_flag=True, 
              help='Flag to normalize the counts in the comparative Pareto plot.')
@click.argument('input_files', nargs=-1)
def compare(input_files, output_dir, labels, img_format, table_format, normalize):
    """
    Compares multiple datasets containing ligand quality metrics.

    Args:
        input_files (tuple): Paths to the input CSV files containing the datasets.
        output_dir (str): Directory to save the comparison results.
        labels (str): Comma-separated labels for the datasets. If not provided, default labels will be used.
    """
    if len(input_files) < 2:
        click.echo("Please provide at least two input files for comparison.")
        return
    if len(input_files) > 3:
        click.echo("Currently, only up to three datasets can be compared.")
        return

    datasets = [pd.read_csv(file) for file in input_files]
    
    if labels:
        label_list = [label.strip() for label in labels.split(',')]
        if len(label_list) != len(datasets):
            click.echo("The number of labels must match the number of input files.")
            return
    else:
        label_list = None

    comparison = DatasetComparison(*datasets, labels=label_list)
    comparison.compare_datasets(output_dir, img_format=img_format, table_format=table_format, normalize=normalize)

if __name__ == "__main__":
    compare()