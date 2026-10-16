"""
inspect
=======

This module provides functionality to inspect the quality of ligands in PDB files. It includes the `LigandInspector` class, which is responsible for evaluating various metrics related to ligand quality.

Its functionality is exposed to the command line via the `inspect` and `plot` subcommands, which take a file containing ligand IDs and an output directory as input, and generate inspection results and plots for the specified ligands.

"""

import pandas as pd
from rich.progress import (
    Progress,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    MofNCompleteColumn,
    TimeRemainingColumn,
    track,
)
import matplotlib.pyplot as plt
import seaborn as sns
import click
import yaml
import os
import tempfile
import shutil
from pathlib import Path
from ligandqc.const import *
from ligandqc.interpret import RMSZInterpreter
from ligandqc.download import LigandDownloader
from ligandqc.log import logger
from ligandqc.rcsb import RCSBInterface
from ligandqc.plot import ParetoPlot
from ligandqc.log import *

global CURRENT_ITEMS
CURRENT_ITEMS = []
PROGRESS = Progress(
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TaskProgressColumn(),
    MofNCompleteColumn(),
    TimeRemainingColumn(),
)


class LigandInspector:
    def __init__(self, **kwargs):
        """
        Initializes the LigandInspector class.

        """

        self.df = pd.DataFrame()
        self.ligand_downloader = LigandDownloader()
        self.rmsz_interpreter = RMSZInterpreter()
        logger.debug("Starting LigandInspector")

    def _filter_df(self, df, only_chain=None, only_SOI=False):
        """
        filters the DataFrame as specified by the input parameters.

        Args:
            df (pd.DataFrame): The DataFrame to be filtered.
            only_chain (str): If set, only the specified chain will be reported.
            only_SOI (bool): If True, only the entries that are subject of investigation will be kept in the output file.

        Returns:
            pd.DataFrame: The filtered DataFrame.

        """

        if only_chain is not None:
            # remove entries that do not match the chain
            df = df[df["Chain ID"] == only_chain]

        if only_SOI:
            # remove entries that are not subject of investigation
            df = df[df["Is Subject of Investigation"] == True]

        # reset the index
        df = df.reset_index(drop=True)
        df.sort_values(by=["PDB ID", "Ligand ID", "Chain ID"], inplace=True)

        return df

    def inspect(
        self,
        ligand_ids=None,
        pdb_ids=None,
        only_chain=None,
        only_SOI=False,
    ):
        """
        wrapper function for async inspection of a batch of pdb-ligand pairs.

        Args:
            ligand_ids (list[str]): A list of ligand IDs.
            pdb_ids (list[str]): A list of PDB IDs.
            only_chain (str): If set, only the specified chain will be reported.
            only_SOI (bool): If True, only the entries that are subject of investigation will be kept in the output file.

        Returns:
            pd.DataFrame: A DataFrame containing the inspection results for each pdb-ligand pair.

        """
        interface = RCSBInterface()
        ligand_ids = list(set(ligand_ids)) if ligand_ids else None
        pdb_ids = list(set(pdb_ids)) if pdb_ids else None

        self.df = interface.get_ligand_quality_data(ligand_ids=ligand_ids, pdb_ids=pdb_ids)

        self.df = self._filter_df(self.df, only_chain=only_chain, only_SOI=only_SOI)

        return self.df

    def evaluate(self, df, ignore_other_entries=False, only_chain=None):
        """
        evaluates the quality of the ligand using the inspection results based on the threshold values defined in the constructor.

        Args:
            df (pd.DataFrame): A DataFrame containing the inspection results.
            ignore_other_entries (bool): If True, only the entries in the input file will be kept in the output file.
            only_chain (str): If set, only the specified chain will be reported.

        Returns:
            pd.DataFrame: A DataFrame containing the evaluation results.
        """

        # evaluate the ligand
        df[f"{METRIC_RSCC}{EVAL_SUFFIX}"] = df.apply(
            self._evaluate_rscc, axis=1
        )
        df[f"{METRIC_RSR}{EVAL_SUFFIX}"] = df.apply(
            self._evaluate_rsr_factor, axis=1
        )
        df[f"{METRIC_STRUCTURAL_OUTLIERS}{EVAL_SUFFIX}"] = df.apply(
            self._evaluate_structural_outliers, axis=1
        )
        df[f"{METRIC_STEREO_CHEMICAL_ERRORS}{EVAL_SUFFIX}"] = df.apply(
            self._evaluate_stereochemical_errors, axis=1
        )
        df[f"{METRIC_AVG_OCCUPANCY}{EVAL_SUFFIX}"] = df.apply(self._evaluate_occupancy, axis=1)
        df[f"{METRIC_ATOMIC_CLASHES}{EVAL_SUFFIX}"] = df.apply(self._evaluate_atomic_clashes, axis=1)
        df[f"{METRIC_RMSZ_BOND_LENGTH}{EVAL_SUFFIX}"] = df.apply(
            self._evaluate_rmsz_bond_length, axis=1
        )
        df[f"{METRIC_RMSZ_BOND_ANGLE}{EVAL_SUFFIX}"] = df.apply(
            self._evaluate_rmsz_bond_angle, axis=1
        )

        return df

    def _evaluate_rmsz_bond_length(self, row):
        """
        evaluates the RMSZ-bond-length of the ligand
        """
        return self.rmsz_interpreter.interpret_rmsz(
            row["Number of bonds"], row[METRIC_RMSZ_BOND_LENGTH]
        )

    def _evaluate_rmsz_bond_angle(self, row):
        """
        evaluates the RMSZ-bond-angle of the ligand
        """
        return self.rmsz_interpreter.interpret_rmsz(
            row["Number of angles"], row[METRIC_RMSZ_BOND_ANGLE]
        )

    def _evaluate_rscc(self, row):
        """
        evaluates the RSCC value of the ligand
        """
        if is_rscc_good(row[METRIC_RSCC]):
            return CATEGORY_GOOD
        if is_rscc_ok(row[METRIC_RSCC]):
            return CATEGORY_OK
        if is_rscc_bad(row[METRIC_RSCC]):
            return CATEGORY_BAD
        if pd.isna(row[METRIC_RSCC]):
            logger.warning(f"{METRIC_RSCC} value is NaN for PDB ID {row['PDB ID']} and ligand {row['Ligand ID']}. Marking as unknown.")
            return CATEGORY_UNKNOWN

        raise ValueError(f"{METRIC_RSCC} value {row[METRIC_RSCC]} is not valid for evaluation.")

    def _evaluate_rsr_factor(self, row):
        """
        evaluates the RSR factor of the ligand
        """
        if is_rsr_good(row[METRIC_RSR]):
            return CATEGORY_GOOD
        if is_rsr_bad(row[METRIC_RSR]):
            return CATEGORY_BAD
        if pd.isna(row[METRIC_RSR]):
            logger.warning(f"{METRIC_RSR} value is NaN for PDB ID {row['PDB ID']} and ligand {row['Ligand ID']}. Marking as unknown.")
            return CATEGORY_UNKNOWN

        raise ValueError(f"{METRIC_RSR} value {row[METRIC_RSR]} is not valid for evaluation.")

    def _evaluate_structural_outliers(self, row):
        """
        evaluates the Structural outliers of the ligand
        """
        if is_structural_outliers_good(row[METRIC_STRUCTURAL_OUTLIERS]):
            return CATEGORY_GOOD
        if is_structural_outliers_ok(row[METRIC_STRUCTURAL_OUTLIERS]):
            return CATEGORY_OK
        if is_structural_outliers_bad(row[METRIC_STRUCTURAL_OUTLIERS]):
            return CATEGORY_BAD
        if pd.isna(row[METRIC_STRUCTURAL_OUTLIERS]):
            logger.warning(f"{METRIC_STRUCTURAL_OUTLIERS} value is NaN for PDB ID {row['PDB ID']} and ligand {row['Ligand ID']}. Marking as unknown.")
            return CATEGORY_UNKNOWN
        
        raise ValueError(f"{METRIC_STRUCTURAL_OUTLIERS} value {row[METRIC_STRUCTURAL_OUTLIERS]} is not valid for evaluation.")

    def _evaluate_stereochemical_errors(self, row):
        """
        evaluates the stereochemical errors of the ligand
        """
        if is_stereo_chemical_errors_good(row[METRIC_STEREO_CHEMICAL_ERRORS]):
            return CATEGORY_GOOD
        if is_stereo_chemical_errors_bad(row[METRIC_STEREO_CHEMICAL_ERRORS]):
            return CATEGORY_BAD
        if pd.isna(row[METRIC_STEREO_CHEMICAL_ERRORS]):
            logger.warning(f"{METRIC_STEREO_CHEMICAL_ERRORS} value is NaN for PDB ID {row['PDB ID']} and ligand {row['Ligand ID']}. Marking as unknown.")
            return CATEGORY_UNKNOWN

        raise ValueError(f"{METRIC_STEREO_CHEMICAL_ERRORS} value {row[METRIC_STEREO_CHEMICAL_ERRORS]} is not valid for evaluation.")

    def _evaluate_occupancy(self, row):
        """
        evaluates the occupancy of the ligand
        """
        if is_avg_occupancy_good(row[METRIC_AVG_OCCUPANCY]):
            return CATEGORY_GOOD
        if is_avg_occupancy_ok(row[METRIC_AVG_OCCUPANCY]):
            return CATEGORY_OK
        if is_avg_occupancy_bad(row[METRIC_AVG_OCCUPANCY]):
            return CATEGORY_BAD
        if pd.isna(row[METRIC_AVG_OCCUPANCY]):
            logger.warning(f"{METRIC_AVG_OCCUPANCY} value is NaN for PDB ID {row['PDB ID']} and ligand {row['Ligand ID']}. Marking as unknown.")
            return CATEGORY_UNKNOWN

        raise ValueError(f"{METRIC_AVG_OCCUPANCY} value {row[METRIC_AVG_OCCUPANCY]} is not valid for evaluation.")

    def _evaluate_atomic_clashes(self, row):
        """
        evaluates the atomic clashes of the ligand
        """
        if is_atomic_clashes_good(row[METRIC_ATOMIC_CLASHES]):
            return CATEGORY_GOOD
        if is_atomic_clashes_ok(row[METRIC_ATOMIC_CLASHES]):
            return CATEGORY_OK
        if is_atomic_clashes_bad(row[METRIC_ATOMIC_CLASHES]):
            return CATEGORY_BAD
        if pd.isna(row[METRIC_ATOMIC_CLASHES]):
            logger.warning(f"{METRIC_ATOMIC_CLASHES} value is NaN for PDB ID {row['PDB ID']} and ligand {row['Ligand ID']}. Marking as unknown.")
            return CATEGORY_UNKNOWN

        raise ValueError(f"{METRIC_ATOMIC_CLASHES} value {row[METRIC_ATOMIC_CLASHES]} is not valid for evaluation.")

    def plot(self, df, filename):
        """
        creates plots of the inspection results using seaborn and matplotlib.

        Args:
            df (pd.DataFrame): A DataFrame containing the inspection results.
            filename (str): The filename to save the plot.

        Returns:
            None

        """

        df.dropna(inplace=True)

        fig = plt.figure(figsize=(18, 14))
        outer = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.9, 1.0], hspace=0.5)

        first_row = outer[0].subgridspec(1, 4, wspace=0.3)
        third_row = outer[2].subgridspec(1, 3, wspace=0.35)

        colors = COLORS_MATPLOTLIB

        # Row 1: four distribution plots
        self._plot_RSR_factor(fig.add_subplot(first_row[0, 0]), df, colors)
        self._plot_RSCC(fig.add_subplot(first_row[0, 1]), df, colors)
        self._plot_RMSZ_bond_length(fig.add_subplot(first_row[0, 2]), df, colors)
        self._plot_RMSZ_bond_angle(fig.add_subplot(first_row[0, 3]), df, colors)

        # Row 2: one full-width horizontal stacked summary
        self._plot_evaluation_summary(fig.add_subplot(outer[1]), df, colors)

        # Row 3: pareto plots
        self._plot_pareto(
            fig.add_subplot(third_row[0, 0]), df, CATEGORY_BAD, colors, reverse=True
        )
        self._plot_pareto(
            fig.add_subplot(third_row[0, 1]), df, CATEGORY_OK, colors, reverse=True
        )
        self._plot_pareto(
            fig.add_subplot(third_row[0, 2]), df, CATEGORY_GOOD, colors, reverse=True
        )

        # # plot the histograms
        # sns.histplot(df, x="Ranking for goodness of fit", bins=n_bins,  ax=axs[3,0])
        # axs[3,0].set_xlim([0, 1])

        # sns.histplot(df, x="Ranking for geometry", bins=n_bins,  ax=axs[3,1])
        # axs[3,1].set_xlim([0, 1])

        # # plot the occurence of the ligand in the PDB as a pie chart
        # pie_data = df.groupby("Ligand occurence in PDB").size()
        # pie_data_df = pd.DataFrame(pie_data)
        # pie_data_df = pie_data_df.rename(columns={0: "count"})
        # pie_data_df.sort_index(inplace=True)
        # axs[3,2].pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, counterclock=False)
        # axs[3,2].set_xlabel("Ligand occurence in PDB")

        plt.savefig(filename)

    def _plot_evaluation_summary(self, ax, df, colors):
        """
        Creates a horizontal stacked bar plot with the percentage of
        good/ok/bad evaluations for selected metrics.
        """
        evaluation_columns = [
            (f"{METRIC_STRUCTURAL_OUTLIERS}{EVAL_SUFFIX}", METRIC_STRUCTURAL_OUTLIERS),
            (f"{METRIC_ATOMIC_CLASHES}{EVAL_SUFFIX}", METRIC_ATOMIC_CLASHES),
            (f"{METRIC_STEREO_CHEMICAL_ERRORS}{EVAL_SUFFIX}", METRIC_STEREO_CHEMICAL_ERRORS),
            (f"{METRIC_AVG_OCCUPANCY}{EVAL_SUFFIX}", METRIC_AVG_OCCUPANCY),
        ]
        threshold_labels = {
            METRIC_STRUCTURAL_OUTLIERS: {
                CATEGORY_GOOD: f"<= {THRESH_STRUCTURAL_OUTLIERS_GOOD}",
                CATEGORY_OK: f"{THRESH_STRUCTURAL_OUTLIERS_GOOD + 1}-{THRESH_STRUCTURAL_OUTLIERS_BAD}",
                CATEGORY_BAD: f"> {THRESH_STRUCTURAL_OUTLIERS_BAD}",
            },
            METRIC_ATOMIC_CLASHES: {
                CATEGORY_GOOD: f"{THRESH_ATOMIC_CLASHES_GOOD}",
                CATEGORY_OK: f"{THRESH_ATOMIC_CLASHES_GOOD + 1}-{THRESH_ATOMIC_CLASHES_BAD}",
                CATEGORY_BAD: f"> {THRESH_ATOMIC_CLASHES_BAD}",
            },
            METRIC_STEREO_CHEMICAL_ERRORS: {
                CATEGORY_GOOD: f"{THRESH_STEREO_CHEMICAL_ERRORS_GOOD}",
                CATEGORY_OK: "N/A",
                CATEGORY_BAD: f"> {THRESH_STEREO_CHEMICAL_ERRORS_BAD}",
            },
            METRIC_AVG_OCCUPANCY: {
                CATEGORY_GOOD: f">= {THRESH_AVG_OCCUPANCY_GOOD}",
                CATEGORY_OK: f"{THRESH_AVG_OCCUPANCY_BAD}-{THRESH_AVG_OCCUPANCY_GOOD}",
                CATEGORY_BAD: f"< {THRESH_AVG_OCCUPANCY_BAD}",
            },
        }
        statuses = [CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD]

        summary_df = pd.DataFrame(
            index=[label for _, label in evaluation_columns], columns=statuses
        )

        for col, label in evaluation_columns:
            ratios = (
                df[col].value_counts(normalize=True).reindex(statuses, fill_value=0)
                * 100
            )
            summary_df.loc[label, statuses] = ratios.values

        summary_df = summary_df.astype(float)

        left = pd.Series(0.0, index=summary_df.index)
        for status in statuses:
            values = summary_df[status]
            ax.barh(
                summary_df.index,
                values,
                left=left,
                color=colors[status],
                alpha=0.85,
                label=status,
            )

            # Put percentage and threshold directly in each visible segment.
            for idx, label in enumerate(summary_df.index):
                width = values.iloc[idx]
                if width < 5:
                    continue
                center_x = left.iloc[idx] + width / 2
                threshold_text = threshold_labels[label][status]
                ax.text(
                    center_x,
                    idx,
                    f"{status} ({width:.1f}%)",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                )

            left = left + values

        ax.set_xlim(0, 100)
        ax.set_xlabel("Percentage (%)")
        ax.set_ylabel("Metric")
        ax.set_title("Evaluation Summary")
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)

    def _plot_RSR_factor(self, ax, df, colors):
        """
        creates a histogram of the RSR factor values.

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        sns.histplot(
            df,
            x=METRIC_RSR,
            binwidth=0.025,
            ax=ax,
            hue=f"{METRIC_RSR}{EVAL_SUFFIX}",
            hue_order=[CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD],
            multiple="stack",
            palette=colors,
            binrange=(0, 1),
        )
        ax.set_xlabel(METRIC_RSR)
        ax.set_ylabel("Count")
        ax.set_xlim([0, 1])
        # ax.axvline(
        #     THRESH_RSR_GOOD,
        #     color=COLORS_MATPLOTLIB[CATEGORY_GOOD],
        #     linestyle="--",
        #     label=f"{METRIC_RSR_SHORT} threshold",
        # )

        self._update_legend(
            ax, df, f"{METRIC_RSR}{EVAL_SUFFIX}", f"{METRIC_RSR_SHORT} {EVAL_TERM}"
        )

    def _update_legend(self, ax, df, key, new_title):
        """
        updates the legend of the plot.

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
            key (str): The key to update the legend for.
            new_title (str): The new title for the legend.

        Returns:
            None

        """
        ratios = {
            i: round(len(df[df[key] == i]) / len(df) * 100, 1)
            for i in [CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD]
        }
        ax.get_legend().set_title(new_title)
        new_labels = {
            CATEGORY_GOOD: f"{CATEGORY_GOOD} ({ratios[CATEGORY_GOOD]}%)",
            CATEGORY_OK: f"{CATEGORY_OK} ({ratios[CATEGORY_OK]}%)",
            CATEGORY_BAD: f"{CATEGORY_BAD} ({ratios[CATEGORY_BAD]}%)",
        }
        for text in ax.get_legend().get_texts():
            text.set_text(new_labels.get(text.get_text(), text.get_text()))

    
    def _plot_RSCC(self, ax, df, colors):
        """
        creates a histogram of the RSCC values.

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """

        sns.histplot(
            df,
            x=METRIC_RSCC,
            binwidth=0.025,
            ax=ax,
            hue=f"{METRIC_RSCC}{EVAL_SUFFIX}",
            hue_order=[CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD],
            multiple="stack",
            palette=colors,
            binrange=(0, 1),
        )
        ax.set_xlabel(METRIC_RSCC)
        ax.set_ylabel("Count")
        ax.set_xlim([0, 1])
        # ax.axvline(
        #     THRESH_RSCC_GOOD,
        #     color=COLORS_MATPLOTLIB[CATEGORY_GOOD],
        #     linestyle="--",
        #     label=f"{METRIC_RSCC_SHORT} threshold",
        # )
        # ax.axvline(
        #     THRESH_RSCC_BAD,
        #     color=COLORS_MATPLOTLIB[CATEGORY_BAD],
        #     linestyle="--",
        #     label=f"{METRIC_RSCC_SHORT} threshold",
        # )

        self._update_legend(
            ax, df, f"{METRIC_RSCC}{EVAL_SUFFIX}", f"{EVAL_TERM}"
        )

    def _plot_RMSZ_bond_length(self, ax, df, colors):
        """
        creates a histogram of the RMSZ values.

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        sns.histplot(
            df,
            x=METRIC_RMSZ_BOND_LENGTH,
            binwidth=0.2,
            binrange=(0, max(df[METRIC_RMSZ_BOND_LENGTH])),
            ax=ax,
            hue=f"{METRIC_RMSZ_BOND_LENGTH}{EVAL_SUFFIX}",
            hue_order=[CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD],
            multiple="stack",
            palette=colors,
        )
        ax.set_xlabel(METRIC_RMSZ_BOND_LENGTH)
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=5)

        self._update_legend(
            ax, df, f"{METRIC_RMSZ_BOND_LENGTH}{EVAL_SUFFIX}", f"{EVAL_TERM}"
        )

    def _plot_RMSZ_bond_angle(self, ax, df, colors):
        """
        creates a histogram of the RMSZ values.

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        sns.histplot(
            df,
            x=METRIC_RMSZ_BOND_ANGLE,
            binwidth=0.2,
            binrange=(0, max(df[METRIC_RMSZ_BOND_ANGLE])),
            ax=ax,
            hue=f"{METRIC_RMSZ_BOND_ANGLE}{EVAL_SUFFIX}",
            hue_order=[CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD],
            multiple="stack",
            palette=colors,
        )
        ax.set_xlabel(METRIC_RMSZ_BOND_ANGLE)
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=5)
        self._update_legend(
            ax, df, f"{METRIC_RMSZ_BOND_ANGLE}{EVAL_SUFFIX}", f"{EVAL_TERM}"
        )

    def _plot_structural_outliers(self, ax, df, colors):
        """
        creates a histogram of the structural outliers values.

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        sns.histplot(
            df,
            x=METRIC_STRUCTURAL_OUTLIERS,
            binwidth=2.5,
            binrange=(0, max(df[METRIC_STRUCTURAL_OUTLIERS])),
            ax=ax,
            hue=f"{METRIC_STRUCTURAL_OUTLIERS}{EVAL_SUFFIX}",
            multiple="stack",
            palette=colors,
        )
        ax.set_xlabel(METRIC_STRUCTURAL_OUTLIERS)
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=30)
        self._update_legend(
            ax, df, f"{METRIC_STRUCTURAL_OUTLIERS}{EVAL_SUFFIX}", f"{EVAL_TERM}"
        )

    def _plot_atomic_clashes(self, ax, df, colors, order):
        """
        creates a pie plot for atomic clashes

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        pie_data = df.groupby(f"{METRIC_ATOMIC_CLASHES}{EVAL_SUFFIX}").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(
            lambda x: colors[x.split(" ")[0]]
        )
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {
            CATEGORY_GOOD: f"{CATEGORY_GOOD} ({THRESH_ATOMIC_CLASHES_GOOD})",
            CATEGORY_OK: f"{CATEGORY_OK} ({THRESH_ATOMIC_CLASHES_GOOD + 1}-{THRESH_ATOMIC_CLASHES_BAD})",
            CATEGORY_BAD: f"{CATEGORY_BAD} (>{THRESH_ATOMIC_CLASHES_BAD})",
        }
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(
            pie_data_df["count"],
            labels=pie_data_df.index,
            autopct="%1.1f%%",
            startangle=90,
            colors=pie_data_df["color"],
            counterclock=False,
        )
        ax.set_xlabel(METRIC_ATOMIC_CLASHES)

    def _plot_stereochemical_errors(self, ax, df, colors, order):
        """
        creates a pie plot for stereochemical errors

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        pie_data = df.groupby(f"{METRIC_STEREO_CHEMICAL_ERRORS}{EVAL_SUFFIX}").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(
            lambda x: colors[x.split(" ")[0]]
        )
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {
            CATEGORY_GOOD: f"{CATEGORY_GOOD} ({THRESH_STEREO_CHEMICAL_ERRORS_GOOD})",
            CATEGORY_BAD: f"{CATEGORY_BAD} (>{THRESH_STEREO_CHEMICAL_ERRORS_BAD})",
        }
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(
            pie_data_df["count"],
            labels=pie_data_df.index,
            autopct="%1.1f%%",
            startangle=90,
            colors=pie_data_df["color"],
            counterclock=False,
        )
        ax.set_xlabel(METRIC_STEREO_CHEMICAL_ERRORS)

    def _plot_average_occupancy(self, ax, df, colors, order):
        """
        creates a pie plot for average occupancy

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """

        # pie plot for the average occupancy
        pie_data = df.groupby(f"{METRIC_AVG_OCCUPANCY}{EVAL_SUFFIX}").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(
            lambda x: colors[x.split(" ")[0]]
        )
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {
            CATEGORY_GOOD: f"{CATEGORY_GOOD} (>={THRESH_AVG_OCCUPANCY_GOOD})",
            CATEGORY_OK: f"{CATEGORY_OK} ({THRESH_AVG_OCCUPANCY_BAD}-{THRESH_AVG_OCCUPANCY_GOOD})",
            CATEGORY_BAD: f"{CATEGORY_BAD} (<{THRESH_AVG_OCCUPANCY_BAD})",
        }
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(
            pie_data_df["count"],
            labels=pie_data_df.index,
            autopct="%1.1f%%",
            startangle=90,
            colors=pie_data_df["color"],
            counterclock=False,
        )
        ax.set_xlabel(METRIC_AVG_OCCUPANCY)

    def _plot_pareto(self, ax, df, category, colors, reverse=False):
        """
        creates a nested donut plot showing the percentage of occurrences of the numbers of good, ok and bad evaluations in the dataset.
        The inner ring shows the number of good evaluations, the middle ring shows the number of ok evaluations and the outer ring shows the number of bad evaluations.

        Args:
            ax (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        pie_data_df = df.copy()
        pie_data_df.loc[:, f"Count of '{category}' {EVAL_TERM_PLURAL}"] = (
            pie_data_df[EVAL_COLUMNS].apply(
                lambda x: (x == category).sum(), axis=1
            )
        )

        p = ParetoPlot(pie_data_df)
        p.plot(
            ax,
            f"Count of '{category}' {EVAL_TERM_PLURAL}",
            sort_by_category=True,
            ascending=not reverse,
            rotate_xticks=False,
        )



@click.command()
@click.option('--ligand-id-file', '-l',
              type=click.Path(),
              default=None,
              show_default=True,
              help="Path to the input file containing ligand IDs.")
@click.option('--output-file', '-o',
              type=click.Path(exists=False),
              required=True,
              help="Path to the output file where results will be saved.")
@click.option('--only-chain', "-oc",
                type=click.STRING,
                help="If set, only the specified chain will be inspected.")
@click.option("--only-SOI", "-os",
                is_flag=True,
                help="If set, only the ligands that are subject of investigation will be inspected.")
@click.option('--logging-level', '-ll',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], case_sensitive=False),
              default='INFO',
              help="Set the logging level; Choices are DEBUG, INFO, WARNING, ERROR, CRITICAL. Default is INFO.")
@click.option("--log-file", "-lf",
              type=click.Path(exists=False),
              default="ligandqc.log",
                help="Path to the log file. Default is ligandqc.log.")
@click.option("--pdb-constraint-file", "-pcf",
              type=click.Path(),
              default=None,
              help="Path to the PDB constraint file. If specified, the results will only show ligands in the specified PDB IDs. The file should contain one PDB ID per line.")

def inspect(ligand_id_file, output_file, only_chain, only_soi, logging_level, log_file, pdb_constraint_file):
    """
    Inspect a batch of ligands or PDB ids and reports their quality metrics and interpretations. You must provide either the ligand ID file or the PDB constraint file, or both. If none is provided, the inspection cannot proceed.
    """
    setup_logger(log_file, logging_level)
    logger.info(f"Inspecting ligands in {ligand_id_file}")

    if ligand_id_file is None:
        ligand_ids = []
    elif os.path.exists(ligand_id_file):
        logger.info(f"Reading input file {ligand_id_file}...")
        with open(ligand_id_file, 'r') as f:
            ligand_ids = f.read().splitlines()
            ligand_ids = [ligand_id.strip() for ligand_id in ligand_ids]
            while "" in ligand_ids:
                ligand_ids.remove("")
        if len(ligand_ids) == 0:
            logger.critical(f"Input file {ligand_id_file} is empty.")
            return
        logger.info(f"Found {len(ligand_ids)} ligands in input file.")

    if pdb_constraint_file:
        logger.info(f"Reading PDB constraint file {pdb_constraint_file}...")
        with open(pdb_constraint_file, 'r') as f:
            pdb_ids = f.read().splitlines()
            pdb_ids = [pdb_id.strip() for pdb_id in pdb_ids]
            while "" in pdb_ids:
                pdb_ids.remove("")
        if len(pdb_ids) == 0:
            logger.critical(f"PDB constraint file {pdb_constraint_file} is empty.")
            return
        logger.info(f"Found {len(pdb_ids)} PDB IDs in constraint file.")
    else:
        logger.info(f"No PDB constraint file provided. Using all PDB IDs.")
        pdb_ids = []

    li = LigandInspector()
    results = li.inspect(ligand_ids=ligand_ids, pdb_ids=pdb_ids, only_chain=only_chain,only_SOI=only_soi)
    
    logger.info(f"Evaluating the provided data...")
    results = li.evaluate(results)

    # Save the results to the specified output file
    logger.info(f"Saving results to {output_file}...")
    if output_file.endswith(".yaml") or output_file.endswith(".yml"):
        with open(output_file, 'w') as f:
            yaml.dump(results.to_dict(orient='records'), f)
    elif output_file.endswith(".json"):
        results.to_json(output_file, orient='records', lines=True)
    elif output_file.endswith(".xlsx"):
        results.to_excel(output_file, index=False)
    elif output_file.endswith(".html"):
        results.to_html(output_file, index=False)
    elif output_file.endswith(".txt"):
        results.to_csv(output_file, sep='\t', index=False)
    elif output_file.endswith(".csv"):
        results.to_csv(output_file, index=False)
    elif output_file.endswith(".md"):
        results.to_markdown(output_file, index=False)
    else:
        click.echo("Unsupported file format. Please use .csv, .xlsx, .json, .yaml, .yml, .html or .txt.")
        return


@click.command()
@click.option('--input-file', '-i', 
              type=click.Path(exists=True),
              required=True,
              help="Path to the input CSV file containing the ligand quality data (output of the inspect command).")
@click.option('--output-file', '-o', 
              required=True, 
              type=click.Path(exists=False),
              help="Path to the output file where the plot will be saved (e.g., .png, .jpg, .pdf).")
def plot(input_file, output_file):
    """
    Plot the results from the input CSV file and save the plot to the specified output file.
    """
    logger.info(f"Plotting results from {input_file}...")
    df = pd.read_csv(input_file)
    li = LigandInspector()
    li.plot(df, output_file)
    logger.info(f"Plot saved to {output_file}")

@click.group()
def main():
    """
    Ligand Inspector CLI tool.
    """
    pass

main.add_command(inspect)
main.add_command(plot)


if __name__ == "__main__":
    main()