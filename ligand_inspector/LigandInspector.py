"""
LigandInspector
===============

This class is responsible for inspecting the quality of the specified ligand and PDB id.

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
from rich.live import Live
from rich.spinner import Spinner
from rich.panel import Panel
from rich.table import Table
from rich.console import Group
from contextlib import suppress
import matplotlib.pyplot as plt
import seaborn as sns
from .RMSZInterpreter import RMSZInterpreter
from .LigandDownloader import LigandDownloader
from .Logger import logger
from .RCSBInterface import RCSBInterface
import json
from collections import defaultdict
from ligand_inspector.plots import ParetoPlot

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

        Args:
            rsr_factor_threshold_good (float): The threshold for the real space R factor, below which the ligand is considered good, default is 0.4. Everything above is considered bad.
            rscc_threshold_good (float): The threshold for the real space correlation coefficient, above which the ligand is considered good, default is 0.95.
            rscc_threshold_ok (float): The threshold for the real space correlation coefficient, above which the ligand is considered ok, default is 0.8. Everything below is considered bad.
            structural_outliers_threshold_good (int): The threshold for the number of Structural outliers, below which the ligand is considered good, default is 5.
            structural_outliers_threshold_ok (int): The threshold for the number of Structural outliers, above which the ligand is considered ok, default is 15. Everything above is considered bad.
            stereochemical_errors_threshold_good (int): The threshold for the number of stereochemical errors, below which the ligand is considered good, default is 0.
            average_occupancy_threshold_good (float): The threshold for the average occupancy, above which the ligand is considered good, default is 0.95.
            average_occupancy_threshold_ok (float): The threshold for the average occupancy, above which the ligand is considered ok, default is 0.5. Everything below is considered bad.
            atomic_clashes_threshold_good (int): The threshold for the number of atomic clashes, below which the ligand is considered good, default is 0.
            atomic_clashes_threshold_ok (int): The threshold for the number of atomic clashes, below which the ligand is considered ok, default is 2. Everything above is bad.

        """

        self.df = pd.DataFrame()
        # evaluation criteria
        self.rsr_factor_threshold_good = kwargs.get("rsr_factor_threshold_good", 0.4)
        self.rscc_threshold_good = kwargs.get("rscc_threshold_good", 0.95)
        self.rscc_threshold_ok = kwargs.get("rscc_threshold_ok", 0.8)
        self.structural_outliers_threshold_good = kwargs.get(
            "structural_outliers_threshold_good", 5
        )
        self.structural_outliers_threshold_ok = kwargs.get(
            "structural_outliers_threshold_ok", 15
        )
        self.stereochemical_errors_threshold_good = kwargs.get(
            "stereochemical_errors_threshold_good", 0
        )
        self.average_occupancy_threshold_good = kwargs.get(
            "average_occupancy_threshold_good", 0.95
        )
        self.average_occupancy_threshold_ok = kwargs.get(
            "average_occupancy_threshold_ok", 0.5
        )
        self.atomic_clashes_threshold_good = kwargs.get(
            "atomic_clashes_threshold_good", 0
        )
        self.atomic_clashes_threshold_ok = kwargs.get("atomic_clashes_threshold_ok", 2)

        self.ligand_downloader = LigandDownloader()
        self.rmsz_interpreter = RMSZInterpreter()
        logger.debug("Starting LigandInspector")

    def filter_df(self, df, only_chain=None, only_SOI=False):
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

        return df

    def inspect_batch(
        self,
        ligand_ids,
        pdb_ids=None,
        only_chain=None,
        intermediate_save=None,
        only_SOI=False,
    ):
        """
        wrapper function for async inspection of a batch of pdb-ligand pairs.

        Args:
            ligand_ids (list[str]): A list of ligand IDs.
            pdb_ids (list[str]): A list of PDB IDs.
            n_workers (int): The number of worker tasks to create.
            only_chain (str): If set, only the specified chain will be reported.
            intermediate_save (str): The path to the intermediate save file. If None, the default is "intermediate.json".
            only_SOI (bool): If True, only the entries that are subject of investigation will be kept in the output file.

        Returns:
            pd.DataFrame: A DataFrame containing the inspection results for each pdb-ligand pair.

        """
        interface = RCSBInterface()
        ligand_ids = list(set(ligand_ids))
        pdb_ids = list(set(pdb_ids)) if pdb_ids else None

        data = interface.get_ligand_quality_data(ligand_ids, pdb_ids=pdb_ids)

        # intermediate save
        if intermediate_save is None:
            intermediate_save = "intermediate.json"

        with open(intermediate_save, "w") as f:
            json.dump(data, f, indent=4)

        # process the data
        self.df = interface.process_ligand_quality_data(
            data, ligand_ids, pdb_ids=pdb_ids
        )

        self.df = self.filter_df(self.df, only_chain=only_chain, only_SOI=only_SOI)

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
        df["Real space correlation coefficient evaluation"] = df.apply(
            self.evaluate_rscc, axis=1
        )
        df["Real space R factor evaluation"] = df.apply(
            self.evaluate_rsr_factor, axis=1
        )
        df["Structural outliers evaluation"] = df.apply(
            self.evaluate_structural_outliers, axis=1
        )
        df["Stereochemical errors evaluation"] = df.apply(
            self.evaluate_stereochemical_errors, axis=1
        )
        df["Average occupancy evaluation"] = df.apply(self.evaluate_occupancy, axis=1)
        df["Atomic clashes evaluation"] = df.apply(self.evaluate_atomic_clashes, axis=1)
        df["RMSZ-bond-length evaluation"] = df.apply(
            self.evaluate_rmsz_bond_length, axis=1
        )
        df["RMSZ-bond-angle evaluation"] = df.apply(
            self.evaluate_rmsz_bond_angle, axis=1
        )
        df["overall evaluation"] = df.apply(
            lambda row: self.overall_evaluation(row), axis=1
        )

        return df

    def evaluate_rmsz_bond_length(self, row):
        """
        evaluates the RMSZ-bond-length of the ligand
        """
        return self.rmsz_interpreter.interpret_rmsz(
            row["Number of bonds"], row["RMSZ-bond-length"]
        )

    def evaluate_rmsz_bond_angle(self, row):
        """
        evaluates the RMSZ-bond-angle of the ligand
        """
        return self.rmsz_interpreter.interpret_rmsz(
            row["Number of angles"], row["RMSZ-bond-angle"]
        )

    def evaluate_rscc(self, row):
        """
        evaluates the RSCC value of the ligand
        """
        if row["Real space correlation coefficient"] >= self.rscc_threshold_good:
            return "good"
        elif row["Real space correlation coefficient"] >= self.rscc_threshold_ok:
            return "ok"
        else:
            return "bad"

    def evaluate_rsr_factor(self, row):
        """
        evaluates the RSR factor of the ligand
        """
        if row["Real space R factor"] <= self.rsr_factor_threshold_good:
            return "good"
        else:
            return "bad"

    def evaluate_structural_outliers(self, row):
        """
        evaluates the Structural outliers of the ligand
        """
        if row["Structural outliers"] <= self.structural_outliers_threshold_good:
            return "good"
        elif row["Structural outliers"] <= self.structural_outliers_threshold_ok:
            return "ok"
        else:
            return "bad"

    def evaluate_stereochemical_errors(self, row):
        """
        evaluates the stereochemical errors of the ligand
        """
        if row["Stereochemical errors"] <= self.stereochemical_errors_threshold_good:
            return "good"
        else:
            return "bad"

    def evaluate_occupancy(self, row):
        """
        evaluates the occupancy of the ligand
        """
        if row["Average occupancy"] >= self.average_occupancy_threshold_good:
            return "good"
        elif row["Average occupancy"] >= self.average_occupancy_threshold_ok:
            return "ok"
        else:
            return "bad"

    def evaluate_atomic_clashes(self, row):
        """
        evaluates the atomic clashes of the ligand
        """
        if row["Atomic clashes"] <= self.atomic_clashes_threshold_good:
            return "good"
        elif row["Atomic clashes"] <= self.atomic_clashes_threshold_ok:
            return "ok"
        else:
            return "bad"

    @staticmethod
    def overall_evaluation(row):
        """
        evaluates the overall quality of the ligand
        """
        evaluation_columns = [i for i in row.index if "evaluation" in i]

        if any([row[i] == "bad" for i in evaluation_columns]):
            return "bad"

        if any([row[i] == "ok" for i in evaluation_columns]):
            return "ok"

        return "good"

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

        colors = {"good": "green", "ok": "yellow", "bad": "orangered"}

        # Row 1: four distribution plots
        self.plot_RSR_factor(fig.add_subplot(first_row[0, 0]), df, colors)
        self.plot_RSCC(fig.add_subplot(first_row[0, 1]), df, colors)
        self.plot_RMSZ_bond_length(fig.add_subplot(first_row[0, 2]), df, colors)
        self.plot_RMSZ_bond_angle(fig.add_subplot(first_row[0, 3]), df, colors)

        # Row 2: one full-width horizontal stacked summary
        self.plot_evaluation_summary(fig.add_subplot(outer[1]), df, colors)

        # Row 3: pareto plots
        self.plot_pareto(
            fig.add_subplot(third_row[0, 0]), df, "bad", colors, reverse=True
        )
        self.plot_pareto(
            fig.add_subplot(third_row[0, 1]), df, "ok", colors, reverse=True
        )
        self.plot_pareto(
            fig.add_subplot(third_row[0, 2]), df, "good", colors, reverse=True
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

        plt.tight_layout()
        plt.savefig(filename)

    def plot_evaluation_summary(self, ax, df, colors):
        """
        Creates a horizontal stacked bar plot with the percentage of
        good/ok/bad evaluations for selected metrics.
        """
        evaluation_columns = [
            ("Structural outliers evaluation", "Structural outliers"),
            ("Atomic clashes evaluation", "Atomic clashes"),
            ("Stereochemical errors evaluation", "Stereochemical errors"),
            ("Average occupancy evaluation", "Average occupancy"),
        ]
        threshold_labels = {
            "Structural outliers": {
                "good": "<= 5",
                "ok": "6-15",
                "bad": "> 15",
            },
            "Atomic clashes": {
                "good": "0",
                "ok": "1-2",
                "bad": "> 2",
            },
            "Stereochemical errors": {
                "good": "0",
                "ok": "N/A",
                "bad": "> 0",
            },
            "Average occupancy": {
                "good": ">= 0.95",
                "ok": "0.5-0.95",
                "bad": "< 0.5",
            },
        }
        statuses = ["good", "ok", "bad"]

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
                if width < 3:
                    continue
                center_x = left.iloc[idx] + width / 2
                threshold_text = threshold_labels[label][status]
                ax.text(
                    center_x,
                    idx,
                    f"{status} ({threshold_text})\n{width:.1f}%",
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

    def plot_RSR_factor(self, ax, df, colors):
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
            x="Real space R factor",
            binwidth=0.025,
            ax=ax,
            hue="Real space R factor evaluation",
            hue_order=["good", "ok", "bad"],
            multiple="stack",
            palette=colors,
            binrange=(0, 1),
        )
        ax.set_xlabel("Real space R factor")
        ax.set_ylabel("Count")
        ax.set_xlim([0, 1])
        ax.vlines(
            self.rsr_factor_threshold_good,
            ymin=0,
            ymax=50,
            color="r",
            linestyle="--",
            label="RSR factor threshold",
        )

        self.update_legend(
            ax, df, "Real space R factor evaluation", "RSR factor evaluation"
        )

    def update_legend(self, ax, df, key, new_title):
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
            for i in ["good", "ok", "bad"]
        }
        ax.get_legend().set_title(new_title)
        new_labels = {"good": f"good ({ratios["good"]}%)",
                      "ok": f"ok ({ratios["ok"]}%)",
                      "bad": f"bad ({ratios["bad"]}%)"}
        for text in ax.get_legend().get_texts():
            text.set_text(new_labels.get(text.get_text(), text.get_text()))

    
    def plot_RSCC(self, ax, df, colors):
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
            x="Real space correlation coefficient",
            binwidth=0.025,
            ax=ax,
            hue="Real space correlation coefficient evaluation",
            hue_order=["good", "ok", "bad"],
            multiple="stack",
            palette=colors,
            binrange=(0, 1),
        )
        ax.set_xlabel("Real space correlation coefficient")
        ax.set_ylabel("Count")
        ax.set_xlim([0, 1])
        ax.vlines(
            self.rscc_threshold_good,
            ymin=0,
            ymax=50,
            color="r",
            linestyle="--",
            label="RSCC threshold",
        )
        ax.vlines(
            self.rscc_threshold_ok,
            ymin=0,
            ymax=50,
            color="orange",
            linestyle="--",
            label="RSCC threshold",
        )

        self.update_legend(
            ax, df, "Real space correlation coefficient evaluation", "RSCC evaluation"
        )

    def plot_RMSZ_bond_length(self, ax, df, colors):
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
            x="RMSZ-bond-length",
            binwidth=0.2,
            binrange=(0, max(df["RMSZ-bond-length"])),
            ax=ax,
            hue="RMSZ-bond-length evaluation",
            hue_order=["good", "ok", "bad"],
            multiple="stack",
            palette=colors,
        )
        ax.set_xlabel("RMSZ-bond-length")
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=5)

        self.update_legend(
            ax, df, "RMSZ-bond-length evaluation", "RMSZ-bond-length evaluation"
        )

    def plot_RMSZ_bond_angle(self, ax, df, colors):
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
            x="RMSZ-bond-angle",
            binwidth=0.2,
            binrange=(0, max(df["RMSZ-bond-angle"])),
            ax=ax,
            hue="RMSZ-bond-angle evaluation",
            hue_order=["good", "ok", "bad"],
            multiple="stack",
            palette=colors,
        )
        ax.set_xlabel("RMSZ-bond-angle")
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=5)
        self.update_legend(
            ax, df, "RMSZ-bond-angle evaluation", "RMSZ-bond-angle evaluation"
        )

    def plot_structural_outliers(self, ax, df, colors):
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
            x="Structural outliers",
            binwidth=2.5,
            binrange=(0, max(df["Structural outliers"])),
            ax=ax,
            hue="Structural outliers evaluation",
            multiple="stack",
            palette=colors,
        )
        ax.set_xlabel("Structural outliers")
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=30)
        self.update_legend(
            ax, df, "Structural outliers evaluation", "Structural outliers evaluation"
        )

    def plot_atomic_clashes(self, ax, df, colors, order):
        """
        creates a pie plot for atomic clashes

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        pie_data = df.groupby("Atomic clashes evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(
            lambda x: colors[x.split(" ")[0]]
        )
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (0)", "ok": "ok (1-2)", "bad": "bad (>2)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(
            pie_data_df["count"],
            labels=pie_data_df.index,
            autopct="%1.1f%%",
            startangle=90,
            colors=pie_data_df["color"],
            counterclock=False,
        )
        ax.set_xlabel("Atomic clashes")

    def plot_stereochemical_errors(self, ax, df, colors, order):
        """
        creates a pie plot for stereochemical errors

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        pie_data = df.groupby("Stereochemical errors evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(
            lambda x: colors[x.split(" ")[0]]
        )
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (0)", "bad": "bad (>0)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(
            pie_data_df["count"],
            labels=pie_data_df.index,
            autopct="%1.1f%%",
            startangle=90,
            colors=pie_data_df["color"],
            counterclock=False,
        )
        ax.set_xlabel("Stereochemical errors")

    def plot_average_occupancy(self, ax, df, colors, order):
        """
        creates a pie plot for average occupancy

        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """

        # pie plot for the average occupancy
        pie_data = df.groupby("Average occupancy evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(
            lambda x: colors[x.split(" ")[0]]
        )
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {
            "good": "good (>0.95)",
            "ok": "ok (0.5-0.95)",
            "bad": "bad (<0.5)",
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
        ax.set_xlabel("Average occupancy evaluation")

    def plot_pareto(self, ax, df, category, colors, reverse=False):
        """
        creates a nested donut plot showing the percentage of occurrences of the numbers of good, ok and bad evaluations in the dataset.
        The inner ring shows the number of good evaluations, the middle ring shows the number of ok evaluations and the outer ring shows the number of bad evaluations.

        Args:
            ax (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.

        Returns:
            None

        """
        evaluation_columns = [
            i for i in df.columns if "evaluation" in i and "overall" not in i
        ]
        pie_data_df = df.copy()
        pie_data_df.loc[:, f"Count of '{category.capitalize()}' evaluations"] = (
            pie_data_df[evaluation_columns].apply(
                lambda x: (x == category).sum(), axis=1
            )
        )

        p = ParetoPlot(pie_data_df)
        p.plot(
            ax,
            f"Count of '{category.capitalize()}' evaluations",
            sort_by_category=True,
            ascending=not reverse,
        )
