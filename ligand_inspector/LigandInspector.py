"""
LigandInspector
===============

This class is responsible for inspecting the quality of the specified ligand and PDB id.

"""

import pandas as pd
from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, MofNCompleteColumn, TimeRemainingColumn, track
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
        self.structural_outliers_threshold_good = kwargs.get("structural_outliers_threshold_good", 5)
        self.structural_outliers_threshold_ok = kwargs.get("structural_outliers_threshold_ok", 15)
        self.stereochemical_errors_threshold_good = kwargs.get("stereochemical_errors_threshold_good", 0)
        self.average_occupancy_threshold_good = kwargs.get("average_occupancy_threshold_good", 0.95)
        self.average_occupancy_threshold_ok = kwargs.get("average_occupancy_threshold_ok", 0.5)
        self.atomic_clashes_threshold_good = kwargs.get("atomic_clashes_threshold_good", 0)
        self.atomic_clashes_threshold_ok = kwargs.get("atomic_clashes_threshold_ok", 2)

        self.ligand_downloader = LigandDownloader()
        self.rmsz_interpreter = RMSZInterpreter()
        logger.debug("Starting LigandInspector")


    
    def filter_df(self, df, only_chain=None):
        """
        filters the DataFrame as specified by the input parameters.

        Args:
            df (pd.DataFrame): The DataFrame to be filtered.
            only_chain (str): If set, only the specified chain will be reported.

        Returns:
            pd.DataFrame: The filtered DataFrame.
        
        """
        
        if only_chain is not None:
            # remove entries that do not match the chain
            df = df[df["Chain ID"]== only_chain]

        return df

    
    def inspect_batch(self, ligand_ids, pdb_ids=None, only_chain=None, intermediate_save=None):
        """
        wrapper function for async inspection of a batch of pdb-ligand pairs.

        Args:
            ligand_ids (list[str]): A list of ligand IDs.
            pdb_ids (list[str]): A list of PDB IDs.
            n_workers (int): The number of worker tasks to create.
            only_chain (str): If set, only the specified chain will be reported.
            intermediate_save (str): The path to the intermediate save file. If None, the default is "intermediate.json".
        
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
        self.df = interface.process_ligand_quality_data(data, ligand_ids, pdb_ids=pdb_ids)
        
        self.df = self.filter_df(self.df, only_chain=only_chain)

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
        df["Real space correlation coefficient evaluation"] = df.apply(self.evaluate_rscc, axis=1)
        df["Real space R factor evaluation"] = df.apply(self.evaluate_rsr_factor, axis=1)
        df["Structural outliers evaluation"] = df.apply(self.evaluate_structural_outliers, axis=1)
        df["Stereochemical errors evaluation"] = df.apply(self.evaluate_stereochemical_errors, axis=1)
        df["Average occupancy evaluation"] = df.apply(self.evaluate_occupancy, axis=1)
        df["Atomic clashes evaluation"] = df.apply(self.evaluate_atomic_clashes, axis=1)
        df["RMSZ-bond-length evaluation"] = df.apply(self.evaluate_rmsz_bond_length, axis=1)
        df["RMSZ-bond-angle evaluation"] = df.apply(self.evaluate_rmsz_bond_angle, axis=1)
        df["overall evaluation"] = df.apply(lambda row: self.overall_evaluation(row), axis=1)

        return df
    
    
    def evaluate_rmsz_bond_length(self, row):
        """
        evaluates the RMSZ-bond-length of the ligand
        """
        return self.rmsz_interpreter.interpret_rmsz(row["Number of bonds"], row["RMSZ-bond-length"])
    
    def evaluate_rmsz_bond_angle(self, row):
        """
        evaluates the RMSZ-bond-angle of the ligand
        """
        return self.rmsz_interpreter.interpret_rmsz(row["Number of angles"], row["RMSZ-bond-angle"])
    
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

        fig, axs = plt.subplots(3, 3, figsize=(14, 14))
        n_bins = len(df) // 10
        if n_bins < 1:
            n_bins = 1

        colors = {"good": "green", "ok": "yellow", "bad": "orangered"}
        order = {"good": 0, "ok": 1, "bad": 2}

        self.plot_RSR_factor(axs[0,0], df, colors)
        self.plot_RSCC(axs[0,1], df, colors)
        self.plot_RMSZ_bond_length(axs[0,2], df, colors)
        self.plot_RMSZ_bond_angle(axs[1,0], df, colors)
        self.plot_structural_outliers(axs[1,1], df, colors)
        self.plot_atomic_clashes(axs[1,2], df, colors, order)
        self.plot_stereochemical_errors(axs[2,0], df, colors, order)
        self.plot_average_occupancy(axs[2,1], df, colors, order)
        self.plot_overall_evaluation(axs[2,2], df, colors, order)
        

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
    
    def plot_RSR_factor(self, ax, df, colors):
        """
        creates a histogram of the RSR factor values.
        
        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
        
        Returns:
            None
        
        """
        sns.histplot(df, x="Real space R factor", binwidth=0.025,  ax=ax, hue="Real space R factor evaluation", multiple="stack", palette=colors, binrange=(0, 1))
        ax.set_xlabel("Real space R factor")
        ax.set_ylabel("Count")
        ax.set_xlim([0, 1])
        ax.vlines(self.rsr_factor_threshold_good, ymin=0, ymax=50, color='r', linestyle='--', label='RSR factor threshold')
    
    def plot_RSCC(self, ax, df, colors):
        """
        creates a histogram of the RSCC values.
        
        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
        
        Returns:
            None
        
        """
        sns.histplot(df, x="Real space correlation coefficient", binwidth=0.025,  ax=ax, hue="Real space correlation coefficient evaluation", multiple="stack", palette=colors, binrange=(0, 1))
        ax.set_xlabel("Real space correlation coefficient")
        ax.set_ylabel("Count")
        ax.set_xlim([0, 1])
        ax.vlines(self.rscc_threshold_good, ymin=0, ymax=50, color='r', linestyle='--', label='RSCC threshold')
        ax.vlines(self.rscc_threshold_ok, ymin=0, ymax=50, color='orange', linestyle='--', label='RSCC threshold')
    
    def plot_RMSZ_bond_length(self, ax, df, colors):
        """
        creates a histogram of the RMSZ values.
        
        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
        
        Returns:
            None
        
        """
        sns.histplot(df, x="RMSZ-bond-length", 
                     binwidth=0.2,
                     binrange=(0, max(df["RMSZ-bond-length"])),
                     ax=ax, 
                     hue="RMSZ-bond-length evaluation",
                     multiple="stack", 
                     palette=colors)
        ax.set_xlabel("RMSZ-bond-length")
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=5)

    def plot_RMSZ_bond_angle(self, ax, df, colors):
        """
        creates a histogram of the RMSZ values.
        
        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
        
        Returns:
            None
        
        """
        sns.histplot(df, x="RMSZ-bond-angle", 
                     binwidth=0.2,
                     binrange=(0, max(df["RMSZ-bond-angle"])),
                     ax=ax, 
                     hue="RMSZ-bond-angle evaluation",
                     multiple="stack", 
                     palette=colors)
        ax.set_xlabel("RMSZ-bond-angle")
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=5)

    def plot_structural_outliers(self, ax, df, colors):
        """
        creates a histogram of the structural outliers values.
        
        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
        
        Returns:
            None
        
        """
        sns.histplot(df, x="Structural outliers", 
                     binwidth=2.5,
                     binrange=(0, max(df["Structural outliers"])),
                     ax=ax, 
                     hue="Structural outliers evaluation",
                     multiple="stack", 
                     palette=colors)
        ax.set_xlabel("Structural outliers")
        ax.set_ylabel("Count")
        ax.set_xlim(xmin=0, xmax=30)
    
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
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (0)", 
                      "ok": "ok (1-2)", 
                      "bad": "bad (>2)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
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
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (0)", "bad": "bad (>0)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
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
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (>0.95)", "ok": "ok (0.5-0.95)", "bad": "bad (<0.5)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        ax.pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
        ax.set_xlabel("Average occupancy evaluation")

    def plot_overall_evaluation(self, ax, df, colors, order):
        """
        creates a pie plot for overall evaluation
        
        Args:
            axs (matplotlib.axes._axes.Axes): The axes to plot on.
            df (pd.DataFrame): A DataFrame containing the inspection results.
        
        Returns:
            None
        
        """
        pie_data = df.groupby("overall evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))

        ax.pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
        ax.set_xlabel("Overall evaluation")