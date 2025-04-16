"""
LigandInspector
===============

This class is responsible for inspecting the quality of the specified ligand and PDB id.

"""

import os
import aiohttp
import pandas as pd
import asyncio
from itertools import batched
from bs4 import BeautifulSoup
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
from rdkit import Chem
from .Logger import logger
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG, MolDrawOptions
from rdkit.Chem import AllChem

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


    def build_url(self, pdb_id, ligand_id):
        """
        Builds the URL for the ligand information.

        Args:
            pdb_id (str): The PDB ID of the protein.
            ligand_id (str): The ligand ID.

        Returns:
            str: The constructed URL.
        """
        return f"https://www.rcsb.org/ligand-validation/{pdb_id}/{ligand_id}"
    
    def inspect_single(self, pdb_id, ligand_id, ignore_other_entries=False, only_chain=None):
        """
        Inspects a single pdb-ligand pair.

        Args:
            pdb_id (str): The PDB ID of the protein.
            ligand_id (str): The ligand ID.
            ignore_other_entries (bool): If True, only the entries in the input file will be kept in the output file.
            only_chain (str): If set, only the specified chain will be reported.

        Returns:
            pd.DataFrame: A DataFrame containing the inspection results.

        """
        df = pd.DataFrame([{"PDB ID": pdb_id, "Ligand ID": ligand_id}])
        df = self.inspect_batch(df, n_workers=1, ignore_other_entries=ignore_other_entries, only_chain=only_chain)
        return df
    
    def filter_df(self, df, ignore_other_entries=False, only_chain=None):
        """
        filters the DataFrame as specified by the input parameters.

        Args:
            df (pd.DataFrame): The DataFrame to be filtered.
            ignore_other_entries (bool): If True, only the entries in the input file will be kept in the output file.
            only_chain (str): If set, only the specified chain will be reported.

        Returns:
            pd.DataFrame: The filtered DataFrame.
        
        """
        if ignore_other_entries:
            # remove entries that are not in the input file
            df = df.loc[df.apply(lambda row: row["Identifier"].startswith(row["PDB ID"] + "_" + row["Ligand ID"]), axis=1)]
        
        if only_chain is not None:
            # remove entries that do not match the chain
            df = df.loc[df.apply(lambda row: row["Identifier"].startswith(row["PDB ID"] + "_" + row["Ligand ID"] + "_" + only_chain), axis=1)]

        return df

    def parse_soup(self, soup):
        """
        Parses the BeautifulSoup object to extract the table data.

        Args:
            soup (BeautifulSoup): The BeautifulSoup object containing the HTML content.
        
        Returns:
            pd.DataFrame: A DataFrame containing the extracted data.
        """

        # Extract the table data
        table = soup.find('table', {'id': 'ligand-validation'})
        assert table is not None
        column_names = []
        for th in table.find_all('th'):
            column_names.append(th.text.strip())
        rows = []
        for tr in table.find_all('tr')[1:]:
            cells = tr.find_all('td')
            row = [cell.text.strip() for cell in cells]
            rows.append(row)
        
        # Create a DataFrame from the extracted data
        df = pd.DataFrame(rows, columns=column_names)
        return df

    def format_df(self, df):
        """
        clean the DataFrame by converting the data types into the correct format.

        Args:
            df (pd.DataFrame): The DataFrame to be cleaned.
        
        Returns:
            pd.DataFrame: The cleaned DataFrame.
        
        """

        # convert percentage columns to float
        percentage_columns = ["Ranking for goodness of fit", 
                              "Ranking for geometry", 
                              "Model completeness"
                              ]
        for col in percentage_columns:
            df[col] = df[col].str.replace('%', '').astype(float) / 100

        # convert float columns to float
        float_columns = ["Real space R factor",
                         "Real space correlation coefficient",
                         "RMSZ-bond-length",
                         "RMSZ-bond-angle",
                         "Average occupancy"
                        ]
        for col in float_columns:
            df[col] = df[col].astype(float)
        
        # convert int columns to int
        int_columns = ["Atomic clashes",
                       "Stereochemical errors",
                       "Outliers of bond length",
                       "Outliers of bond angle"
                      ]
        for col in int_columns:
            df[col] = df[col].str.replace('-', '0').astype(int)
        
        # sort the columns
        df = df[["PDB ID", "Ligand ID", "Identifier"] + percentage_columns + float_columns + int_columns + ["Ligand sdf"]]

        # sort the entries according to the PDB ID and Ligand ID
        df = df.sort_values(by=["PDB ID", "Ligand ID"])
        # reset the index
        df = df.reset_index(drop=True)
        return df
    
    def inspect_batch(self, df, n_workers=10, ignore_other_entries=False, only_chain=None):
        """
        wrapper function for async inspection of a batch of pdb-ligand pairs.

        Args:
            df (pd.DataFrame): A DataFrame containing the columns "PDB ID" and "Ligand ID".
            n_workers (int): The number of worker tasks to create.
            ignore_other_entries (bool): If True, only the entries in the input file will be kept in the output file.
            only_chain (str): If set, only the specified chain will be reported.
        
        Returns:
            pd.DataFrame: A DataFrame containing the inspection results for each pdb-ligand pair.
        
        """

        asyncio.run(self.inspect_batch_async(df, n_workers=n_workers))
        self.df = self.format_df(self.df)
        
        self.df = self.filter_df(self.df, ignore_other_entries=ignore_other_entries, only_chain=only_chain)

        return self.df

        

    async def worker(self, worker_id, q, session, live, task):
        """
        worker function for downloading files asynchronously.

        Args:
            worker_id (int): The ID of the worker.
            q (asyncio.Queue): The queue containing the URLs and file paths.
            session (aiohttp.ClientSession): The aiohttp session.
            live (rich.live.Live): The rich live object for displaying progress.
            task (rich.progress.Task): The task object for progress tracking.
        """
        while True:
            row = await q.get()
            item_id = f"{row['PDB ID']}_{row['Ligand ID']}"
            CURRENT_ITEMS.append(item_id)
            CURRENT_ITEMS.sort()
            state = None

            # update live
            self.update_live(live)

            # download the table with ligand metrics from RCSB
            async with session.get(row['url']) as response:
                if response.status == 200:
                    content = await response.text()
                    state = "success"
                else:
                    content = None
                    logger.error(f"Failed to download the ligand metrics from RCSB for {item_id}. Status code: {response.status}. Please check if the system is available under the following link: {row['url']}")
                    state = "failed"
                    q.task_done()
                    CURRENT_ITEMS.remove(item_id)
                    self.update_live(live)
                    PROGRESS.update(task, advance=1)
                    continue
                             

            soup = BeautifulSoup(content, 'html.parser')
            try:
                df = self.parse_soup(soup)
                df["PDB ID"] = row['PDB ID']
                df["Ligand ID"] = row['Ligand ID']
            except AssertionError:
                logger.error(f"Failed to obtain the table with ligand metrics from RCSB for {item_id}. Please check if the ligand validation report is available for this PDB ID and ligand ID.")
                state = "failed"
                q.task_done()
                CURRENT_ITEMS.remove(item_id)
                self.update_live(live)
                PROGRESS.update(task, advance=1)
                continue

            # download the ligand sdf file
            try:
                ligand_sdf = await self.ligand_downloader.download_ligand_async(row["Ligand ID"], session)
                df.loc[:, "Ligand sdf"] = ligand_sdf
            except ValueError:
                state = "failed"
                ligand_sdf = None

            
            # wrapping up
            self.df = pd.concat([self.df, df], ignore_index=True)

            # if state == "success":
            #     logger.info(f"{item_id} webscraping successful {u'\N{check mark}'}")
                
            q.task_done()
            CURRENT_ITEMS.remove(item_id)
            self.update_live(live)
            PROGRESS.update(task, advance=1)

    def update_live(self, live):
        """
        Updates the live object with the current state of the items being processed.

        Args:
            live (rich.live.Live): The rich live object for displaying progress.
        """
        grid = Table.grid()
        n_columns = 5
        for i in range(n_columns):
            grid.add_column()

        for b in batched(CURRENT_ITEMS, n_columns):
            spinners = [Spinner("dots", text=f"Processing {i}...     ") for i in b]
            grid.add_row(*spinners)

        group = Group(
            Panel(grid, expand=False),
            PROGRESS
        )
        live.update(group)
        
    async def inspect_batch_async(self, df, n_workers=10):
        """
        Inspects a batch of pdb-ligand pairs using asyncio and aiohttp.

        Args:
            df (pd.DataFrame): A DataFrame containing the columns "PDB ID" and "Ligand ID".
            n_workers (int): The number of worker tasks to create.
        
        Returns:
            pd.DataFrame: A DataFrame containing the inspection results for each pdb-ligand pair.
        
        """

        if "PDB ID" not in df.columns or "Ligand ID" not in df.columns:
            raise ValueError("DataFrame must contain 'PDB ID' and 'Ligand ID' columns.")
        
        if df.empty:
            raise ValueError("DataFrame is empty.")
        
        # build the urls
        df["url"] = df.apply(lambda row: self.build_url(row["PDB ID"], row["Ligand ID"]), axis=1)

        q = asyncio.Queue()

        for _, row in df.iterrows():
            await q.put(row)

        async with aiohttp.ClientSession() as session:
            workers = []
            n_workers = min(n_workers, len(df))
            with Live(refresh_per_second=10) as live:
                task1 = PROGRESS.add_task("Web scraping ... ", total=len(df))
                for i in range(n_workers):
                    task = asyncio.create_task(self.worker(i, q, session, live, task1))
                    workers.append(task)

                await q.join()

        with suppress(asyncio.CancelledError):
            for t in workers:
                t.cancel()

            await asyncio.gather(*workers, return_exceptions=True)

    def get_bond_count(self, mol):
        """
        Get the number of bonds between heavy atoms in the molecule.
        
        Parameters:
            mol (object): The molecule object.
            
        Returns:
            int: Number of bonds in the molecule.
        """
        if mol is None:
            return None
        
        bond_count = 0
        bonds = mol.GetBonds()
        for bond in bonds:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetSymbol() != "H" and atom2.GetSymbol() != "H":
                bond_count += 1

        return bond_count
    
    def visualize_bond_count(self, mol, filepath):
        """
        Visualize the molecule and annotate the bonds.
        """
        # remove all hydrogens
        mol = Chem.rdmolops.RemoveHs(mol)

        # Generate 2D coordinates for the molecule
        AllChem.Compute2DCoords(mol)

        
        # Draw the molecule with bond annotations
        options = MolDrawOptions()
        options.addBondIndices = True
        drawer = MolDraw2DSVG(600, 600)
        drawer.SetDrawOptions(options)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        with open(filepath, "w") as f:
            f.write(svg)

    def get_angle_count(self, mol):
        """
        Get the number of angles in the molecule.
        
        Parameters:
            mol (object): The molecule object.
            
        Returns:
            int: Number of angles between atoms in the molecule.
        """
        if mol is None:
            return None

        # an angle is defined by 3 atoms
        angle_thruples = []
        atoms = mol.GetAtoms()

        # iterate through all atoms and find the angles
        for atom in atoms:
            if atom.GetSymbol() == "H":
                continue
            atom_idx = atom.GetIdx()
            
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == "H":
                    continue
                neighbor_idx = neighbor.GetIdx()

                second_neighbors = neighbor.GetNeighbors()
                for second_neighbor in second_neighbors:
                    if second_neighbor.GetSymbol() == "H":
                        continue
                    second_neighbor_idx = second_neighbor.GetIdx()
                    if second_neighbor_idx == atom_idx:
                        continue
                    
                    new_thruple = tuple(sorted([atom_idx, neighbor_idx, second_neighbor_idx]))
                    if new_thruple not in angle_thruples:
                        angle_thruples.append(new_thruple)
        return len(angle_thruples)

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

        # check Ligand occurence in PDB
        df.loc[:,"primary entry"] = df.apply(lambda row: row["Identifier"].startswith(row["PDB ID"]), axis=1)
        df["secondary entry"] = df.apply(lambda row: not row["primary entry"], axis=1)

        df.loc[:, "Ligand occurence in PDB"] = df.apply(lambda row: df.groupby("PDB ID").get_group(row["PDB ID"])["secondary entry"].sum() +1, axis=1)
        if ignore_other_entries:
            df = df[df["primary entry"]]

        df.loc[:, "Chain ID"] = df.apply(lambda row: row["Identifier"].split("_")[2], axis=1)
        if only_chain is not None:
            df = df[df["Chain ID"] == only_chain]

        df = df.drop(columns=["primary entry", "secondary entry"])
        df.loc[:,"Structural outliers"] = df.loc[:,["Outliers of bond length", "Outliers of bond angle"]].sum(axis=1)

        df.loc[:, "Number of bonds"] = df.apply(lambda row: self.get_bond_count(self.load_mol_from_sdf_block(row["Ligand sdf"])), axis=1)

        df.loc[:, "Number of angles"] = df.apply(lambda row: self.get_angle_count(self.load_mol_from_sdf_block(row["Ligand sdf"])), axis=1)

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
    
    def load_mol_from_sdf_block(self, sdf_block):
        """
        Loads a molecule from an SDF block.

        Args:
            sdf_block (str): The SDF block containing the molecule data.
        
        Returns:
            rdkit.Chem.Mol: The loaded molecule.
        """
        if sdf_block is None:
            return None
        
        RDLogger.DisableLog('rdApp.*')
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_block)
        mol = next(suppl)
        RDLogger.EnableLog('rdApp.*')
        return mol
    
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
        ax.set_xlim(xmin=0)

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
        ax.set_xlim(xmin=0)

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
        ax.set_xlim(xmin=0)
    
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