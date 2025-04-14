"""
LigandInspector
===============

This class is responsible for inspecting the quality of the specified ligand and PDB id.

"""

import os
import aiohttp
import pandas as pd
import asyncio
import requests
from bs4 import BeautifulSoup
from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, MofNCompleteColumn, TimeRemainingColumn, track
from rich.live import Live
from rich.spinner import Spinner
from rich.panel import Panel
from rich.table import Table
from contextlib import suppress
import matplotlib.pyplot as plt
import seaborn as sns
from .RMSZInterpreter import RMSZInterpreter
from .LigandDownloader import LigandDownloader
from rdkit import Chem


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

        # url = self.build_url(pdb_id, ligand_id)
        # response = requests.get(url)
        # soup = BeautifulSoup(response.content, 'html.parser')
        # try:
        #     df = self.parse_soup(soup)
        # except AssertionError:
        #     raise ValueError(f"Could not find the table for PDB ID: {pdb_id} and Ligand ID: {ligand_id}, url: {url}")
        # df["PDB ID"] = pdb_id
        # df["Ligand ID"] = ligand_id

        # df = self.format_df(df)
        

        # df = self.filter_df(df, ignore_other_entries=ignore_other_entries, only_chain=only_chain)

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
                    state = "failed"

            soup = BeautifulSoup(content, 'html.parser')
            try:
                df = self.parse_soup(soup)
                df["PDB ID"] = row['PDB ID']
                df["Ligand ID"] = row['Ligand ID']
            except AssertionError:
                state = "failed"

            # download the ligand sdf file
            ligand_sdf = await self.ligand_downloader.download_ligand_async(row["Ligand ID"], session)
            df.loc[:, "Ligand sdf"] = ligand_sdf
            

            # wrapping up
            self.df = pd.concat([self.df, df], ignore_index=True)
            if state == "success":
                live.console.log(item_id + u'\N{check mark}')
            else:
                live.console.log(item_id + u'\N{cross mark}')

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
        table = Table.grid()
        for i in CURRENT_ITEMS:
            table.add_row(Spinner("dots", text=f"Processing {i} ..."))
        # print(f"Current items: {CURRENT_ITEMS}")
        table.add_row(PROGRESS)
        live.update(Panel(table))
        
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
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_block)
        mol = next(suppl)
        return mol
    
    def evaluate_rmsz_bond_length(self, row):
        """
        evaluates the RMSZ-bond-length of the ligand
        """
        mol = self.load_mol_from_sdf_block(row["Ligand sdf"])
        return self.rmsz_interpreter.interpret_molecule(mol, row["RMSZ-bond-length"], mode="bond")
    
    def evaluate_rmsz_bond_angle(self, row):
        """
        evaluates the RMSZ-bond-angle of the ligand
        """
        mol = self.load_mol_from_sdf_block(row["Ligand sdf"])
        return self.rmsz_interpreter.interpret_molecule(mol, row["RMSZ-bond-angle"], mode="angle")
    
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

        fig, axs = plt.subplots(4, 3, figsize=(15, 15))
        n_bins = len(df) // 10
        if n_bins < 1:
            n_bins = 1

        df = self.evaluate(df)


        colors = {"good": "green", "ok": "yellow", "bad": "orangered"}
        order = {"good": 0, "ok": 1, "bad": 2}

        sns.histplot(df, x="Real space R factor", bins=n_bins,  ax=axs[0,0], hue="Real space R factor evaluation", multiple="stack", palette=colors)
        axs[0,0].set_xlim([0, 1])

        sns.histplot(df, x="Real space correlation coefficient", bins=n_bins,  ax=axs[0,1], hue="Real space correlation coefficient evaluation", multiple="stack", palette=colors)
        axs[0,1].set_xlim([0, 1])

        # sns.scatterplot(data=df, x="Ranking for goodness of fit", y="Ranking for geometry", ax=axs[2,0], hue="Real space correlation coefficient evaluation")
        # axs[2,0].set_xlim([0, 1])
        # axs[2,0].set_ylim([0, 1])

        sns.histplot(df, x="RMSZ-bond-length", bins=n_bins,  ax=axs[0,2], multiple="stack")

        sns.histplot(df, x="RMSZ-bond-angle", bins=n_bins,  ax=axs[1,0],multiple="stack")

        sns.histplot(df, x="Structural outliers", bins=n_bins,  ax=axs[1,1], hue="Structural outliers evaluation", multiple="stack", palette=colors)

        # pie plot for atomic clashes
        pie_data = df.groupby("Atomic clashes evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (0)", "ok": "ok (1-2)", "bad": "bad (>2)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        axs[1,2].pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
        axs[1,2].set_xlabel("Atomic clashes")

        # pie plot for the steriochemical errors
        pie_data = df.groupby("Stereochemical errors evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (0)", "bad": "bad (>0)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        axs[2,0].pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
        axs[2,0].set_xlabel("Stereochemical errors")

        # pie plot for the average occupancy
        pie_data = df.groupby("Average occupancy evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        print(pie_data_df)
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))
        labels_map = {"good": "good (>0.95)", "ok": "ok (0.5-0.95)", "bad": "bad (<0.5)"}
        pie_data_df = pie_data_df.rename(index=labels_map)

        axs[2,1].pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
        axs[2,1].set_xlabel("Average occupancy evaluation")
        # sns.histplot(df, x="Average occupancy", bins=n_bins,  ax=axs[2,1], hue="Real space correlation coefficient evaluation", multiple="stack")

        # pie plot for the overall evaluation
        pie_data = df.groupby("overall evaluation").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.loc[:, "color"] = pie_data_df.index.map(lambda x: colors[x.split(" ")[0]])
        pie_data_df.sort_index(inplace=True, key=lambda x: x.map(order))

        axs[2,2].pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, colors=pie_data_df["color"], counterclock=False)
        axs[2,2].set_xlabel("Overall evaluation")
        

        # plot the histograms
        sns.histplot(df, x="Ranking for goodness of fit", bins=n_bins,  ax=axs[3,0])
        axs[3,0].set_xlim([0, 1])

        sns.histplot(df, x="Ranking for geometry", bins=n_bins,  ax=axs[3,1])
        axs[3,1].set_xlim([0, 1])

        # plot the occurence of the ligand in the PDB as a pie chart
        pie_data = df.groupby("Ligand occurence in PDB").size()
        pie_data_df = pd.DataFrame(pie_data)
        pie_data_df = pie_data_df.rename(columns={0: "count"})
        pie_data_df.sort_index(inplace=True)
        axs[3,2].pie(pie_data_df["count"], labels=pie_data_df.index, autopct='%1.1f%%', startangle=90, counterclock=False)
        axs[3,2].set_xlabel("Ligand occurence in PDB")


        plt.tight_layout()
        plt.savefig(filename)