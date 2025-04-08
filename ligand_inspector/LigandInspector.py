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
    def __init__(self):
        self.df = pd.DataFrame()

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

        url = self.build_url(pdb_id, ligand_id)
        response = requests.get(url)
        soup = BeautifulSoup(response.content, 'html.parser')
        try:
            df = self.parse_soup(soup)
        except AssertionError:
            raise ValueError(f"Could not find the table for PDB ID: {pdb_id} and Ligand ID: {ligand_id}, url: {url}")
        df["PDB ID"] = pdb_id
        df["Ligand ID"] = ligand_id

        df = self.format_df(df)

        df = self.filter_df(df, ignore_other_entries=ignore_other_entries, only_chain=only_chain)

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
        df = df[["PDB ID", "Ligand ID","Identifier"] + percentage_columns + float_columns + int_columns]

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
                self.df = pd.concat([self.df, df], ignore_index=True)
            except AssertionError:
                state = "failed"

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

        
            
