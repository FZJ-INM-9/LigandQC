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

class LigandInspector:
    def __init__(self):
        pass

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
    
    def inspect_single(self, pdb_id, ligand_id):
        """
        Inspects a single pdb-ligand pair.

        Args:
            pdb_id (str): The PDB ID of the protein.
            ligand_id (str): The ligand ID.

        Returns:
            pd.DataFrame: A DataFrame containing the inspection results.

        """

        url = self.build_url(pdb_id, ligand_id)
        response = requests.get(url)
        soup = BeautifulSoup(response.content, 'html.parser')
        df = self.parse_soup(soup)
        df["PDB ID"] = pdb_id
        df["Ligand ID"] = ligand_id

        df = self.format_df(df)

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


        return df