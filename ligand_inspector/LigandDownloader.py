"""
LigandDownloader
================

This class is responsible for downloading ligand data from the RCSB PDB database.
It downloads the idealized ligand in sdf format using asynchronous http requests.

"""

import os
import asyncio
import aiohttp
from .Logger import logger

class LigandDownloader:
    """
    A class to download ligand data from the RCSB PDB database.
    """

    def __init__(self):
        pass

    def build_url(self, ligand_id):
        """
        Build the URL for downloading the ligand data.
        
        Args:
            ligand_id (str): The ID of the ligand to download.
        
        Returns:
            str: The URL for downloading the ligand data.
        """
        return f"http://ligand-expo.rcsb.org/reports/{ligand_id[0]}/{ligand_id}/{ligand_id}_ideal.sdf"
    
    async def download_ligand_async(self, ligand_id, session):
        """
        Download the ligand data asynchronously.
        
        Args:
            ligand_id (str): The ID of the ligand to download.
            session (aiohttp.ClientSession): The aiohttp session to use for the request.
        
        Returns:
            str: The ideal ligand sdf file content as a string.
        """
        url = self.build_url(ligand_id)
        async with session.get(url) as response:
            if response.status == 200:
                return await response.text()
            else:
                logger.error(f"Failed to download ligand {ligand_id}. Status code: {response.status}. Please check if the sdf file is available under the following link: {url}")
                raise ValueError(f"Failed to download ligand {ligand_id}. Status code: {response.status}")
        
    async def download_single_ligand_async(self, ligand_id):
        """
        Download a single ligand asynchronously.
        
        Args:
            ligand_id (str): The ID of the ligand to download.
        
        Returns:
            str: The ideal ligand sdf file content as a string.
        """
        async with aiohttp.ClientSession() as session:
            return await self.download_ligand_async(ligand_id, session)
            
    def download_single_ligand(self, ligand_id):
        """
        Download the ligand data synchronously.
        
        Args:
            ligand_id (str): The ID of the ligand to download.
        
        Returns:
            str: The ideal ligand sdf file content as a string.
        """
        return asyncio.run(self.download_single_ligand_async(ligand_id))

    async def download_multiple_ligands_async(self, ligand_ids):
        """
        Download multiple ligands asynchronously.
        
        Args:
            ligand_ids (list): A list of ligand IDs to download.
        
        Returns:
            list: A list of ideal ligand sdf file content as strings in the same order as ligand_ids.
        """
        async with aiohttp.ClientSession() as session:
            tasks = [self.download_ligand_async(ligand_id, session) for ligand_id in ligand_ids]
            return await asyncio.gather(*tasks)
    
    def download_multiple_ligands(self, ligand_ids):
        """
        Download multiple ligands synchronously.
        
        Args:
            ligand_ids (list): A list of ligand IDs to download.
        
        Returns:
            list: A list of ideal ligand sdf file content as strings in the same order as ligand_ids.
        """
        return asyncio.run(self.download_multiple_ligands_async(ligand_ids))