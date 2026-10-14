"""
LigandDownloader
================

This class is responsible for downloading ligand data from the RCSB PDB database.
It downloads the idealized ligand in sdf format using asynchronous http requests.
For that it uses the api of rcsb.org, which downloads the ideal ligand sdf file from the Chemical Component Dictionary (CCD) of the RCSB PDB database.

The API works as follows (on the example of ligands EOF and 4WI):

.. code-block:: text
    https://download.rcsb.org/batch/ccd/EOF_ideal.sdf:4WI_ideal.sdf

this downloads a zip file containing the ideal sdf files of the ligands EOF_ideal.sdf and 4WI_ideal.sdf.
The zip file is then
"""

import os
import asyncio
import aiohttp
from .log import logger
import tempfile
from itertools import batched
import zipfile


class LigandDownloader:
    """
    A class to download ligand data from the RCSB PDB database.
    """

    def __init__(self):
        pass

    def build_urls(self, ligand_ids):
        """
        Build the URL for downloading the ligand data.

        Args:
            ligand_ids (list): A list of ligand IDs to download.

        Returns:
            list: The URLs for downloading the ligand data in batches of 20 batches each.
        """
        urls = [
            "https://download.rcsb.org/batch/ccd/"
            + ":".join(f"{ligand_id}_ideal.sdf" for ligand_id in batch)
            for batch in batched(ligand_ids, 20)
        ]
        return urls

    async def download_ligands_async(self, ligand_ids, session):
        """
        Download the ligand data asynchronously.

        Args:
            ligand_ids (str or list): The ID or list of IDs of the ligand(s) to download.
            session (aiohttp.ClientSession): The aiohttp session to use for the request.

        Returns:
            str: The ideal ligand sdf file content as a string.
        """
        if isinstance(ligand_ids, str):
            ligand_ids = [ligand_ids]

        urls = self.build_urls(ligand_ids)

        for url in urls:
            logger.info(f"Downloading ligand(s) {ligand_ids} from {url}")

            async with session.get(url) as response:
                if response.status == 200:
                    # create temporary directory to store the zip file
                    tmp_dir = tempfile.mkdtemp()
                    tmp_file = os.path.join(tmp_dir, "ligands.zip")
                    with open(tmp_file, "wb") as f:
                        f.write(await response.read())

                    # extract the zip file and read the sdf files
                    with zipfile.ZipFile(tmp_file, "r") as zip_ref:
                        zip_ref.extractall(tmp_dir)

                    # read the sdf files
                    sdf_contents = {}
                    for ligand_id in ligand_ids:
                        print(list(os.listdir(tmp_dir)))
                        sdf_file = os.path.join(
                            tmp_dir, f"{ligand_id.lower()}_ideal.sdf"
                        )
                        if os.path.exists(sdf_file):
                            with open(sdf_file, "r") as f:
                                sdf_contents[ligand_id] = f.read()
                        else:
                            logger.error(
                                f"SDF file for ligand {ligand_id} not found in the zip archive."
                            )
                            sdf_contents[ligand_id] = None

                    return sdf_contents
                else:
                    logger.error(
                        f"Failed to download ligand {ligand_id}. Status code: {response.status}. Please check if the sdf file is available under the following link: {url}"
                    )
                    # raise ValueError(f"Failed to download ligand {ligand_id}. Status code: {response.status}")
                    return None

    async def download_single_ligand_async(self, ligand_id):
        """
        Download a single ligand asynchronously.

        Args:
            ligand_id (str): The ID of the ligand to download.

        Returns:
            str: The ideal ligand sdf file content as a string.
        """
        async with aiohttp.ClientSession() as session:
            ligand_dict = await self.download_ligands_async(ligand_id, session)
            return ligand_dict.get(ligand_id, None) if ligand_dict else None

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
            tasks = [
                self.download_ligands_async(ligand_id, session)
                for ligand_id in ligand_ids
            ]
            results = await asyncio.gather(*tasks)
            # extract the sdf content from the returned dictionaries
            return [
                res.get(ligand_id, None) for res, ligand_id in zip(results, ligand_ids)
            ]

    def download_multiple_ligands(self, ligand_ids):
        """
        Download multiple ligands synchronously.

        Args:
            ligand_ids (list): A list of ligand IDs to download.

        Returns:
            list: A list of ideal ligand sdf file content as strings in the same order as ligand_ids.
        """
        return asyncio.run(self.download_multiple_ligands_async(ligand_ids))
