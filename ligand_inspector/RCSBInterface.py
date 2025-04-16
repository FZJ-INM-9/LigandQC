"""
RCSBInterface
=============

This module provides an interface to the RCSB PDB database for retrieving ligand information.
It handles all API calls and data processing.

"""

from rcsbapi.data import DataQuery
from rcsbapi.search import AttributeQuery
import pandas as pd
from .Logger import logger
from itertools import batched
from tqdm.auto import tqdm

class RCSBInterface:
    """
    RCSBInterface is a class that interacts with the RCSB PDB database to retrieve ligand information.
    """

    def __init__(self):
        pass

    def get_all_pdbids(self, ligand_ids):
        """
        Retrieves all PDB IDs associated with a given ligand ID.

        Args:
            ligand_id (list[str]): A list of ligand IDs to search for.

        Returns:
            list: A list of PDB IDs associated with the ligand ID.
        """
        logger.info(f"Fetching all PDB IDs for provided ligand_ids...")
        query = AttributeQuery(
            attribute="rcsb_chem_comp_container_identifiers.comp_id",
            operator="in",
            value=ligand_ids
        )
        results = list(query())
        logger.info("Done fetching PDB IDs.")
        return results
    
    
    def get_ligand_quality_data(self, ligand_ids, pdb_ids=None):
        """
        Retrieves the ligand quality data for a given list of ligand IDs and PDB IDs.
        If PDB IDs are provided, only the data for those PDB IDs will be retrieved.
        If no PDB IDs are provided, all data for the ligand IDs will be retrieved.

        Args:
            ligand_ids (list): A list of ligand IDs to search for.
            pdb_ids (list, optional): A list of PDB IDs to filter the results by. Defaults to None.

        Returns:
            pd.DataFrame: A DataFrame containing the ligand quality data.
        
        """
        if pdb_ids is None:
            logger.info("Fetching all PDB IDs for the given ligand IDs.")
            pdb_ids = []
            pdb_ids.extend(self.get_all_pdbids(ligand_ids))
        
        logger.info(f"Fetching ligand quality data for ligand IDs: {ligand_ids} and PDB IDs: {pdb_ids}")
        data = []
        batch_size = 500
        n_batches = (len(pdb_ids) + batch_size - 1) // batch_size
        logger.info(f"Processing {n_batches} batches of PDB IDs.")
        for batch in tqdm(batched(pdb_ids, batch_size), total=n_batches):
            query = DataQuery(
                input_type="entries",
                input_ids=batch,
                return_data_list=[
                    "nonpolymer_entities.pdbx_entity_nonpoly.comp_id",
                    "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_entity_instance_container_identifiers",
                    "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score",
                    ])
            query.exec()
            data += query.get_response()['data']['entries']
        
        logger.info("Processing the response data...")
        processed_data = []
        for i in data:
            for j in i["nonpolymer_entities"]:
                for k in j["nonpolymer_entity_instances"]:
                    this_data = {}
                    this_data["PDB ID"] = i["rcsb_id"]
                    this_data["Ligand ID"] = k["rcsb_nonpolymer_entity_instance_container_identifiers"]["comp_id"]
                    this_data["Chain ID"] = k["rcsb_nonpolymer_entity_instance_container_identifiers"]["auth_asym_id"]
                    this_data["Residue ID"] = k["rcsb_nonpolymer_entity_instance_container_identifiers"]["auth_seq_id"]
                    this_data["Ranking for goodness of fit"] = k["rcsb_nonpolymer_instance_validation_score"][0]["ranking_model_fit"]
                    this_data["Ranking for geometry"] = k["rcsb_nonpolymer_instance_validation_score"][0]["ranking_model_geometry"]
                    this_data["Model completeness"] = k["rcsb_nonpolymer_instance_validation_score"][0]["completeness"]
                    this_data["Real space R factor"] = k["rcsb_nonpolymer_instance_validation_score"][0]["RSR"]
                    this_data["Real space correlation coefficient"] = k["rcsb_nonpolymer_instance_validation_score"][0]["RSCC"]
                    this_data["RMSZ-bond-length"] = k["rcsb_nonpolymer_instance_validation_score"][0]["mogul_bonds_RMSZ"]
                    this_data["RMSZ-bond-angle"] = k["rcsb_nonpolymer_instance_validation_score"][0]["mogul_angles_RMSZ"]
                    this_data["Average occupancy"] = k["rcsb_nonpolymer_instance_validation_score"][0]["average_occupancy"]
                    this_data["Atomic clashes"] = k["rcsb_nonpolymer_instance_validation_score"][0]["intermolecular_clashes"]
                    this_data["Stereochemical errors"] = k["rcsb_nonpolymer_instance_validation_score"][0]["stereo_outliers"]
                    this_data["Outliers of bond length"] = k["rcsb_nonpolymer_instance_validation_score"][0]["mogul_bond_outliers"]
                    this_data["Outliers of bond angle"] = k["rcsb_nonpolymer_instance_validation_score"][0]["mogul_angle_outliers"]
                    this_data["Ligand occurence in PDB"] = len(pdb_ids)
                    this_data["Structural outliers"] = this_data["Outliers of bond length"] + this_data["Outliers of bond angle"]
                    this_data["Number of bonds"] = k["rcsb_nonpolymer_instance_validation_score"][0]["num_mogul_bonds_RMSZ"]
                    this_data["Number of angles"] = k["rcsb_nonpolymer_instance_validation_score"][0]["num_mogul_angles_RMSZ"]
                    this_data["Number of atoms"] = k["rcsb_nonpolymer_instance_validation_score"][0]["natoms_eds"]
                    processed_data.append(this_data)
        df = pd.DataFrame(processed_data)
        df = df.drop_duplicates()
        df = df[df["Ligand ID"].isin(ligand_ids)]
        df = df[df["PDB ID"].isin(pdb_ids)]
        df = df.reset_index(drop=True)
        return df