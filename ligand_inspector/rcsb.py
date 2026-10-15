"""
rcsb
====

This utility module provides an interface to the RCSB PDB database for retrieving ligand information.
It handles all API calls and data processing.

"""

from rcsbapi.data import DataQuery
from rcsbapi.search import AttributeQuery
import pandas as pd
from .log import logger
from itertools import batched
from tqdm.auto import tqdm
from ligand_inspector.const import *


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
            value=ligand_ids,
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
            list: A list of dictionaries containing the ligand quality data.

        """
        if pdb_ids is None:
            logger.info("Fetching all PDB IDs for the given ligand IDs.")
            pdb_ids = []
            pdb_ids.extend(self.get_all_pdbids(ligand_ids))

        logger.info(
            f"Fetching ligand quality data for {len(ligand_ids)} ligand IDs and {len(pdb_ids)} PDB IDs..."
        )
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
                    "exptl.method",
                ],
            )
            query.exec()
            data += query.get_response()["data"]["entries"]

        return data

    def process_ligand_quality_data(self, data, ligand_ids, pdb_ids=None):
        """
        Processes the ligand quality data retrieved from the RCSB PDB database.

        Args:
            data (list): A list of dictionaries containing the ligand quality data.
            ligand_ids (list): A list of ligand IDs to filter the results by.
            pdb_ids (list): A list of PDB IDs to filter the results by.

        Returns:
            pd.DataFrame: A DataFrame containing the processed ligand quality data.
        """
        logger.info("Wrangling the API response data...")
        processed_data = []
        key_map = {
            "Ranking for goodness of fit": "ranking_model_fit",
            "Ranking for geometry": "ranking_model_geometry",
            "Model completeness": "completeness",
            METRIC_RSR: "RSR",
            METRIC_RSCC: "RSCC",
            METRIC_RMSZ_BOND_LENGTH: "mogul_bonds_RMSZ",
            METRIC_RMSZ_BOND_ANGLE: "mogul_angles_RMSZ",
            METRIC_AVG_OCCUPANCY: "average_occupancy",
            METRIC_ATOMIC_CLASHES: "intermolecular_clashes",
            METRIC_STEREO_CHEMICAL_ERRORS: "stereo_outliers",
            "Outliers of bond length": "mogul_bond_outliers",
            "Outliers of bond angle": "mogul_angle_outliers",
            "Number of bonds": "num_mogul_bonds_RMSZ",
            "Number of angles": "num_mogul_angles_RMSZ",
            "Number of atoms": "natoms_eds",
            "Is Subject of Investigation": "is_subject_of_investigation",
        }

        for i in data:
            if i["nonpolymer_entities"] is None:
                logger.warning(
                    f"No nonpolymer entities found for PDB ID: {i['rcsb_id']}. Skipping this entry."
                )
                continue
            for j in i["nonpolymer_entities"]:

                if j["nonpolymer_entity_instances"] is None:
                    continue

                for k in j["nonpolymer_entity_instances"]:
                    this_data = {}
                    this_data["PDB ID"] = i["rcsb_id"]
                    this_data["Ligand occurence in PDB"] = "not evaluated."
                    this_data["Ligand ID"] = k[
                        "rcsb_nonpolymer_entity_instance_container_identifiers"
                    ]["comp_id"]
                    this_data["Chain ID"] = k[
                        "rcsb_nonpolymer_entity_instance_container_identifiers"
                    ]["auth_asym_id"]
                    this_data["Residue ID"] = k[
                        "rcsb_nonpolymer_entity_instance_container_identifiers"
                    ]["auth_seq_id"]
                    if i["exptl"] is None:
                        this_data["Experimental Method"] = None
                    else:
                        this_data["Experimental Method"] = i["exptl"][0]["method"]

                    if k["rcsb_nonpolymer_instance_validation_score"] is None:
                        for key, value in key_map.items():
                            this_data[key] = None
                    else:

                        for key, value in key_map.items():
                            if (
                                value
                                in k["rcsb_nonpolymer_instance_validation_score"][0]
                            ):
                                this_data[key] = k[
                                    "rcsb_nonpolymer_instance_validation_score"
                                ][0][value]
                            else:
                                this_data[key] = None
                    this_data["Is Subject of Investigation"] = (
                        this_data["Is Subject of Investigation"] == "Y"
                    )

                    if (this_data["Outliers of bond length"] is not None) and (
                        this_data["Outliers of bond angle"] is not None
                    ):
                        this_data[METRIC_STRUCTURAL_OUTLIERS] = (
                            this_data["Outliers of bond length"]
                            + this_data["Outliers of bond angle"]
                        )
                    else:
                        this_data[METRIC_STRUCTURAL_OUTLIERS] = None
                    processed_data.append(this_data)
        df = pd.DataFrame(processed_data)
        df = df.drop_duplicates()
        df = df[df["Ligand ID"].isin(ligand_ids)]
        if pdb_ids is not None:
            df = df[df["PDB ID"].isin(pdb_ids)]
        df = df.reset_index(drop=True)
        return df
