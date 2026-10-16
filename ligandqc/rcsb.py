"""
rcsb
====

This utility module provides an interface to the RCSB PDB database for retrieving ligand information.
It handles all API calls and data processing.

"""

import click
from rcsbapi.data import DataQuery, ALL_STRUCTURES
from rcsbapi.search import AttributeQuery
import pandas as pd
from ligandqc.log import logger
from itertools import batched
from tqdm.auto import tqdm
from ligandqc.const import *
import re
from collections import defaultdict


class RCSBInterface:
    """
    RCSBInterface is a class that interacts with the RCSB PDB database to retrieve ligand information.
    """

    def __init__(self):
        pass

    def get_pdb_ids(self, ligand_ids):
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
            pdb_ids.extend(self.get_pdb_ids(ligand_ids))

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

    def get_all_ids(self, ligand_id_output_file, pdb_id_output_file, mapping_output_file):
        """
        Creates a dataset of ligands from the PDB with a mapping of PDB IDs to ligand IDs.

        Args:
            ligand_id_output_file (str): path to the output file for ligand IDs
            pdb_id_output_file (str): path to the output file for PDB IDs
            mapping_output_file (str): path to the output file for mapping of PDB IDs to ligand IDs
        
        Returns:
            dict: a dictionary containing the ligand IDs as keys and the PDB IDs in which they occur as values
        """
        # get all ligand IDs present in any structure
        logger.info("Fetching all PDB IDs, that contain ligands...")
        q = DataQuery(
            input_type="entries",
            input_ids=ALL_STRUCTURES,
            return_data_list=["nonpolymer_entities.pdbx_entity_nonpoly.comp_id",
                              "exptl.method",
                              "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.is_subject_of_investigation"]
        )
        q.exec(progress_bar=True)

        entries_with_ligand = []
        for i in q.get_response()["data"]["entries"]:
            if i['exptl'] is not None:
                experimental_method = i['exptl'][0]['method']
                if experimental_method != "X-RAY DIFFRACTION":
                    continue
            if i['nonpolymer_entities'] is not None:
                entries_with_ligand.append(i)

        logger.info(f"Found {len(entries_with_ligand)} entries with ligands in the Protein Data Bank.")

        ligand_id_pattern = re.compile(r"^[0-9A-Z]{3}$")
        invalid_ligand_ids = set()
        mapping = defaultdict(lambda: [])
        for i in entries_with_ligand:
            for j in i["nonpolymer_entities"]:
                lig_id = j["pdbx_entity_nonpoly"]["comp_id"]
                if not ligand_id_pattern.match(lig_id):
                    invalid_ligand_ids.add(lig_id)
                    continue
                
                nonpolymer_entity_instances = j['nonpolymer_entity_instances']

                # check if any of the nonpolymer_entity_instances are subject of investigation
                if nonpolymer_entity_instances is None:
                    continue

                for k in nonpolymer_entity_instances:
                    if k['rcsb_nonpolymer_instance_validation_score'] is None:
                        continue
                    for l in k['rcsb_nonpolymer_instance_validation_score']:
                        if l.get('is_subject_of_investigation', "N") == "Y":
                            # add the PDB ID to the mapping if it is not already present
                            pdb_id = i['rcsb_id']
                            mapping[lig_id].append(pdb_id)
        
        logger.info(f"Found {len(mapping)} valid ligand IDs.")

        # export to files
        ligand_ids = list(mapping.keys())
        ligand_ids.sort()
        with open(ligand_id_output_file, 'w') as f:
            f.write("\n".join(ligand_ids))
        logger.info(f"Exported {len(ligand_ids)} ligand IDs to {ligand_id_output_file}.")

        # export PDB IDs
        pdb_ids = set()
        for i in mapping.values():
            pdb_ids.update(i)
        pdb_ids = list(pdb_ids)
        pdb_ids.sort()
        with open(pdb_id_output_file, 'w') as f:
            f.write("\n".join(pdb_ids))
        logger.info(f"Exported {len(pdb_ids)} PDB IDs to {pdb_id_output_file}.")

        # export mapping
        mapping_condensed = defaultdict(lambda: [])
        for k, v in mapping.items():
            counter = defaultdict(lambda: 0)
            for i in v:
                counter[i] += 1
            for i in counter.keys():
                mapping_condensed[k].append(f"{i} ({counter[i]})")
        data = [{"Ligand ID": k, "PDB IDs": ";".join(sorted(v))} for k, v in mapping_condensed.items()]
        df = pd.DataFrame(data)
        df.sort_values(by="Ligand ID", inplace=True)
        df.reset_index(drop=True, inplace=True)
        df.to_csv(mapping_output_file, index=False)
        logger.info(f"Exported mapping to {mapping_output_file}.")

        if invalid_ligand_ids:
            logger.warning(f"{len(invalid_ligand_ids)} invalid ligand IDs found: {', '.join(sorted(invalid_ligand_ids))}")

        return mapping
                


        

        
@click.command()
@click.option('--ligand-id-output-file',
              "-ol",
              type=click.Path(exists=False),
              help="Path to the output file for ligand IDs.",
              default="ligand_ids.txt",
              show_default=True)
@click.option('--pdb-id-output-file',
              "-op",
              type=click.Path(exists=False),
              help="Path to the output file for PDB IDs.",
              default="pdb_ids.txt",
              show_default=True)
@click.option('--mapping-output-file',
              "-om",
              type=click.Path(exists=False),
              help="Path to the output file for the mapping.",
              default="mapping.csv",
              show_default=True)
def get_all_ids(ligand_id_output_file, pdb_id_output_file, mapping_output_file):
    """
    Scrapes the entire PDB for ligand ids and creates a mapping of PDB IDs to ligand IDs.

    Saves three files:
    1. A list of all ligand IDs found in the PDB as Subjects of Investigation (SOI)
    2. A list of all PDB IDs that contain those ligand ids
    3. A CSV file that maps ligand IDs to the PDB IDs and how many times they occur in each PDB ID.

    The default output files are saved in the current working directory.

    """
    interface = RCSBInterface()
    interface.get_all_ids(ligand_id_output_file, pdb_id_output_file, mapping_output_file)

if __name__ == "__main__":
    get_all_ids()