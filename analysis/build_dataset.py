"""
Build Dataset
-------------

This script builds a dataset for ligand inpection of the entire PDB.

"""

import os
import pandas as pd
import click
import logging
from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, MofNCompleteColumn, TimeRemainingColumn, track
from rich.logging import RichHandler
import re
from rcsbapi.data import DataQuery, ALL_STRUCTURES
from collections import defaultdict

dataset_logger = logging.getLogger("dataset_logger")
if dataset_logger.hasHandlers():
    dataset_logger.handlers.clear()
dataset_logger.propagate = False
dataset_logger.setLevel(logging.INFO)
handler = RichHandler()
handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
dataset_logger.addHandler(handler)

class DatasetBuilder:
    def __init__(self):
        pass

    def build_dataset2(self, input_file, output_file):
        """
        builds the dataset by reading the input file and writing the output file

        Args:
            input_file (str): path to the input file
            output_file (str): path to the output file
        
        Returns:
            None
        """

        data = []
        ligand_id = None
        pdb_id = None
        release_status = None
        ligand_keyword = "_chem_comp.id"
        pdb_keyword = "_chem_comp.pdbx_model_coordinates_db_code"
        release_status_keyword = "_chem_comp.pdbx_release_status"
        new_entry_keyword = "data_"
        max_entries = (10 + 26)**3
        pdb_id_pattern = re.compile(r"^\d[0-9A-Z]{3}$")
        ligand_id_pattern = re.compile(r"^[0-9A-Z]{3}$")

        progress = Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
        )

        with open(input_file, 'r') as f:
            with progress as p:
                current_progress = 0
                task = progress.add_task("Processing file...", total=max_entries)
                
                for line in f:
                    if line.startswith(new_entry_keyword):
                        if ligand_id and not ligand_id_pattern.match(ligand_id):
                            dataset_logger.error(f"Invalid ligand ID found: '{ligand_id}'. Skipping entry.")
                            ligand_id = None
                            pdb_id = None
                        if pdb_id and not pdb_id_pattern.match(pdb_id):
                            dataset_logger.error(f"Invalid PDB ID found for ligand ID '{ligand_id}': '{pdb_id}'. Skipping entry.")
                            pdb_id = None
                            ligand_id = None
                        if ligand_id and release_status != "REL":
                            dataset_logger.warning(f"Obsolete entry found: '{ligand_id}'. Skipping entry.")
                            ligand_id = None
                            pdb_id = None
                        if ligand_id and pdb_id:
                            data.append([pdb_id, ligand_id])
                            advancement = self.determine_progress(ligand_id, current_progress)
                            current_progress += advancement
                            self.update_progress(p, task, advancement)
                        ligand_id = None
                        pdb_id = None
                        continue
                    if line.startswith(ligand_keyword):
                        ligand_id = line.split()[1].strip()
                        continue
                    if line.startswith(pdb_keyword):
                        pdb_id = line.split()[1].strip()
                        pdb_id = pdb_id.strip('"')
                        pdb_id = pdb_id.split(",")[0].upper()
                        continue
                    if line.startswith(release_status_keyword):
                        release_status = line.split()[1].strip()
                    
                
                # Final entry
                if ligand_id and pdb_id:
                    data.append([pdb_id, ligand_id])
                    advancement = max_entries - current_progress
                    current_progress += advancement
                    self.update_progress(p, task, advancement)
        df = pd.DataFrame(data, columns=['PDB ID', 'Ligand ID'])
        df.drop_duplicates(inplace=True)
        df.reset_index(drop=True, inplace=True)
        df.to_csv(output_file, index=False)
        dataset_logger.info(f"Dataset built with {len(df)} entries.")
        dataset_logger.info(f"Dataset saved to {output_file}.")
    
    def determine_progress(self, ligand_id, current_progress):
        """
        Determines the amount of progress to add by checking the ligand_id

        Args:
            ligand_id (str): the current ligand_id
            current_progress (int): the current progress

        Returns:
            int: the amount of progress to add
        
        """
        ligand_id = ligand_id.upper()
        options = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        n_options = len(options)
        first_position = options.index(ligand_id[0])
        second_position = options.index(ligand_id[1])
        third_position = options.index(ligand_id[2])
        current_position = (first_position * n_options**2) + (second_position * n_options) + third_position
        advancement = current_position - current_progress
        return advancement

    def update_progress(self, progress, task, advancement):
        """
        Updates the progress bar with the current ligand_id

        Args:
            progress (Progress): the progress bar
            task (Task): the task to update
            advancement (int): the amount of progress to add
        
        Returns:
            None
        """
        progress.update(task, advance=advancement)

    def build_dataset(self, ligand_id_output_file, pdb_id_output_file, mapping_output_file):
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
        dataset_logger.info("Fetching all PDB IDs, that contain ligands...")
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

        dataset_logger.info(f"Found {len(entries_with_ligand)} entries with ligands in the Protein Data Bank.")

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
        
        dataset_logger.info(f"Found {len(mapping)} valid ligand IDs.")

        # export to files
        ligand_ids = list(mapping.keys())
        ligand_ids.sort()
        with open(ligand_id_output_file, 'w') as f:
            f.write("\n".join(ligand_ids))
        dataset_logger.info(f"Exported {len(ligand_ids)} ligand IDs to {ligand_id_output_file}.")

        # export PDB IDs
        pdb_ids = set()
        for i in mapping.values():
            pdb_ids.update(i)
        pdb_ids = list(pdb_ids)
        pdb_ids.sort()
        with open(pdb_id_output_file, 'w') as f:
            f.write("\n".join(pdb_ids))
        dataset_logger.info(f"Exported {len(pdb_ids)} PDB IDs to {pdb_id_output_file}.")

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
        dataset_logger.info(f"Exported mapping to {mapping_output_file}.")

        if invalid_ligand_ids:
            dataset_logger.warning(f"{len(invalid_ligand_ids)} invalid ligand IDs found: {', '.join(sorted(invalid_ligand_ids))}")

        return mapping
                


        

        
@click.command()
@click.argument('ligand_id_output_file', 
                type=click.Path(exists=False))
@click.argument('pdb_id_output_file', 
                type=click.Path(exists=False))
@click.argument('mapping_output_file',
                type=click.Path(exists=False))
# @click.option('--log_file', default='ligand_inspector.log', help='Path to the log file.')
# @click.option('--logging_level', default='INFO', help='Logging level.')
def main(ligand_id_output_file, pdb_id_output_file, mapping_output_file):
    """
    Main function to build the dataset.

    Args:
        input_file (str): path to the input file
        output_file (str): path to the output file

    Returns:
        None
    """
    builder = DatasetBuilder()
    builder.build_dataset(ligand_id_output_file, pdb_id_output_file, mapping_output_file)

if __name__ == "__main__":
    main()