"""
Command Line Interface for Ligand Inspector
===========================================

This is the main entry point for the Ligand Inspector command line interface.
Two commands are provided:
1. single: Inspect a single ligand.
2. batch: Inspect a batch of ligands from a CSV file.

Usage:
    ligand_inspector single <pdb_id> <ligand_id> <output_file>
    ligand_inspector batch <input_file> <output_file> --workers <n_workers>

"""


import click
import pandas as pd
import yaml
import logging
from .inspect import *
import os
from . import __version__
from ligand_inspector.report import report





@click.group()
def main():
    """Command Line Interface for Ligand Inspector."""
    click.echo(f"""
▗▖   ▗▄▄▄▖ ▗▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄     ▗▄▄▄▖▗▖  ▗▖ ▗▄▄▖▗▄▄▖ ▗▄▄▄▖ ▗▄▄▖▗▄▄▄▖▗▄▖ ▗▄▄▖ 
▐▌     █  ▐▌   ▐▌ ▐▌▐▛▚▖▐▌▐▌  █      █  ▐▛▚▖▐▌▐▌   ▐▌ ▐▌▐▌   ▐▌     █ ▐▌ ▐▌▐▌ ▐▌
▐▌     █  ▐▌▝▜▌▐▛▀▜▌▐▌ ▝▜▌▐▌  █      █  ▐▌ ▝▜▌ ▝▀▚▖▐▛▀▘ ▐▛▀▀▘▐▌     █ ▐▌ ▐▌▐▛▀▚▖
▐▙▄▄▖▗▄█▄▖▝▚▄▞▘▐▌ ▐▌▐▌  ▐▌▐▙▄▄▀    ▗▄█▄▖▐▌  ▐▌▗▄▄▞▘▐▌   ▐▙▄▄▖▝▚▄▄▖  █ ▝▚▄▞▘▐▌ ▐▌
                                                                                
created by: Christian Bergeest, Forschungszentrum Jülich GmbH, Germany, 2025
version: {__version__}               

""")

main.add_command(report)
# @main.command()
# @click.argument('ligand_id')
# @click.argument('output_file',
#                 type=click.Path(exists=False),
#                 required=True)
# @click.option('--pdb-id', '-pdb',
#                 type=click.STRING,
#                 multiple=True,
#                 help="If set, the specified PDB ID will be used. Otherwise, the ligand ID will be used to search for all PDB IDs. Can be used multiple times to specify multiple PDB IDs.")
# @click.option("--only-chain", "-oc",
#                 type=click.STRING,
#                 help="If set, only the specified chain will be inspected.")
# @click.option('--logging-level', '-ll',
#               type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], case_sensitive=False),
#               default='ERROR',
#               help="Set the logging level; Choices are DEBUG, INFO, WARNING, ERROR, CRITICAL. Default is INFO.")
# @click.option("--log-file", "-lf",
#               type=click.Path(exists=False),
#               default="ligand_inspector.log",
#                 help="Path to the log file. Default is ligand_inspector.log.")
# def single(pdb_id, ligand_id, output_file, only_chain, logging_level, log_file):
#     """
#     Inspect the quality of the ligand <LIGAND_ID> inside <PDB_ID>.pdb from the Protein Data Bank. The results will be saved to <OUTPUT_FILE>. Supported formats are .csv, .xlsx, .json, .yaml, .html, and .txt.
#     """
#     setup_logger(log_file, logging_level)

#     logger.info(f"Inspecting {pdb_id} {ligand_id}")
#     li = LigandInspector()
#     logger.info(f"Downloading information from PDB...")
#     df = li.inspect_single(pdb_id, ligand_id, only_chain=only_chain)
#     logger.info(f"Evaluating the provided data...")
#     df = li.evaluate(df)

#     # Save the results to the specified output file
#     logger.info(f"Saving results to {output_file}...")
#     if output_file.endswith(".yaml") or output_file.endswith(".yml"):
#         with open(output_file, 'w') as f:
#             yaml.dump(df.iloc[0].to_dict(), f)
#     elif output_file.endswith(".json"):
#         df.to_json(output_file, orient='records', lines=True)
#     elif output_file.endswith(".xlsx"):
#         df.to_excel(output_file, index=False)
#     elif output_file.endswith(".html"):
#         df.to_html(output_file, index=False)
#     elif output_file.endswith(".txt"):
#         df.to_csv(output_file, sep='\t', index=False)
#     elif output_file.endswith(".csv"):
#         df.to_csv(output_file, index=False)
#     elif output_file.endswith(".md"):
#         df.to_markdown(output_file, index=False)
#     else:
#         click.echo("Unsupported file format. Please use .csv, .xlsx, .json, .yaml, .yml, .html or .txt.")
#         return
#     click.echo(f"Inspection results saved to {output_file}")


@main.command()
@click.argument('input_file')
@click.argument('output_file')
def plot(input_file, output_file):
    """
    Plot the results from <INPUT_FILE> and save to <OUTPUT_FILE>.
    """
    logger.info(f"Plotting results from {input_file}...")
    df = pd.read_csv(input_file)
    li = LigandInspector()
    li.plot(df, output_file)
    logger.info(f"Plot saved to {output_file}")


@main.command()
@click.argument('ligand_id_file')
@click.argument('output_directory')
def download(ligand_id_file, output_directory):
    """
    Download the ligands from <LIGAND_ID_FILE> and save to <OUTPUT_DIRECTORY>.
    """
    downloader = LigandDownloader()

    with open(ligand_id_file, 'r') as f:
        ligand_ids = f.read().splitlines()
        ligand_ids = [ligand_id.strip() for ligand_id in ligand_ids]
        while "" in ligand_ids:
            ligand_ids.remove("")
        
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    ligands = downloader.download_multiple_ligands(ligand_ids)

    for ligand_data in ligands:
        # Get the ligand ID from the data
        if ligand_data is None:
            continue
        ligand_id = ligand_data.split()[0].strip()

        ligand_file = os.path.join(output_directory, f"{ligand_id}.sdf")
        with open(ligand_file, 'w') as f:
            f.write(ligand_data)
        logger.info(f"Downloaded {ligand_id} to {ligand_file}")

    logger.info(f"Downloaded {len(ligand_ids)} ligands to {output_directory}")

    
if __name__ == "__main__":
    main()