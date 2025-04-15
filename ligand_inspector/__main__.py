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
from .LigandInspector import *
from . import __version__



def setup_logger(log_file, logging_level):
    # set logging level
    logging_level = logging_level.upper()
    if logging_level not in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
        click.echo("Invalid logging level. Please use DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
        return
    # set logging level
    logger.setLevel(logging_level)
    # set log file
    handler = logging.FileHandler(log_file)
    handler.setLevel(logging_level)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)



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

@main.command()
@click.argument('pdb_id',
                type=click.STRING)
@click.argument('ligand_id')
@click.argument('output_file',
                type=click.Path(exists=False),
                required=True)
@click.option('--ignore-other-entries', "-ign",
              is_flag=True,
              default=False,
              help="If set, entries in the scraped table that are not in the input file will be removed from the reported file.")
@click.option("--only-chain", "-oc",
                type=click.STRING,
                help="If set, only the specified chain will be inspected.")
@click.option('--logging-level', '-ll',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], case_sensitive=False),
              default='ERROR',
              help="Set the logging level; Choices are DEBUG, INFO, WARNING, ERROR, CRITICAL. Default is INFO.")
@click.option("--log-file", "-lf",
              type=click.Path(exists=False),
              default="ligand_inspector.log",
                help="Path to the log file. Default is ligand_inspector.log.")
@click.option("--keep-sdf", "-ks",
              is_flag=True,
              default=False,
              help="If set, the content of the SDF file will be kept in the output csv file in the column 'Ligand sdf'.")
def single(pdb_id, ligand_id, output_file, ignore_other_entries, only_chain, logging_level, log_file, keep_sdf):
    """
    Inspect the quality of the ligand <LIGAND_ID> inside <PDB_ID>.pdb from the Protein Data Bank. The results will be saved to <OUTPUT_FILE>. Supported formats are .csv, .xlsx, .json, .yaml, .html, and .txt.
    """
    setup_logger(log_file, logging_level)

    logger.info(f"Inspecting {pdb_id} {ligand_id}")
    li = LigandInspector()
    logger.info(f"Downloading information from PDB...")
    df = li.inspect_single(pdb_id, ligand_id, ignore_other_entries=ignore_other_entries, only_chain=only_chain)
    logger.info(f"Evaluating the provided data...")
    df = li.evaluate(df)
    # drop the 'Ligand sdf' column
    if not keep_sdf:
        df = df.drop(columns=['Ligand sdf'], errors='ignore')

    # Save the results to the specified output file
    logger.info(f"Saving results to {output_file}...")
    if output_file.endswith(".yaml") or output_file.endswith(".yml"):
        with open(output_file, 'w') as f:
            yaml.dump(df.iloc[0].to_dict(), f)
    elif output_file.endswith(".json"):
        df.to_json(output_file, orient='records', lines=True)
    elif output_file.endswith(".xlsx"):
        df.to_excel(output_file, index=False)
    elif output_file.endswith(".html"):
        df.to_html(output_file, index=False)
    elif output_file.endswith(".txt"):
        df.to_csv(output_file, sep='\t', index=False)
    elif output_file.endswith(".csv"):
        df.to_csv(output_file, index=False)
    elif output_file.endswith(".md"):
        df.to_markdown(output_file, index=False)
    else:
        click.echo("Unsupported file format. Please use .csv, .xlsx, .json, .yaml, .yml, .html or .txt.")
        return
    click.echo(f"Inspection results saved to {output_file}")

@main.command()
@click.argument('input_file')
@click.argument('output_file')
@click.option('--ignore-other-entries', "-ign",
              is_flag=True,
              help="If set, entries in the scraped table that are not in the input file will be removed from the reported file.")
@click.option('--workers', default=10, help='Number of workers for batch inspection.')
@click.option('--only-chain', "-oc",
                type=click.STRING,
                help="If set, only the specified chain will be inspected.")
@click.option('--logging-level', '-ll',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], case_sensitive=False),
              default='INFO',
              help="Set the logging level; Choices are DEBUG, INFO, WARNING, ERROR, CRITICAL. Default is INFO.")
@click.option("--log-file", "-lf",
              type=click.Path(exists=False),
              default="ligand_inspector.log",
                help="Path to the log file. Default is ligand_inspector.log.")
@click.option("--keep-sdf", "-ks",
              is_flag=True,
              default=False,
              help="If set, the content of the SDF file will be kept in the output csv file in the column 'Ligand sdf'.")
def batch(input_file, output_file, ignore_other_entries, workers, only_chain, logging_level, log_file, keep_sdf):
    """
    Inspect a batch of ligands as specified in <INPUT_FILE>. The input file should be a CSV file with columns 'PDB ID' and 'Ligand ID'. The results will be saved to <OUTPUT_FILE>. Supported formats are .csv, .xlsx, .json, .yaml, .html, and .txt.
    """
    setup_logger(log_file, logging_level)
    logger.info(f"Inspecting ligands in {input_file}")

    if not os.path.exists(input_file):
        logger.critical(f"Input file {input_file} does not exist.")
        return
    
    if not os.path.isfile(input_file):
        logger.critical(f"Input file {input_file} is not a file.")
        return
    
    logger.info(f"Reading input file {input_file}...")
    df = pd.read_csv(input_file)

    li = LigandInspector()
    logger.info(f"Downloading information from PDB...")
    results = li.inspect_batch(df, n_workers=workers, ignore_other_entries=ignore_other_entries, only_chain=only_chain)
    
    logger.info(f"Evaluating the provided data...")
    results = li.evaluate(results)

    # drop the 'Ligand sdf' column
    if not keep_sdf:
        results = results.drop(columns=['Ligand sdf'], errors='ignore')

    # Save the results to the specified output file
    logger.info(f"Saving results to {output_file}...")
    if output_file.endswith(".yaml") or output_file.endswith(".yml"):
        with open(output_file, 'w') as f:
            yaml.dump(results.to_dict(orient='records'), f)
    elif output_file.endswith(".json"):
        results.to_json(output_file, orient='records', lines=True)
    elif output_file.endswith(".xlsx"):
        results.to_excel(output_file, index=False)
    elif output_file.endswith(".html"):
        results.to_html(output_file, index=False)
    elif output_file.endswith(".txt"):
        results.to_csv(output_file, sep='\t', index=False)
    elif output_file.endswith(".csv"):
        results.to_csv(output_file, index=False)
    elif output_file.endswith(".md"):
        results.to_markdown(output_file, index=False)
    else:
        click.echo("Unsupported file format. Please use .csv, .xlsx, .json, .yaml, .yml, .html or .txt.")
        return

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

if __name__ == "__main__":
    main()