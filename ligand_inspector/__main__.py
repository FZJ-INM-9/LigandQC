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

from .LigandInspector import *

@click.group()
def main():
    """Command Line Interface for Ligand Inspector."""
    pass

@main.command()
@click.argument('pdb_id',
                type=click.STRING)
@click.argument('ligand_id')
@click.argument('output_file')
@click.option('--ignore-other-entries', "-ign",
              is_flag=True,
              default=False,
              help="If set, entries in the scraped table that are not in the input file will be removed from the reported file.")
@click.option("--only-chain", "-oc",
                type=click.STRING,
                help="If set, only the specified chain will be inspected.")
def single(pdb_id, ligand_id, output_file, ignore_other_entries, only_chain):
    """
    Inspect the quality of the ligand <LIGAND_ID> inside <PDB_ID>.pdb from the Protein Data Bank. The results will be saved to <OUTPUT_FILE>. Supported formats are .csv, .xlsx, .json, .yaml, .html, and .txt.
    """
    li = LigandInspector()
    df = li.inspect_single(pdb_id, ligand_id, ignore_other_entries=ignore_other_entries, only_chain=only_chain)
    df = li.evaluate(df)
    # drop the 'Ligand sdf' column
    df = df.drop(columns=['Ligand sdf'], errors='ignore')
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
def batch(input_file, output_file, ignore_other_entries, workers, only_chain):
    """
    Inspect a batch of ligands as specified in <INPUT_FILE>. The input file should be a CSV file with columns 'PDB ID' and 'Ligand ID'. The results will be saved to <OUTPUT_FILE>. Supported formats are .csv, .xlsx, .json, .yaml, .html, and .txt.
    """
    df = pd.read_csv(input_file)
    li = LigandInspector()
    results = li.inspect_batch(df, n_workers=workers, ignore_other_entries=ignore_other_entries, only_chain=only_chain)
    results = li.evaluate(results)
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
    click.echo(f"Batch inspection results saved to {output_file}")

if __name__ == "__main__":
    main()