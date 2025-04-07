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
def single(pdb_id, ligand_id, output_file):
    """
    Inspect the quality of the ligand <LIGAND_ID> inside <PDB_ID>.pdb from the Protein Data Bank. The results will be saved to <OUTPUT_FILE>.csv
    """
    li = LigandInspector()
    df = li.inspect_single(pdb_id, ligand_id)
    df.to_csv(output_file, index=False)
    click.echo(f"Inspection results saved to {output_file}")

@main.command()
@click.argument('input_file')
@click.argument('output_file')
@click.option('--workers', default=10, help='Number of workers for batch inspection.')
def batch(input_file, output_file, workers):
    """
    Inspect a batch of ligands from a CSV file.
    """
    df = pd.read_csv(input_file)
    li = LigandInspector()
    results = li.inspect_batch(df, n_workers=workers)
    results.to_csv(output_file, index=False)
    click.echo(f"Batch inspection results saved to {output_file}")

if __name__ == "__main__":
    main()