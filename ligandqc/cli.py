"""
Command Line Interface for Ligand Inspector
===========================================

This is the main entry point for the Ligand Inspector command line interface.

"""


import click
from ligandqc import __version__
from ligandqc.inspect import *
from ligandqc.download import download
from ligandqc.report import report
from ligandqc.inspect import inspect, plot
from ligandqc.compare import compare
from ligandqc.rcsb import get_all_ids





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
main.add_command(inspect)
main.add_command(compare)
main.add_command(download)
main.add_command(plot)
main.add_command(get_all_ids)

    
if __name__ == "__main__":
    main()