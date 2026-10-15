"""
Command Line Interface for Ligand Inspector
===========================================

This is the main entry point for the Ligand Inspector command line interface.

"""


import click
from ligand_inspector import __version__
from ligand_inspector.inspect import *
from ligand_inspector.download import download
from ligand_inspector.report import report
from ligand_inspector.inspect import inspect, plot
from ligand_inspector.compare import compare





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

    
if __name__ == "__main__":
    main()