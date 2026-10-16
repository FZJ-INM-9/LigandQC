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
888      d8b                                 888  .d88888b.   .d8888b.  
888      Y8P                                 888 d88P" "Y88b d88P  Y88b 
888                                          888 888     888 888    888 
888      888  .d88b.   8888b.  88888b.   .d88888 888     888 888        
888      888 d88P"88b     "88b 888 "88b d88" 888 888     888 888        
888      888 888  888 .d888888 888  888 888  888 888 Y8b 888 888    888 
888      888 Y88b 888 888  888 888  888 Y88b 888 Y88b.Y8b88P Y88b  d88P 
88888888 888  "Y88888 "Y888888 888  888  "Y88888  "Y888888"   "Y8888P"  
                  888                                   Y8b             
             Y8b d88P                                                   
              "Y88P"                                                    
              
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