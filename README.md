# :mag: Webscraping Ligand Quality

In this project, we will scrape the ligand quality data from the Protein Data Bank using Python.

> [!WARNING]
> This project is under active development and is not yet ready for production use. The code is not fully tested and may contain bugs. Use at your own risk.

# :inbox_tray: Installation

## From Source
First, clone this repository:

```bash
git clone https://jugit.fz-juelich.de/inm-9/ligand-quality-scraping.git
cd ligand-quality-scraping
```

Then, create a new python virtual environment using your favorite tool and install the required packages.
The following command will create a new virtual environment using micromamba and install the packages using poetry:

```bash
micromamba create -n ligand_quality -c conda-forge python==3.12 poetry
micromamba activate ligand_quality
poetry install
```

# Usage

Once the virtual environment is activated, you can run the script using the following command:

```bash
ligand_inspector --help
```

This will show you the available subcommands, arguments and options. There are two modes of operation: single and batch.
In the single mode, you can inspect a single ligand by providing its three letter ligand id and the PDB id of the structure it is bound to. The output will be saved to the provided csv file path.

```bash
ligand_inspector single 5REJ T1V output.csv
```

In the batch mode, you can provide a csv file containing two columns ("PDB ID" and "Ligand ID"). The script will read the file and inspect all rows. The output will be saved to the provided csv file path.

```bash
ligand_inspector batch input.csv output.csv
```

