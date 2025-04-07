# Webscraping Ligand Quality

In this project, we will scrape the ligand quality data from the Protein Data Bank using Python.

## Installation

### From Source
First, clone this repository:

```bash
git clone ...
cd ligand_quality
```

Then, create a new python virtual environment using your favorite tool and install the required packages using poetry.
The following command will create a new virtual environment using micromamba and install the packages:

```bash
micromamba create -n ligand_quality -c conda-forge python==3.12 poetry
micromamba activate ligand_quality
poetry install
```


