# :mag: Webscraping Ligand Quality

In this project, we will scrape the ligand quality data from the Protein Data Bank using Python.

> [!WARNING] :construction: Construction Warning :construction:
> This project is under active development and is not yet ready for production use. The code is not fully tested and may contain bugs. Use at your own risk.

# :inbox_tray: Installation
To install the package, use one of the following methods.

## :card_file_box: From Source
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

## :rocket: Package Manager
No package manager support yet.

# :joystick: Usage

Once the virtual environment is activated, you can use the tool as follows:

```bash
ligand_inspector --help
```

This will show you the available subcommands, arguments and options. There are two modes of operation: single and batch.
In the single mode, you can inspect a single ligand by providing the PDB ID of the structure it is bound to and the three letter Ligand ID. The output will be saved to the provided csv file path.

```bash
ligand_inspector single 5REJ T1V output.csv
```

In the batch mode, you can provide a csv file containing two columns ("PDB ID" and "Ligand ID"). The script will read the file and inspect all rows. The output will be saved to the provided csv file path.

```bash
ligand_inspector batch input.csv output.csv
```

# :brain: Rationale
The goal of this project is to provide a simple and easy to use tool to inspect the quality of ligands in protein structures. The tool will scrape the data from the Protein Data Bank and provide a summary of the quality of the ligand. The tool will also provide a way to save the data to a csv file for further analysis.
The tool is intended to be used by researchers and scientists who want to quickly inspect the quality of ligands in protein structures.
Currently, in the Protein Data Bank there are multiple metrics that are used to assess the quality of ligands in protein structures. These metrics are not always easy to find and compare, especially for large datasets.
The metrics provided are:
- Real Space R-factor (RSR)
- Real Space Correlation Coefficient (RSCC)
- RMSZ-bond-length
- RMSZ-bond-angle
- Outliers of bond length
- Outliers of bond angle
- Average occupancy
- Atomic clashes
- Stereo-chemical errors

The above mentioned metrics alone are not sufficient to assess the quality of a ligand. Together they can provide a good overview.
In this tool, the ligand stuctures are categorized into "Good", "OK" and "Bad" based rule of thumbs for each metric. An overall category is assigned based on the worst metric category obtained. "Outliers of bond length" and "Outliers of bond angle" have been summarized into "Structural outliers".

The rule of thumbs go as follows:

| metric                  | Good       | OK            | Bad        | citation                                  |
|-------------------------|------------|---------------|------------|-------------------------------------------|
| RSR                     | <= 0.4     | -             | > 0.4      | https://doi.org/10.1107/S2059798318002541 |
| RSCC                    | >= 0.95    | >= 0.8        | < 0.8      | https://doi.org/10.1107/S2059798318002541 |
| RMSZ-bond-length        | *          | <= 2.0        | > 2.0      | [Link to Mogul FAQ](https://support.ccdc.cam.ac.uk/support/solutions/articles/103000306395-what-is-the-z-score-listed-in-the-summary-statistics-for-each-3d-bond-and-valence-angle-search-result)<br>accessed 2024-04-14                             |
| RMSZ-bond-angle         | *          | <= 2.0        | > 2.0      | [Link to Mogul FAQ](https://support.ccdc.cam.ac.uk/support/solutions/articles/103000306395-what-is-the-z-score-listed-in-the-summary-statistics-for-each-3d-bond-and-valence-angle-search-result)<br>accessed 2024-04-14                             |
| Structural Outliers     | <= 5       | <= 15         | > 15       | -                                         |
| Average occupancy       | >= 0.95    | >= 0.5        | < 0.5      | -                                         |
| Atomic clashes          | 0          | <= 2          | > 2        | -                                         |
| Stereo-chemical errors  | 0          | -             | > 0        | https://doi.org/10.1016/j.str.2021.10.003 |

## :rotating_light: RMSZ values
The Z-value $Z_i$ describes the deviation of a value $X_i$ from the mean of a population $\mu$ in units of standard deviation $\sigma$.
The Z-value is calculated as follows:

```math
Z_i = \frac{X_i - \mu}{\sigma}
```

Where $N$ is the number of values, the $RMSZ$ value is calculated as follows:

```math
RMSZ = \sqrt{\frac{1}{N} \sum_{i=1}^{N} Z_i^2}
```

The Z-values are normalized to have a mean of 0, therefore the RMSZ value is equivalent to the *standard deviation* of the Z-values. 
Applying this to the RMSZ-bond-length and RMSZ-bond-angle metrics, we can see that the RMSZ value for the ligand is a measure of the deviation of the bond lengths/angles from the mean bond lengths/angles of the specific type.
To evaluate the provided RMSZ value, we need to know if the sample size (number of bonds/angles in a ligand) is large enough to prove statistical deviation from the underlying general population.
This is done by performing a series of **F-tests for equality of variances**. 

The F-test is a statistical test that compares the variances of two populations. The null hypothesis is that the variances are equal, and the alternative hypothesis is that the variances are not equal. The F-test returns a p-value, which is the probability of observing the data given that the null hypothesis is true. If the p-value is less than a certain threshold (usually 5%), we reject the null hypothesis and conclude that the variances are not equal.
For each ligand, the number of bonds/angles was determined by analysing the molecular structure of the ligand as obtained from the RCSB Ligand Expo using rdkit. 
A treshold-RMSZ value is determined, below which no statistical deviation can be proven.

Using a bisectional search, the treshold-RMSZ value is determined, below which no statistical deviation can be proven.
The size of the true population was chosen to be 1E7, which is large enough to be considered infinite compared to the number of bonds/angles typically found in a ligand which we assume to be in the range of 5-200.


# :alembic: Example
```bash
ligand_inspector single 7TLL 4WI output.yaml --only-chain "A" --ignore-other-entries
```

yields the following output:

```yaml
Atomic clashes: 0
Atomic clashes evaluation: good
Average occupancy: 1.0
Average occupancy evaluation: good
Chain ID: A
Identifier: 7TLL_4WI_A_401
Ligand ID: 4WI
Ligand occurence in PDB: 7
Model completeness: 1.0
Outliers of bond angle: 5
Outliers of bond length: 3
PDB ID: 7TLL
RMSZ-bond-angle: 1.23
RMSZ-bond-angle evaluation: ok
RMSZ-bond-length: 1.43
RMSZ-bond-length evaluation: ok
Ranking for geometry: 0.33
Ranking for goodness of fit: 0.71
Real space R factor: 0.09
Real space R factor evaluation: good
Real space correlation coefficient: 0.908
Real space correlation coefficient evaluation: ok
Stereochemical errors: 0
Stereochemical errors evaluation: good
Structural outliers: 8
Structural outliers evaluation: ok
overall evaluation: ok
```