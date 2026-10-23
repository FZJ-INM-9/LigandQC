.. ligandqc documentation master file, created by
   sphinx-quickstart on Tue May 12 14:06:04 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _repository: https://github.com/FZJ-INM-9/LigandQC
.. _documentation: https://ligandqc.readthedocs.io/en/latest/

##############################################
Welcome to **ligandqc** documentation!
##############################################

LigandQC is a python package with a command line interface (CLI) to check the quality 
of one or more ligands in protein structures via the RSCB API. The provided quality metrics are
categorized into "Good", "OK" and "Bad" based on rule of thumbs or statistical considerations for 
each metric.

The CLI provides commands to retrieve and evaluate the quality metrics as well as generate plots 
and reports from the generated data. The tool can be used to download ligand structures from the
Cambridge Structural Database (CSD) and compare the quality metrics and dataset overlap of up to 
three datasets at once. A utility command is provided to analyze the entire PDB as a whole.

The tool is intended to be used by researchers and scientists who want to perform analyses or 
simulations on a dataset of protein-ligand complexes and want to quickly inspect the quality of 
the ligands in their curated dataset.
Tutorials and example analyses on how to use the python package or the CLI are provided in the 
repository_ and the documentation_.


|:sparkles:| Features
=====================

* Scrape ligand quality data from the Protein Data Bank using a command line interface (CLI)
* Filter for Subject of Interest (SOI) ligands in a PDB file
* Filter for specific chains in a PDB file
* Save the data to a csv file for further analysis
* Create plots of the ligand quality data across multiple PDB files
* Download the ideal 3D coordinates of the ligand from the Cambridge Structural Database (CSD)
* Generate a report of the ligand quality data for a single ligand across multiple PDB files
* Compare the ligand quality data of one or more datasets


.. toctree::
   :hidden:
   :caption: About

   overview
   contributors
   license
   citation
   references

.. toctree::
   :hidden:
   :caption: Getting Started

   installation
   tutorials
   
.. toctree::
   :hidden:
   :caption: Help and Support
   
   faq
   how_to_get_help
   known_issues

.. toctree::
   :hidden:
   :caption: Source
   
   cli
   api/modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`