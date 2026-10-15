.. ligand_inspector documentation master file, created by
   sphinx-quickstart on Tue May 12 14:06:04 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##############################################
Welcome to **ligand_inspector** documentation!
##############################################

**ligand_inspector** is a python package to inspect the quality metrics of single PDB files or a dataset of ligands.

|:sparkles:| Features

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