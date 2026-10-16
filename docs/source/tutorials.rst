.. _tutorials:

=========
Tutorials
=========

Welcome to the `ligandqc` tutorials. This project is built to inspect the
quality of ligand structures against standard quality metrics, compare the outcome
across datasets, and export useful summaries for downstream analysis.

The project is organized around a few core modules:

* :mod:`ligandqc.inspect` - fetches ligand-quality data and evaluates it
* :mod:`ligandqc.download` - downloads the ideal ligand structure files
* :mod:`ligandqc.report` - generates summary tables and export reports
* :mod:`ligandqc.compare` - compares two or more datasets
* :mod:`ligandqc.plot` - Venn and Pareto plotting helpers

This tutorial covers the main workflow: inspect data, evaluate it, generate a
report, and compare datasets.

.. ATTENTION::
   Before starting, make sure the project is installed. See :doc:`installation`.
   You can also check the CLI with ``ligandqc --help``.


Quick overview of the project
=================================

The code follows a simple pattern:

1. Collect ligand IDs and optionally PDB IDs.
2. Query the Protein Data Bank for ligand quality data.
3. Evaluate each metric with the project-specific thresholds.
4. Save a CSV, YAML, JSON, HTML, or Excel result file.
5. Generate a report or plots.
6. Compare multiple datasets if needed.

The main public command-line entry points are:

* ``ligandqc inspect``
* ``ligandqc download``
* ``ligandqc plot``
* ``ligandqc report``
* ``ligandqc compare``


1. Inspect the quality of a single ligand in a single pdb file
==============================================================

In this case, we want to explore the quality of a single ligand occurance on a specific protein chain in a known pdb file. The ``only_SOI`` flag makes sure that the ligand of the specified ID is "subject of interest". This can be important as there are several cases, where ligands are present in PDB files as a crystallization agent or for other reasons. This flag makes sure that the ligand in question was flagged by the publishers of the file as "subject of interest".



.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            cat > ligand_ids.txt << EOF
            4WI
            EOF

            cat > pdb_ids.txt << EOF
            7TLL
            EOF

            ligandqc inspect \
                --ligand-id-file ligand_ids.txt \
                --output-file output.yaml \
                --pdb-constraint-file pdb_ids.txt \
                --only-chain "A" \
                --only-SOI

        This fetches the relevant ligand-quality metrics, evaluates them, and saves the
        result to ``output.yaml``. Changing the output extension will switch the export
        format, for example to CSV or HTML.

    .. tab-item:: Python
        :sync: python

        The core class is ``LigandInspector`` in :mod:`ligandqc.inspect`.
        It performs the data collection and metric evaluation.
        
        .. code-block:: python

            from ligandqc.inspect import LigandInspector

            ligand_ids = ["4WI"]
            pdb_ids = ["7TLL"]

            inspector = LigandInspector()
            df = inspector.inspect(
                ligand_ids=ligand_ids,
                pdb_ids=pdb_ids,
                only_chain="A",
                only_SOI=True,
            )
            df = inspector.evaluate(df)

            with open("output.yaml", "w") as f:
                yaml.dump(df.to_dict('records'), f)


which returns this file:

.. card:: output yaml content

    .. literalinclude:: ./_static/tutorial_files/tut1/output.yaml
        :language: yaml


2. Inspect the quality of a single ligand across all available PDB files
========================================================================

For this tutorial, we investigate the quality of Nirmatrelvir (4WI) across all available files in the Protein Data Bank. Notice how there are no pdb files specified. If that is the case, ``ligandqc`` searches the protein databank for all PDB files that contain this ligand.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            cat > ligands.txt << EOF
            4WI
            EOF

            ligandqc inspect \
                --ligand-id-file ligands.txt \
                --output-file quality.csv \
                --only-SOI

            ligandqc plot \
                --input-file quality.csv \
                --output-file quality.png

        The output file contains one row per ligand/PDB pair, with a column for each
        quality metric and the evaluation status for each metric.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            from ligandqc.inspect import LigandInspector

            ligand_ids = ["4WI"]

            inspector = LigandInspector()
            df = inspector.inspect(ligand_ids=ligand_ids, 
                                   only_SOI=True)
            df = inspector.evaluate(df)

            df.to_csv("quality.csv", index=False)

            inspector.plot(df, "quality.png")

        This is the main analysis loop used by the project.


so the generated data is 

.. card:: Generated Data

    .. literalinclude:: ./_static/tutorial_files/tut2/quality.csv

    .. figure:: ./_static/tutorial_files/tut2/quality.png
        :width: 600px

        Plot of ligand quality for Nirmatrelvir across different files in the PDB.



Generate a report
-----------------

The report module in :mod:`ligandqc.report` converts the evaluated data
into summary tables. It is designed to produce a compact overview of the measure
counts and percentages across all evaluation columns.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            ligandqc report \
                --input quality.csv \
                --output terminal

            ligandqc report \
                --input quality.csv \
                --output .


        This prints a summary table with the category distribution and the metric breakdown to the terminal, and the second command writes csv and html files containing the tables to the current directory.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            import pandas as pd
            from ligandqc.report import ReportGenerator

            df = pd.read_csv("quality.csv")
            report = ReportGenerator(df)
            report.generate_report(".")


which writes the following table to the terminal and/or the current directory.

.. card::

    .. literalinclude:: ./_static/tutorial_files/tut2/terminal_output.txt


3. Inspect a dataset of PDB IDs where you don't know the ligand id
==================================================================

A very common case that can happen when authors do a redocking study and their published data does not contain the ligand information but only the PDB ids. Because it's apparently very clear which ligands they mean. (ficticious scenario).

In that case, no ligand ids are provided, but the pdb ids instead must be.
If none of the two are provided, the program will fail.

In such cases, we recommend to put the ``only_SOI`` flag, as otherwise, crystallizing agents like DMS and ions or metals end up in the dataset.


.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            cat > pdb_ids.txt << EOF
            5REK
            7TLL
            5REJ
            EOF

            ligandqc inspect \
                --pdb-constraint-file pdb_ids.txt \
                --output-file ligands_quality.csv \
                --only-SOI
            
            ligandqc plot \
                --input-file ligands_quality.csv \
                --output-file ligands_quality.png

        The output file contains one row per ligand/PDB pair, with a column for each
        quality metric and the evaluation status for each metric.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            from ligandqc.inspect import LigandInspector

            pdb_ids = ["5REK", "5REJ", "7TLL"]

            inspector = LigandInspector()
            df = inspector.inspect(pdb_ids=pdb_ids, only_SOI=True)
            df = inspector.evaluate(df)

            df.to_csv("quality.csv", index=False)

            inspector.plot(df, "quality.png")


which should result in a file like this:

.. card::

    .. literalinclude:: ./_static/tutorial_files/tut3/quality.csv

and a plot like this:

.. card::

    .. figure:: ./_static/tutorial_files/tut3/quality.png
        :width: 600px

        Plot of ligand quality for the dataset of PDB files.



4. Inspect the entire PDB
=========================

A less common case where you would want to inspect the quality metrics of the entire Protein Data Bank. Here, we first generate files containing all ligand ids and all pdb ids found in the RCSB, then we use the typical commands, that we already know.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            ligandqc get-all-ids \
                -ol ligands.txt \
                -op pdb-ids.txt \
                -om mapping.txt 
            
            ligandqc inspect \
                --input-file ligands.txt \
                --output-file quality.csv \
                --only-SOI
            
            ligandqc plot \
                --input-file quality.csv \
                --output-file quality.png
    

.. card::

    .. figure:: ./_static/tutorial_files/tut4/quality.png
        :width: 600px

        Plot of ligand quality for the entire PDB.



5. Download ligand coordinates
==============================

The download command is used when you want the ideal 3D coordinates of a ligand.
This is useful for structural inspection or downstream modeling work.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            echo 4WI > ligands.txt

            ligandqc download \
                --ligand-id-file ligands.txt \
                --output-directory .

        This saves a structure file for the requested ligand IDs to the current directory. In this case it's only the one of Nirmatrelvir.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            from ligandqc.LigandDownloader import LigandDownloader

            downloader = LigandDownloader()
            downloader.download_ligand("4WI", "4WI.sdf")

which grabs the sdf file:

.. card::

    .. literalinclude:: ./_static/tutorial_files/tut5/4WI.sdf

6. Compare multiple datasets
============================

The compare workflow is designed for comparing two or more datasets. It creates a
Venn diagram of overlap, comparative Pareto plots, and summary tables. In this example, we are comparing to datasets of ligands for SARS-CoV-2 Mpro, which partially overlap.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            cat > dataset_a_ligand_ids.txt << EOFF
            T6M
            3XI
            T2Y
            83W
            T1V
            4W8
            HER
            OF9
            0BO
            H6F
            NOL
            EOFF

            cat > dataset_b_ligand_ids.txt << EOFF
            4W8
            HER
            OF9
            0BO
            H6F
            NOL
            T71
            MWN
            KMF
            8UI
            VHJ
            T9P
            Y48
            XTP
            H6L
            N2X
            SVL
            YVG
            N1U
            VHM
            8I7
            EOFF

            ligandqc inspect \
                --input-file dataset_a_ligand_ids.txt \
                --output-file dataset_a_quality.csv \
                --only-SOI

            ligandqc inspect \
                --input-file dataset_b_ligand_ids.txt \
                --output-file dataset_b_quality.csv \
                --only-SOI

            ligandqc compare \
                dataset_a_quality.csv dataset_b_quality.csv \
                --output-dir . \
                --labels dataset_a,dataset_b

        This saves a set of comparison plots and tables to the output directory.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            import pandas as pd
            from ligandqc.compare import DatasetComparison

            df_a = pd.read_csv("dataset_a.csv")
            df_b = pd.read_csv("dataset_b.csv")

            comparison = DatasetComparison(df_a, df_b, labels=["dataset_a", "dataset_b"])
            comparison.compare_datasets("comparison_results")

This will generate the following output plots in the specified output directory.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut6/venn_diagram_ligand_ids.png
        :width: 300px

        Venn diagram comparison of dataset overlap for ligand ids only.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut6/venn_diagram_pdb_id_ligand_id_pairs.png
        :width: 300px

        Venn diagram comparison of dataset overlap for PDB id, ligand id pairs

.. card::
    
    .. figure:: ./_static/tutorial_files/tut6/comparative_pareto_plot_good.png
        :width: 300px

        Comparative pareto plot for number of good evaluations.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut6/comparative_pareto_plot_ok.png
        :width: 300px

        Comparative pareto plot for number of ok evaluations.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut6/comparative_pareto_plot_bad.png
        :width: 300px

        Comparative pareto plot for number of bad evaluations.




Recommended workflow
====================

For most real-world use, the recommended sequence is:

1. Prepare a ligand list.
2. Run inspection.
3. Evaluate the output.
4. Generate a report.
5. Create plots for quick visual checks.
6. Compare with another dataset if needed.

This keeps the workflow reproducible and makes it easy to compare metric quality
across different ligand sets.

For additional support, see :doc:`cli`, :doc:`faq`, :doc:`known_issues`, and
:doc:`how_to_get_help`. 