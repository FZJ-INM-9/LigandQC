.. _tutorials:

=========
Tutorials
=========

Welcome to the `ligand_inspector` tutorials. This project is built to inspect the
quality of ligand structures against standard quality metrics, compare the outcome
across datasets, and export useful summaries for downstream analysis.

The project is organized around a few core modules:

* :mod:`ligand_inspector.inspect` - fetches ligand-quality data and evaluates it
* :mod:`ligand_inspector.download` - downloads the ideal ligand structure files
* :mod:`ligand_inspector.report` - generates summary tables and export reports
* :mod:`ligand_inspector.compare` - compares two or more datasets
* :mod:`ligand_inspector.plot` - Venn and Pareto plotting helpers

This tutorial covers the main workflow: inspect data, evaluate it, generate a
report, and compare datasets.

.. ATTENTION::
   Before starting, make sure the project is installed. See :doc:`installation`.
   You can also check the CLI with ``ligand_inspector --help``.


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

* ``ligand_inspector inspect``
* ``ligand_inspector download``
* ``ligand_inspector plot``
* ``ligand_inspector report``
* ``ligand_inspector compare``


1. Inspect one ligand
=====================

The core class is ``LigandInspector`` in :mod:`ligand_inspector.inspect`.
It performs the data collection and metric evaluation.

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

            ligand_inspector inspect \
                -i ligand_ids.txt \
                -o output.yaml \
                --pdb-constraint-file pdb_ids.txt \
                --only-chain "A"

        This fetches the relevant ligand-quality metrics, evaluates them, and saves the
        result to ``output.yaml``. Changing the output extension will switch the export
        format, for example to CSV or HTML.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            from ligand_inspector.inspect import LigandInspector

            ligand_ids = ["4WI"]
            pdb_ids = ["7TLL"]

            inspector = LigandInspector()
            df = inspector.inspect_batch(
                ligand_ids,
                pdb_ids=pdb_ids,
                only_chain="A",
                only_SOI=True,
            )
            df = inspector.evaluate(df)

            df.to_csv("output.yaml", index=False)

        The resulting dataframe contains the raw metric values and evaluation columns such
        as ``Real space R factor evaluation`` and ``overall evaluation``.


which returns this file

.. card:: output yaml content

    .. literalinclude:: ./_static/tutorial_files/tut1/output.yaml
        :language: yaml


2. Inspect a dataset with a plot
================================

For larger studies, provide a multiple ligand IDs and optionally a PDB
constraint file. For this tutorial, we investigate the quality of Nirmatrelvir across different files.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            cat > ligands.txt << EOF
            4WI
            EOF

            ligand_inspector inspect \
                -i ligands.txt \
                -o ligands_quality.csv \
                --only-SOI
            
            ligand_inspector plot \
                -i ligands_quality.csv \
                -o ligands_quality.png

        The output file contains one row per ligand/PDB pair, with a column for each
        quality metric and the evaluation status for each metric.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            from ligand_inspector.inspect import LigandInspector

            ligand_ids = ["4WI"]

            inspector = LigandInspector()
            df = inspector.inspect_batch(ligand_ids, only_SOI=True)
            df = inspector.evaluate(df)

            df.to_csv("ligands_quality.csv", index=False)

            inspector.plot(df, "ligands_quality.png")

        This is the main analysis loop used by the project.


so the generated data is 

.. card:: Generated Data

    .. literalinclude:: ./_static/tutorial_files/tut2/ligands_quality.csv

    .. figure:: ./_static/tutorial_files/tut2/ligands_quality.png
        :width: 600px

        Plot of ligand quality for Nirmatrelvir across different files in the PDB.



Generate a report
-----------------

The report module in :mod:`ligand_inspector.report` converts the evaluated data
into summary tables. It is designed to produce a compact overview of the measure
counts and percentages across all evaluation columns.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            ligand_inspector report \
                --input output.csv \
                --output terminal


        This prints a terminal summary table with the category distribution and the metric
        breakdown, and the second command writes the report tables to disk.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            import pandas as pd
            from ligand_inspector.report import ReportGenerator

            df = pd.read_csv("output.csv")
            report = ReportGenerator(df)
            report.generate_report("report_output")


which writes the following table to the terminal

.. card::

    .. literalinclude:: ./_static/tutorial_files/tut2/terminal_output.txt


3. Download ligand coordinates
==============================

The download command is used when you want the ideal 3D coordinates of a ligand.
This is useful for structural inspection or downstream modeling work.

.. tab-set::
    :sync-group: tutorial

    .. tab-item:: CLI
        :sync: cli

        .. code-block:: bash

            echo 4WI > ligands.txt

            ligand_inspector download \
                -i ligands.txt \
                -o .

        This saves a structure file for the requested ligand IDs to the current directory. In this case it's only the one of Nirmatrelvir.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            from ligand_inspector.LigandDownloader import LigandDownloader

            downloader = LigandDownloader()
            downloader.download_ligand("4WI", "4WI.sdf")



4. Compare multiple datasets
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

            ligand_inspector inspect \
                -i dataset_a_ligand_ids.txt \
                -o dataset_a_quality.csv \
                --only-SOI

            ligand_inspector inspect \
                -i dataset_b_ligand_ids.txt \
                -o dataset_b_quality.csv \
                --only-SOI

            ligand_inspector compare \
                dataset_a_quality.csv dataset_b_quality.csv \
                --output-dir . \
                --labels dataset_a,dataset_b

        This saves a set of comparison plots and tables to the output directory.

    .. tab-item:: Python
        :sync: python

        .. code-block:: python

            import pandas as pd
            from ligand_inspector.compare import DatasetComparison

            df_a = pd.read_csv("dataset_a.csv")
            df_b = pd.read_csv("dataset_b.csv")

            comparison = DatasetComparison(df_a, df_b, labels=["dataset_a", "dataset_b"])
            comparison.compare_datasets("comparison_results")

This will generate the following output plots

.. card::
    
    .. figure:: ./_static/tutorial_files/tut4/venn_diagram_ligand_ids.png
        :width: 300px

        Venn diagram comparison of dataset overlap for ligand ids only.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut4/venn_diagram_pdb_id_ligand_id_pairs.png
        :width: 300px

        Venn diagram comparison of dataset overlap for PDB id, ligand id pairs

.. card::
    
    .. figure:: ./_static/tutorial_files/tut4/comparative_pareto_plot_good.png
        :width: 300px

        Comparative pareto plot for number of good evaluations.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut4/comparative_pareto_plot_ok.png
        :width: 300px

        Comparative pareto plot for number of ok evaluations.

.. card::
    
    .. figure:: ./_static/tutorial_files/tut4/comparative_pareto_plot_bad.png
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