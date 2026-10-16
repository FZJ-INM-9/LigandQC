.. _installation:

============
Installation
============

This page describes ways to install the `ligandqc` python package.

From scratch
============

Use this approach if you want a local development setup from source.

1. Clone the repository.

   .. code-block:: bash

      git clone <your-repository-url>
      cd ligandqc

2. Create and activate a virtual environment using the tool of your choice (e.g. micromamba or conda) with the basic dependencies (check out the mentioned file)

   .. code-block:: bash

      micromamba create -f micromamba.yaml
      

3. Install build tooling and the package in editable mode.

   .. code-block:: bash

      poetry install

4. Verify the installation.

   .. code-block:: bash

      micromamba activate ligandqc
      enunciate --help


Via package manager
===================

Install `ligandqc` via a package manager.

Via pip
-------

Install the latest released version:

.. code-block:: bash

   pip install ligandqc

Via poetry
----------

... to be added ...


Via conda
---------

... to be added ...