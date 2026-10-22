.. _overview:

########
Overview
########


*******************
|:brain:| Rationale
*******************

The goal of this project is to provide a simple and easy to use tool to inspect the quality of ligands in protein structures. The tool will scrape the data from the Protein Data Bank and provide a summary of the quality of the ligand. The tool will also provide a way to save the data to a csv file for further analysis.
The tool is intended to be used by researchers and scientists who want to quickly inspect the quality of ligands in protein structures.
Currently, in the Protein Data Bank there are multiple metrics that are used to assess the quality of ligands in protein structures. These metrics are not always easy to find and compare, especially for large datasets.
The metrics provided are:

* Real Space R-factor (RSR)
* Real Space Correlation Coefficient (RSCC)
* RMSZ Bond Length
* RMSZ Bond Angle
* Outliers of bond length
* Outliers of bond angle
* Average occupancy
* Atomic clashes
* Stereo-chemical errors


The above mentioned metrics alone are not sufficient to assess the quality of a ligand. Together they can provide a good overview.
In this tool, the ligand stuctures are categorized into "Good", "OK" and "Bad" based rule of thumbs for each metric. An overall category is assigned based on the worst metric category obtained. The outlier rates for bond lengths and bond angles are evaluated independently as percentages of the total bond and angle counts.

The rule of thumbs go as follows:


.. list-table::
    :header-rows: 1

    * - metric
      - |CATEGORY_GOOD|
      - |CATEGORY_OK|
      - |CATEGORY_BAD|
      - reference
    * - RSR                        
      - :math:`\leq` |THRESH_RSR_GOOD|
      - non-existent      
      - :math:`>` |THRESH_RSR_GOOD|    
      - :footcite:t:`Smart.2018`
    * - RSCC
      - :math:`\geq` |THRESH_RSCC_GOOD|
      - :math:`\geq` |THRESH_RSCC_BAD|
      - :math:`<` |THRESH_RSCC_BAD|
      - :footcite:t:`Smart.2018`
    * - RMSZ Bond Length
      - see :ref:`rmsz_rationale`         
      - :math:`\leq` |THRESH_RMSZ_BOND_LENGTH_BAD|  
      - :math:`>` |THRESH_RMSZ_BOND_LENGTH_BAD|    
      - :footcite:t:`CCDC.2024` 
    * - RMSZ Bond Angle            
      - see :ref:`rmsz_rationale`         
      - :math:`\leq` |THRESH_RMSZ_BOND_ANGLE_BAD|
      - :math:`>` |THRESH_RMSZ_BOND_ANGLE_BAD|
      - :footcite:t:`CCDC.2024` 
    * - Outliers of bond length
      - :math:`\leq` |THRESH_OUTLIERS_BOND_LENGTH_GOOD|
      - :math:`\leq` |THRESH_OUTLIERS_BOND_LENGTH_BAD|
      - :math:`>` |THRESH_OUTLIERS_BOND_LENGTH_BAD|
      - none                    
    * - Outliers of bond angle
      - :math:`\leq` |THRESH_OUTLIERS_BOND_ANGLE_GOOD|
      - :math:`\leq` |THRESH_OUTLIERS_BOND_ANGLE_BAD|
      - :math:`>` |THRESH_OUTLIERS_BOND_ANGLE_BAD|
      - none                      
    * - Average Occupancy
      - :math:`\geq` |THRESH_AVG_OCCUPANCY_GOOD|
      - :math:`\geq` |THRESH_AVG_OCCUPANCY_BAD|
      - :math:`<` |THRESH_AVG_OCCUPANCY_BAD|
      - none                       
    * - Atomic Clashes
      - :math:`=` |THRESH_ATOMIC_CLASHES_GOOD|
      - :math:`\leq` |THRESH_ATOMIC_CLASHES_BAD|
      - :math:`>` |THRESH_ATOMIC_CLASHES_BAD| 
      - none                       
    * - Stereochemical Errors
      - :math:`=` |THRESH_STEREO_CHEMICAL_ERRORS_GOOD|
      - :math:`\leq` |THRESH_STEREO_CHEMICAL_ERRORS_BAD|
      - :math:`>` |THRESH_STEREO_CHEMICAL_ERRORS_BAD|
      - :footcite:t:`Shao.2022` 


.. _rmsz_rationale:

|:rotating_light:| RMSZ values
==============================

The Z-value :math:`Z_i` describes the deviation of a value :math:`X_i` from the mean of a population :math:`\mu` in units of standard deviation :math:`\sigma`.
The Z-value is calculated as follows:

.. math::

    Z_i = \frac{X_i - \mu}{\sigma}


Where :math:`N` is the number of values, the :math:`RMSZ` value is calculated as follows:

.. math::

    RMSZ = \sqrt{\frac{1}{N} \sum_{i=1}^{N} Z_i^2}


The Z-values are normalized to have a mean of 0, therefore the RMSZ value is equivalent to the 
*standard deviation* of the Z-values. (see proof below)

.. math::

    s_z = \sqrt{\frac{k}{k-1}} \text{RMSZ}


Applying this to the RMSZ Bond Length and RMSZ Bond Angle 
metrics, we can see that the RMSZ value for the ligand is a measure of the deviation of the bond
lengths/angles from the mean bond lengths/angles of the specific type.
To evaluate the provided RMSZ value, we need to know if the sample size (number of bonds/angles 
in a ligand) is large enough to prove statistical deviation from the underlying general population.
This is done by performing a series of **F-tests for equality of variances**. 

The F-test is a statistical test that compares the variances of two populations. The null hypothesis 
is that the variances are equal, and the alternative hypothesis is that the variances are not equal. 
The F-test returns a p-value, which is the probability of observing the data given that the null 
hypothesis is true. If the p-value is less than a certain threshold (usually 5%), the null 
hypothesis is rejected and we conclude that the variances are not equal.
In this software the number of bonds/angles are scraped from the RCSB API along with the RMSZ values. 
Using a bisectional search algorithm, the treshold-RMSZ value is determined, below which no 
statistical deviation from the reference distribution can be proven.
The size of the true population was set to be 1E7, which is large enough to be considered infinite
compared to the number of bonds/angles typically found in a ligand which we assume to be typically in the 
range of 5-200.

.. figure:: ./_static/rmsz_categorization.png
    :width: 400px


.. figure:: ./_static/rmsz_ok_threshold.png
    :width: 400px


|:abacus:| Mathematical Proof
=============================


To derive :math:`s_z^2 = \frac{k}{k-1} \text{RMSZ}^2`, we start from the fundamental 
definitions of sample variance (:math:`s_z^2`), root mean square (:math:`\text{RMSZ}`), and sample mean (:math:`\bar{z}`).

Definitions
-----------

.. math::

    s_z^2 &= \frac{1}{k-1} \sum_{i=1}^{k} (z_i - \bar{z})^2 \\
    \text{RMSZ} &= \sqrt{\frac{1}{k} \sum_{i=1}^{k} z_i^2} \implies \text{RMSZ}^2 = \frac{1}{k} \sum_{i=1}^{k} z_i^2 \\
    \bar{z} &= \frac{1}{k} \sum_{i=1}^{k} z_i \implies \sum_{i=1}^{k} z_i = k\bar{z}


Step 1: Expand the sum of squared deviations
--------------------------------------------

.. math::

    \sum_{i=1}^{k} (z_i - \bar{z})^2 &= \sum_{i=1}^{k} \left( z_i^2 - 2z_i\bar{z} + \bar{z}^2 \right) \\
    &= \sum_{i=1}^{k} z_i^2 - 2\bar{z}\sum_{i=1}^{k} z_i + \sum_{i=1}^{k} \bar{z}^2


Step 2: Substitute the sample mean terms
----------------------------------------

Substitute :math:`\sum_{i=1}^{k} z_i = k\bar{z}` and :math:`\sum_{i=1}^{k} \bar{z}^2 = k\bar{z}^2`:

.. math::

    \sum_{i=1}^{k} (z_i - \bar{z})^2 &= \sum_{i=1}^{k} z_i^2 - 2\bar{z}(k\bar{z}) + k\bar{z}^2 \\
    &= \sum_{i=1}^{k} z_i^2 - k\bar{z}^2


Step 3: Substitute back into the sample variance definition
-----------------------------------------------------------


.. math::

    s_z^2 = \frac{1}{k-1} \left( \sum_{i=1}^{k} z_i^2 - k\bar{z}^2 \right)


Step 4: Substitute :math:`\text{RMSZ}^2` into the equation
----------------------------------------------------------

Since :math:`\sum_{i=1}^{k} z_i^2 = k \cdot \text{RMSZ}^2`:

.. math::

    s_z^2 &= \frac{1}{k-1} \left( k \cdot \text{RMSZ}^2 - k\bar{z}^2 \right) \\
    s_z^2 &= \frac{k}{k-1} \left( \text{RMSZ}^2 - \bar{z}^2 \right)


Step 5: Apply zero-mean assumption
----------------------------------

Assuming the z-scores are centered at zero (:math:`\bar{z} = 0`):

.. math::

    s_z^2 = \frac{k}{k-1} \text{RMSZ}^2


.. footbibliography::