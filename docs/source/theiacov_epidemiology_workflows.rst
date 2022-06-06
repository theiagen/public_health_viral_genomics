==============================
TheiaCoV Epidemiology Series
==============================

Genomic Epidemiology has been an important approach in the effort to mitigate disease transmission. An often-critical step is generating phylogenetic trees to explore the genetic relationship between viruses on a local, national or global scale. The TheiaCoV Epidemiology Series facilitates this for SARS-CoV-2 by generating files for visualization of phylogenetic trees with accompanying metadata.

The TheiaCoV Epidemiology Series contains two WDL workflows, TheiaCoV_Augur_Prep and TheiaCoV_Augur_Run. These must be run sequentially to first prepare each individual sample for running Augur, and second to run Augur itself on the set of samples, generating the phylogenetic tree files with accompanying metadata. The outputs from these workflows can be vizualised in `Auspice <https://docs.nextstrain.org/projects/auspice/en/latest/>`_ and `UShER <https://github.com/yatisht/usher>`_.

TheiaCoV_Augur_Prep
================
The TheiaCoV_Augur_Prep workflow was written to prepare individual sample assemblies and their metadata for subsequently running the TheiaCoV_Augur_Run analysis. 

**Input:** The TheiaCoV_Augur_Prep workflow takes assembly FASTA files and metadata formatted in a datatable. FASTA files may be generated with one of the TheiaCoV Characterization workflows and should adhere to quality control guidelines, (e.g. `QC guidelines produced by PH4GE <https://github.com/pha4ge/pipeline-resources/blob/udubs-qc-guidance-dev/docs/qc-solutions.md#gisaid-assembly-acceptance-criteria>`_). The metadata can be uploaded to Terra as TSV file, formatted as in this `example <https://docs.google.com/spreadsheets/d/1PF1u3R-ZGm53UiVsTlIcpg9Qk2dUJgtx/edit#gid=253517867>`_.

**Action:** Using BASH commands, assembly files are de-identified and their metadata are reformatted into individual augur_metadata.tsv files for use with the Augur software.

More details for using the TheiaCoV_Augur_Prep workflow can be found in tables of `input requirements <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_prep_required_inputs.csv>`_.

TheiaCoV_Augur_Run
===============
The TheiaCoV_Augur_Run workflow was written to process an array of assembly/consensus files (FASTA format) and and array of sample metadata files (TSV format) using a modified version of The Broad Institute's sarscov2_nextstrain WDL workflow to create an Auspice JSON file; output from the modified sarscov2_nextstrain workflow will also be used to infer SNP distances and create a static PDF report.

Upon initiating a TheiaCoV_Augur_Run run, the input assembly/consensus file array and the associated metadata file array will be used to generate a JSON file that is compatible with phylogenetic tree building software. This JSON can then be used in Auspice or Nextstrain to view the phylogenetic tree. This phylogenetic tree can be used in genomic epidemiological analysis to visualize the genetic relatedness of a set of samples. The associated metadata can then be used to add context to the phylogenetic visualization.

.. note::
   Input FASTAs with highly discordant quality metrics may result in inaccurate inference of genetic relatedness.
  There must be some sequence diversity among the set of input assemblies to be analyzed. If insufficient diversity is present, it may be necessary to add a more divergent sequence to the set.

  You may generate phylogenies multiple times, running the workflows, assessing results and amending inputs to generate a final tree with suffient diversity and high-quality data.
Required User Inputs
********************
Download CSV: :download:`TheiaCoV_Augur_Run_required_inputs.csv <tables/theiacov_workflows/theiacov_augur_prep_required_inputs.csv>`

.. csv-table::
   :file: tables/theiacov_workflows/theiacov_augur_run_required_inputs.csv
   :widths: 20, 20, 20, 40
   :header-rows: 1

|
