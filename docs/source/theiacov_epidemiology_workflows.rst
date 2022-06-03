
TheiaCoV_Augur_Prep
================
The TheiaCoV_Augur_Prep workflow was written to process consensus assemblies (FASTA format) and the associated metadata in preparation for running the TheiaCoV_Augur_Run. Input assemblies should be of similar quality (percent reference coverage, number of ambiguous bases, etc.). Inputs with highly discordant quality metrics may result in inaccurate inference of genetic relatedness.

.. note::
  There must be some sequence diversity in the input set of assemblies to be analyzed. As a rule of thumb, the smaller the input set, the more sequence diversity will be required to make any sort of genomic inference. If a small (~10) set of viral genomic assemblies is used as the input then it may be necessary to add one significantly divergent assembly.

Upon initiating a TheiaCoV_Augur_Prep run, input assembly/consensus files and associated metadata will be used to produce the array of assembly/consensus files and the array of metadata files to be used as inputs for the TheiaCoV_Augur_Run workflow.

Metadata files are prepared with the Augur_Prep workflow by using BASH commands to first de-identify, and then to parse the headers of the input assembly files.

Required User Inputs
********************
Download CSV: :download:`TheiaCoV_Augur_Prep_required_inputs.csv <tables/theiacov_workflows/theiacov_augur_prep_required_inputs.csv>`

.. csv-table::
   :file: tables/theiacov_workflows/theiacov_augur_prep_required_inputs.csv
   :widths: 20, 20, 20, 40
   :header-rows: 1
|

TheiaCoV_Augur_Run
===============
The TheiaCoV_Augur_Run workflow was written to process an array of assembly/consensus files (FASTA format) and and array of sample metadata files (TSV format) using a modified version of The Broad Institute's sarscov2_nextstrain WDL workflow to create an Auspice JSON file; output from the modified sarscov2_nextstrain workflow will also be used to infer SNP distances and create a static PDF report.

Upon initiating a TheiaCoV_Augur_Run run, the input assembly/consensus file array and the associated metadata file array will be used to generate a JSON file that is compatible with phylogenetic tree building software. This JSON can then be used in Auspice or Nextstrain to view the phylogenetic tree. This phylogeneic tree can be used in genomic epidemiological analysis to visualize the genetic relatedness of a set of samples. The associated metadata can then be used to add context to the phylogenetic visualization.

Required User Inputs
********************
Download CSV: :download:`TheiaCoV_Augur_Run_required_inputs.csv <tables/theiacov_workflows/theiacov_augur_prep_required_inputs.csv>`

.. csv-table::
   :file: tables/theiacov_workflows/theiacov_augur_run_required_inputs.csv
   :widths: 20, 20, 20, 40
   :header-rows: 1

|
