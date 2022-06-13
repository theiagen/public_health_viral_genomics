==============================
TheiaCoV Epidemiology Series
==============================

Genomic Epidemiology is an important approach in the effort to mitigate disease transmission. An often-critical step is generating phylogenetic trees to explore the genetic relationship between viruses on a local, national or global scale. The TheiaCoV Epidemiology Series facilitates this for SARS-CoV-2 by generating files for visualization of phylogenetic trees with accompanying metadata.

The TheiaCoV Epidemiology Series contains two WDL workflows, TheiaCoV_Augur_Prep and TheiaCoV_Augur_Run. These must be run sequentially to first prepare each individual sample for running Augur, and second to run Augur itself on the set of samples, generating the phylogenetic tree files with accompanying metadata. The outputs from these workflows can be vizualised in `Auspice <https://docs.nextstrain.org/projects/auspice/en/latest/>`_ and `UShER <https://github.com/yatisht/usher>`_. Helpful resources from the CDC provide an `introduction to NextStrain <https://www.cdc.gov/amd/training/covid-toolkit/module3-1.html>`_ (which includes Auspice), guides to NextStrain `interactive trees <https://www.cdc.gov/amd/training/covid-toolkit/module3-4.html>`_ and an `introduction of UShER <https://www.cdc.gov/amd/training/covid-toolkit/module3-3.html>`_. There is also a video about `how to read trees <https://www.cdc.gov/amd/training/covid-toolkit/module1-3.html>`_ if this is new to you.

TheiaCoV_Augur_Prep
=====================
The TheiaCoV_Augur_Prep workflow was written to prepare individual sample assemblies and their metadata for subsequently running the TheiaCoV_Augur_Run analysis. 

**Input:** The TheiaCoV_Augur_Prep workflow takes assembly FASTA files and metadata formatted in a datatable. FASTA files may be generated with one of the TheiaCoV Characterization workflows and should adhere to quality control guidelines, (e.g. `QC guidelines produced by PH4GE <https://github.com/pha4ge/pipeline-resources/blob/udubs-qc-guidance-dev/docs/qc-solutions.md#gisaid-assembly-acceptance-criteria>`_). The metadata can be uploaded to Terra as TSV file, formatted as in this `example <https://docs.google.com/spreadsheets/d/1PF1u3R-ZGm53UiVsTlIcpg9Qk2dUJgtx/edit#gid=253517867>`_.

**Action:** Using BASH commands, assembly files are de-identified and their metadata are reformatted into individual augur_metadata.tsv files which are later used with the Augur software in TheiaCoV_Augur_Run_.

More details for using the TheiaCoV_Augur_Prep workflow can be found in tables of `required inputs <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_prep_required_inputs.csv>`_, `optional inputs <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_prep_optional_inputs.csv>`_ and `outputs <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_prep_outputs.csv>`_.

TheiaCoV_Augur_Run
====================
The TheiaCoV_Augur_Run workflow takes a set of assembly/concensus files (FASTA format) and sample metadata files (TSV format) that have been reformatted using TheiaCoV_Augur_Prep_ and runs Augur via a a modified version of The Broad Institute's sarscov2_nextstrain WDL workflow to generate the phylogenetic tree files with accompanying metadata. It will also be used to infer SNP distances and create a static PDF report.

**Input**: The TheiaCoV_Augur_Run workflow takes in a *set* of SARS-CoV-2 FASTA and metadata files. If running the workflow via Terra, individual samples will need to be added to a set level datatable before running the workflow. Input FASTAs should meet QA metrics. Sets of FASTAs with highly discordant quality metrics may result in inaccurate inference of genetic relatedness. There must be some sequence diversity among the set of input assemblies. If insufficient diversity is present, it may be necessary to add a more divergent sequence to the set. 

**Action:** The TheiaCoV_Augur_Run workflow uses the inputs to generate a phylogenetic tree in JSON format that is compatible with phylogenetic tree visualization software. 

**Output:** The output JSON is intended to be uploaded to `Auspice <https://clades.nextstrain.org/>`_ to view the phylogenetic tree. This provides a visualization of the genetic relationships between your set of samples. The metadata_merged output can also be uploaded to add context to the phylogenetic visualization. The combined_assemblies output can be uploaded to `UShER <https://genome.ucsc.edu/cgi-bin/hgPhyloPlace>`_ to view the samples on a global tree of representative sequences from the public repositories.

.. note::
   You may generate phylogenies multiple times, running the TheiaCov_Augur_Run workflow, assessing results and amending inputs to generate a final tree with suffient diversity and high-quality data of interest.

More details for using the TheiaCoV_Augur_Run workflow can be found in tables of `required inputs <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_run_required_inputs.csv>`_, `optional inputs <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_run_optional_inputs.csv>`_ and `outputs <https://github.com/theiagen/public_health_viral_genomics/blob/main/docs/source/tables/theiacov_workflows/theiacov_augur_run_outputs.csv>`_.

.. toggle-header::
    :header: **References**

        When publishing work using TheiaCoV_Augur workflows, please reference the following:
        
         **Nextstrain:** Hadfield J, Megill C, Bell SM, Huddleston J, Potter B, Callender C, Sagulenko P, Bedford T, Neher RA. Nextstrain: real-time tracking of pathogen evolution. Bioinformatics. 2018 Dec 1;34(23):4121-3.
      
        To publish work using inferences from UShER, please references
        
         **UShER:** Turakhia Y, Thornlow B, Hinrichs AS, De Maio N, Gozashti L, Lanfear R, Haussler D, Corbett-Detig R. Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic. Nature Genetics. 2021 Jun;53(6):809-16.