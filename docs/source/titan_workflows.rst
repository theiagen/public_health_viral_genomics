======================
Titan Workflow Series
======================

The Titan series is a collection of WDL workflows developed for the genomic characterization and genomic epidemiology of viral samples to support public health decision-making. As of today, May 04th, 2021, these workflows are specific to SARS-CoV-2 analysis, but development is underway to allow for other viral pathogens of concern.


Titan Workflows for Genomic Characterization
--------------------------------------------
Genomic characterization, *i.e.* generating consensus assemblies (FASTA format) from next-generation sequencing (NGS) read data (FASTQ format) to assign samples with relevant nomenclature designation (e.g. PANGO lineage and NextClade clades) is an increasingly critical function to public health laboratories around the world.

The Titan Series includes four separate WDL workflows that process NGS read data from four different sequencing approaches (Illumina paired-end, Illumina single-end, Oxford Nanopore Technology (ONT), and Clear Labs) to generate consensus assemblies, produce relevant quality-control metrics for both the input, read data and the generated assembly, and assign samples with a lineage and clade designation using Pangolin and NextClade, respectively.

All four Titan workflows for genomic characterization will generate a SARS-CoV-2 assembly by mapping input-read data to the Wu Han-1 reference genome, removing primer reads from that alignment, and then calling the consensus assembly based on the primer-trimmed alignment. These consensus assemblies are then fed into the Pangolin and NextClade CLI tools for lineage and clade assignments.

The major difference between each of these Titan workflows is in how the read mapping, primer trimming, and consensus genome calling is performed. More information on the technical details of these processes and information on how to utilize and apply these workflows for public health investigations is available below.


A series of introductory training videos that provide conceptual overviews of methodologies and walkthrough tutorials on how to utilize the Titan workflows through Terra are available on the Theiagen Genomics YouTube page:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/watch?v=zP9I1r6TNrw&list=PLU47xRg_MKJrtyoFwqGiywl7lQj6vq8Uz&index=2" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
   </div>




Titan_Illumina_PE
=================


Titan_Illumina_SE
=================

Titan_ClearLabs
=================

Titan_ONT
=========
