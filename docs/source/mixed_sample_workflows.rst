==============================
Additional TheiaCoV workflows 
==============================

TheiaCoV_WWVC
===============
The TheiaCoV WasteWater Variant Calling workflow performs variant calling on SARS-CoV-2 in waster water samples and identifies mutations in the Spike gene associated with known VOCs and VUIs. It undertakes read trimming and taxonomic assignment, alignment of reads to a reference genome then a modified version of the `WasteWaterVariantCalling WDL Worfklow produced by Colorado Department of Public Health and Environment (CDPHE) <https://github.com/CDPHE/WasteWaterVariantCalling>`_.

**Input:**
* Reads- Illumina PE
* Primer bed file
* spike bed file
* spike annotations

**Actions:** Runs basic QC (fastq-scan), trimming (trimmomatic), and taxonomic ID (Kraken2) on illumina PE reads, 
aligns reads to optional ref/Wuhan-1 ref with bwa, 
primer seq's are removed with iVar trim, 
WasteWaterVariantCalling wf: performs variant calling and variant filtration using Freebayes, uses bcftools view to pull out VOC/VUI-associated spike mutations

.. toggle-header::
    :header: **References**

        When publishing work using TheiaCoV_WWVC, please reference the following:

        **fastq-scan** Petit RA, III. 2020. fastq-scan. Output FASTQ summary statistics in JSON format. https://github.com/rpetit3/fastq-scan.

        **trimmomatic** Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014 Aug 1;30(15):2114-20.
        
        **Kraken2** Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome biology. 2019 Dec;20(1):1-3.
        
        **BWA** Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. bioinformatics. 2009 Jul 15;25(14):1754-60.

        **iVar** Grubaugh ND, Gangavarapu K, Quick J, Matteson NL, De Jesus JG, Main BJ, Tan AL, Paul LM, Brackney DE, Grewal S, Gurfield N. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome biology. 2019 Dec;20(1):1-9.

        **bcftools** Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb;10(2):giab008.
        
        **freebayes** Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907. 2012 Jul 17.
