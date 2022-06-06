==============================
Additional TheiaCoV workflows 
==============================

TheiaCoV_WWVC
===============
The TheiaCoV WasteWater Variant Calling workflow performs variant calling on SARS-CoV-2 in waster water samples and identifies mutations in the Spike gene associated with known VOCs and VUIs. It is a modified version of the `WasteWaterVariantCalling WDL Worfklow produced by Colorado Department of Public Health and Environment (CDPHE) <https://github.com/CDPHE/WasteWaterVariantCalling>`_.

**Input:**
* Reads- Illumina PE
* Primer bed file
* spike bed file
* spike annotations

**Actions:** Runs basic QC (fastq-scan), trimming (trimmomatic), and taxonomic ID (Kraken2) on illumina PE reads, aligns with bwa, consensus_call.primer_trim, WasteWaterVariantCalling wf

**Outputs:**
