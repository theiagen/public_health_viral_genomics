version 1.0

import "wf_titan_clearlabs.wdl" as assembly

workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], File]] inputSamples
  }
  scatter (sample in inputSamples) {
    call assembly.titan_clearlabs {
      input:
        samplename = sample.left[0],
        clear_lab_fastq  = sample.right
    }
  }
  output {
    Array[String]	seq_platform	=	titan_clearlabs.seq_platform

    Array[File]	dehosted_reads	=	titan_clearlabs.dehosted_reads

    Array[String]	kraken_version	=	titan_clearlabs.kraken_version
    Array[Float]	kraken_human	=	titan_clearlabs.kraken_human
    Array[Float]	kraken_sc2	=	titan_clearlabs.kraken_sc2
    Array[String]	kraken_report	=	titan_clearlabs.kraken_report
    Array[Float]	kraken_human_dehosted	=	titan_clearlabs.kraken_human_dehosted
    Array[Float]	kraken_sc2_dehosted	=	titan_clearlabs.kraken_sc2_dehosted
    Array[String]	kraken_report_dehosted	=	titan_clearlabs.kraken_report_dehosted

    Array[File]	aligned_bam	=	titan_clearlabs.aligned_bam
    Array[File]	aligned_bai	=	titan_clearlabs.aligned_bai
    Array[File]	variants_from_ref_vcf	=	titan_clearlabs.variants_from_ref_vcf
    Array[String]	artic_version	=	titan_clearlabs.artic_version
    Array[File]	assembly_fasta	=	titan_clearlabs.assembly_fasta
    Array[Int]	number_N	=	titan_clearlabs.number_N
    Array[Int]	assembly_length_unambiguous	=	titan_clearlabs.assembly_length_unambiguous
    Array[Int]	number_Degenerate	=	titan_clearlabs.number_Degenerate
    Array[Int]	number_Total	=	titan_clearlabs.number_Total
    Array[Float]	pool1_percent	=	titan_clearlabs.pool1_percent
    Array[Float]	pool2_percent	=	titan_clearlabs.pool2_percent
    Array[Float]	percent_reference_coverage	=	titan_clearlabs.percent_reference_coverage
    Array[String]	assembly_method	=	titan_clearlabs.assembly_method

    Array[File]	consensus_stats	=	titan_clearlabs.consensus_stats
    Array[File]	consensus_flagstat	=	titan_clearlabs.consensus_flagstat
    Array[Float]	meanbaseq_trim	=	titan_clearlabs.meanbaseq_trim
    Array[Float]	meanmapq_trim	=	titan_clearlabs.meanmapq_trim
    Array[Float]	assembly_mean_coverage	=	titan_clearlabs.assembly_mean_coverage
    Array[String]	samtools_version	=	titan_clearlabs.samtools_version

    Array[String]	pango_lineage	=	titan_clearlabs.pango_lineage
    Array[Float]	pangolin_aLRT	=	titan_clearlabs.pangolin_aLRT
    Array[String]	pangolin_version	=	titan_clearlabs.pangolin_version
    Array[File]	pango_lineage_report	=	titan_clearlabs.pango_lineage_report
    Array[String]	pangolin_docker	=	titan_clearlabs.pangolin_docker

    Array[File]	nextclade_json	=	titan_clearlabs.nextclade_json
    Array[File]	auspice_json	=	titan_clearlabs.auspice_json
    Array[File]	nextclade_tsv	=	titan_clearlabs.nextclade_tsv
    Array[String]	nextclade_clade	=	titan_clearlabs.nextclade_clade
    Array[String]	nextclade_aa_subs	=	titan_clearlabs.nextclade_aa_subs
    Array[String]	nextclade_aa_dels	=	titan_clearlabs.nextclade_aa_dels
    Array[String]	nextclade_version	=	titan_clearlabs.nextclade_version

    Array[File]	vadr_alerts_list	=	titan_clearlabs.vadr_alerts_list
    Array[Int]	vadr_num_alerts	=	titan_clearlabs.vadr_num_alerts
    Array[String]	vadr_docker	=	titan_clearlabs.vadr_docker
  }
}
