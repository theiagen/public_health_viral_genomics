version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_amplicon_metrics.wdl" as amplicon_metrics
import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_qc_utils.wdl" as qc_utils

workflow titan_clearlabs {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }

  input {
    String  samplename
    File  clear_lab_fastq
    String  seq_method  = "ONT via Clear Labs WGS"
    String? artic_primer_version  = "V3"
    String  pangolin_docker_image = "staphb/pangolin:3.0.5-pangolearn-2021-06-05"
    Int?  normalise  = 20000
  }
  call qc_utils.fastqc_se as fastqc_se_raw {
    input:
      read1 = clear_lab_fastq
  }
  call read_clean.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = clear_lab_fastq
  }
  call qc_utils.fastqc_se as fastqc_se_clean {
    input:
      read1 = ncbi_scrub_se.read1_dehosted
  }
  call taxon_ID.kraken2 as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted
  }
  call medaka.consensus {
    input:
      samplename = samplename,
      filtered_reads = ncbi_scrub_se.read1_dehosted,
      artic_primer_version = artic_primer_version,
      normalise = normalise
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = consensus.sorted_bam
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = consensus.trim_sorted_bam
  }
  call taxon_ID.pangolin3 {
    input:
      samplename = samplename,
      fasta = consensus.consensus_seq,
      docker = pangolin_docker_image
  }
  call taxon_ID.kraken2 as kraken2_raw {
    input:
      samplename = samplename,
      read1 = clear_lab_fastq
  }
  call taxon_ID.nextclade_one_sample {
    input:
      genome_fasta = consensus.consensus_seq
  }
  call amplicon_metrics.bedtools_cov {
    input:
      bamfile = consensus.trim_sorted_bam,
      baifile = consensus.trim_sorted_bai
  }
  call ncbi.vadr {
    input:
      genome_fasta = consensus.consensus_seq,
  }
  output {
    String titan_gc_version = "Titan_ClearLabs v1.4.4"
    String	seq_platform	=	seq_method

    File	dehosted_reads	=	ncbi_scrub_se.read1_dehosted

    Int fastqc_raw = fastqc_se_raw.number_reads
    Int fastqc_clean = fastqc_se_clean.number_reads
    
  	String	kraken_version	=	kraken2_raw.version
  	Float	kraken_human	=	kraken2_raw.percent_human
  	Float	kraken_sc2	=	kraken2_raw.percent_sc2
  	String	kraken_report	=	kraken2_raw.kraken_report
    Float	kraken_human_dehosted	=	kraken2_dehosted.percent_human
  	Float	kraken_sc2_dehosted	=	kraken2_dehosted.percent_sc2
  	String	kraken_report_dehosted	=	kraken2_dehosted.kraken_report

  	File	aligned_bam	=	consensus.trim_sorted_bam
  	File	aligned_bai	=	consensus.trim_sorted_bai
  	File	variants_from_ref_vcf	=	consensus.medaka_pass_vcf
  	String	artic_version	=	consensus.artic_pipeline_version
  	File	assembly_fasta	=	consensus.consensus_seq
  	Int	number_N	=	consensus.number_N
  	Int	assembly_length_unambiguous	=	consensus.number_ATCG
  	Int	number_Degenerate	=	consensus.number_Degenerate
  	Int	number_Total	=	consensus.number_Total
  	Float	pool1_percent	=	consensus.pool1_percent
  	Float	pool2_percent	=	consensus.pool2_percent
  	Float	percent_reference_coverage	=	consensus.percent_reference_coverage
  	String	assembly_method	=	consensus.artic_pipeline_version

  	File	consensus_stats	=	stats_n_coverage.stats
  	File	consensus_flagstat	=	stats_n_coverage.flagstat
  	Float	meanbaseq_trim	=	stats_n_coverage_primtrim.meanbaseq
  	Float	meanmapq_trim	=	stats_n_coverage_primtrim.meanmapq
  	Float	assembly_mean_coverage	=	stats_n_coverage_primtrim.depth
  	String	samtools_version	=	stats_n_coverage.samtools_version

    String	pango_lineage	=	pangolin3.pangolin_lineage
    String	pangolin_conflicts	=	pangolin3.pangolin_conflicts
    String pangolin_notes = pangolin3.pangolin_notes
  	String	pangolin_version	=	pangolin3.version
  	File	pango_lineage_report	=	pangolin3.pango_lineage_report
  	String	pangolin_docker	=	pangolin3.pangolin_docker
  	String	pangolin_usher_version	=	pangolin3.pangolin_usher_version

  	File	nextclade_json	=	nextclade_one_sample.nextclade_json
  	File	auspice_json	=	nextclade_one_sample.auspice_json
  	File	nextclade_tsv	=	nextclade_one_sample.nextclade_tsv
  	String	nextclade_clade	=	nextclade_one_sample.nextclade_clade
  	String	nextclade_aa_subs	=	nextclade_one_sample.nextclade_aa_subs
  	String	nextclade_aa_dels	=	nextclade_one_sample.nextclade_aa_dels
  	String	nextclade_version	=	nextclade_one_sample.nextclade_version

  	File	vadr_alerts_list	=	vadr.alerts_list
  	Int	vadr_num_alerts	=	vadr.num_alerts
  	String	vadr_docker	=	vadr.vadr_docker
  }
}
