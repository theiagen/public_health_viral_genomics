version 1.0

import "wf_read_QC_trim.wdl" as read_qc
import "../tasks/task_alignment.wdl" as align
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_versioning.wdl" as versioning
import "../tasks/task_qc_utils.wdl" as qc_utils

workflow theiacov_illumina_pe {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File read2_raw
    File primer_bed
    String nextclade_dataset_name = "sars-cov-2"
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2022-01-18T12:00:00Z"
  }
  call read_qc.read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = read1_raw,
      read2_raw = read2_raw
  }
  call align.bwa {
    input:
      samplename = samplename,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean
  }
  call consensus_call.primer_trim {
    input:
      samplename = samplename,
      primer_bed = primer_bed,
      bamfile = bwa.sorted_bam
  }
  call consensus_call.variant_call {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam
  }
  call consensus_call.consensus {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam
  }
  call qc_utils.consensus_qc {
    input:
      assembly_fasta = consensus.consensus_seq
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = bwa.sorted_bam
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam
  }
  call taxon_ID.pangolin3 {
    input:
      samplename = samplename,
      fasta = consensus.consensus_seq
  }
  call taxon_ID.nextclade_one_sample {
    input:
      genome_fasta = consensus.consensus_seq,
      dataset_name = nextclade_dataset_name,
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
  }
  call taxon_ID.nextclade_output_parser_one_sample {
    input:
      nextclade_tsv = nextclade_one_sample.nextclade_tsv
  }
  call ncbi.vadr {
    input:
      genome_fasta = consensus.consensus_seq,
      assembly_length_unambiguous = consensus_qc.number_ATCG
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_illumina_pe_version = version_capture.phvg_version
    String theiacov_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String  seq_platform = seq_method
    # Read QC
    File read1_dehosted = read_QC_trim.read1_dehosted
    File read2_dehosted = read_QC_trim.read2_dehosted
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    Int fastq_scan_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    Int fastq_scan_clean1 = read_QC_trim.fastq_scan_clean1
    Int fastq_scan_clean2 = read_QC_trim.fastq_scan_clean2
    String fastq_scan_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String trimmomatic_version = read_QC_trim.trimmomatic_version
    String bbduk_docker = read_QC_trim.bbduk_docker
    String kraken_version = read_QC_trim.kraken_version
    Float kraken_human = read_QC_trim.kraken_human
    Float kraken_sc2 = read_QC_trim.kraken_sc2
    String kraken_report = read_QC_trim.kraken_report
    Float kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment
    String bwa_version = bwa.bwa_version
    String samtools_version = bwa.sam_version
    String assembly_method = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
    File aligned_bam = primer_trim.trim_sorted_bam
    File aligned_bai = primer_trim.trim_sorted_bai
    Float primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String ivar_version_primtrim = primer_trim.ivar_version
    String samtools_version_primtrim   = primer_trim.samtools_version
    String primer_bed_name = primer_trim.primer_bed_name
    File ivar_tsv = variant_call.sample_variants_tsv
    File ivar_vcf = variant_call.sample_variants_vcf
    String ivar_variant_version = variant_call.ivar_version
    # Assembly QC
    File assembly_fasta = consensus.consensus_seq
    String ivar_version_consensus = consensus.ivar_version
    String samtools_version_consensus = consensus.samtools_version
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Alignment QC
    File consensus_stats = stats_n_coverage.stats
    File consensus_flagstat = stats_n_coverage.flagstat
    Float meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    Float assembly_mean_coverage = stats_n_coverage_primtrim.depth
    Float s_gene_mean_coverage = stats_n_coverage_primtrim.s_gene_depth
    String samtools_version_stats = stats_n_coverage.samtools_version
    # Lineage Assignment
    String pango_lineage = pangolin3.pangolin_lineage
    String pangolin_conflicts = pangolin3.pangolin_conflicts
    String pangolin_notes = pangolin3.pangolin_notes
    String pangolin_assignment_version = pangolin3.pangolin_assignment_version
    File pango_lineage_report = pangolin3.pango_lineage_report
    String pangolin_docker = pangolin3.pangolin_docker
    String pangolin_versions = pangolin3.pangolin_versions
    # Clade Assigment
    File nextclade_json = nextclade_one_sample.nextclade_json
    File auspice_json = nextclade_one_sample.auspice_json
    File nextclade_tsv = nextclade_one_sample.nextclade_tsv
    String nextclade_version = nextclade_one_sample.nextclade_version
    String nextclade_docker = nextclade_one_sample.nextclade_docker
    String  nextclade_aa_subs = nextclade_output_parser_one_sample.nextclade_aa_subs
    String  nextclade_aa_dels = nextclade_output_parser_one_sample.nextclade_aa_dels
    String  nextclade_clade = nextclade_output_parser_one_sample.nextclade_clade
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String vadr_num_alerts = vadr.num_alerts
    String vadr_docker = vadr.vadr_docker
  }
}
