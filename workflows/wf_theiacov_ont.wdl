version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_qc_utils.wdl" as qc_utils
import "../tasks/task_versioning.wdl" as versioning

workflow theiacov_ont {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on ONT NGS platforms."
  }
  input {
    String samplename
    String seq_method = "OXFORD_NANOPORE"
    File primer_bed
    File demultiplexed_reads
    Int? normalise = 200
    String nextclade_dataset_name = "sars-cov-2"
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2022-01-18T12:00:00Z"
    Int? max_length = 700
    Int? min_length = 400
  }
  call qc_utils.fastq_scan_se as fastq_scan_raw_reads {
    input:
      read1 = demultiplexed_reads
  }
  call read_clean.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = demultiplexed_reads
  }
  call medaka.read_filtering {
    input:
      demultiplexed_reads = ncbi_scrub_se.read1_dehosted,
      samplename = samplename,
      min_length = min_length,
      max_length = max_length
  }
  call qc_utils.fastq_scan_se as fastq_scan_clean_reads {
    input:
      read1 = read_filtering.filtered_reads
  }
  call taxon_ID.kraken2 as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted
  }
  call medaka.consensus {
    input:
      samplename = samplename,
      filtered_reads = read_filtering.filtered_reads,
      primer_bed = primer_bed,
      normalise = normalise
  }
  call qc_utils.consensus_qc {
    input:
      assembly_fasta = consensus.consensus_seq
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
      fasta = consensus.consensus_seq
  }
  call taxon_ID.kraken2 as kraken2_raw {
    input:
      samplename = samplename,
      read1 = demultiplexed_reads
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
    String theiacov_ont_version = version_capture.phvg_version
    String theiacov_ont_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC
    File reads_dehosted = ncbi_scrub_se.read1_dehosted
    Int num_reads_raw = fastq_scan_raw_reads.read1_seq
    Int num_reads_clean = fastq_scan_clean_reads.read1_seq
    String fastq_scan_version = fastq_scan_clean_reads.version
    String kraken_version = kraken2_raw.version
    Float kraken_human = kraken2_raw.percent_human
    Float kraken_sc2 = kraken2_raw.percent_sc2
    String kraken_report = kraken2_raw.kraken_report
    Float kraken_human_dehosted = kraken2_dehosted.percent_human
    Float kraken_sc2_dehosted = kraken2_dehosted.percent_sc2
    String kraken_report_dehosted = kraken2_dehosted.kraken_report
    # Read Alignment
    File aligned_bam = consensus.trim_sorted_bam
    File aligned_bai = consensus.trim_sorted_bai
    File variants_from_ref_vcf = consensus.medaka_pass_vcf
    String artic_version = consensus.artic_pipeline_version
    String artic_docker = consensus.artic_pipeline_docker
    String medaka_reference = consensus.medaka_reference
    String primer_bed_name = consensus.primer_bed_name
    File assembly_fasta = consensus.consensus_seq
    String assembly_method = consensus.artic_pipeline_version
    # Assembly QC
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
    String samtools_version = stats_n_coverage.samtools_version
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
    String nextclade_aa_subs = nextclade_output_parser_one_sample.nextclade_aa_subs
    String nextclade_aa_dels = nextclade_output_parser_one_sample.nextclade_aa_dels
    String nextclade_clade = nextclade_output_parser_one_sample.nextclade_clade
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String vadr_num_alerts = vadr.num_alerts
    String vadr_docker = vadr.vadr_docker
  }
}
