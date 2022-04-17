version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_qc_utils.wdl" as qc_utils
import "../tasks/task_versioning.wdl" as versioning

workflow hivgc_ont {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on ONT NGS platforms."
  }
  input {
    String samplename
    String seq_method = "OXFORD_NANOPORE"
    File primer_bed
    File reference_genome
    File demultiplexed_reads
    File mutation_db
    Int normalise = 200
    Int max_length = 600
    Int min_length = 50
  }
  call qc_utils.fastq_scan_se as fastq_scan_raw_reads {
    input:
      read1 = demultiplexed_reads
  }
  call medaka.read_filtering as read_filtering {
    input:
      demultiplexed_reads = demultiplexed_reads,
      samplename = samplename,
      min_length = min_length,
      max_length = max_length
  }
  call qc_utils.fastq_scan_se as fastq_scan_clean_reads {
    input:
      read1 = read_filtering.filtered_reads
  }
  call medaka.consensus as consensus {
    input:
      samplename = samplename,
      filtered_reads = read_filtering.filtered_reads,
      primer_bed = primer_bed,
      reference_genome = reference_genome,
      normalise = normalise
  }
  call taxon_ID.quasitools_one_sample as quasitools {
    input:
      samplename = samplename,
      read1 = read_filtering.filtered_reads,
      mutation_db = mutation_db
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
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String hivgc_ont_version = version_capture.phvg_version
    String hivgc_ont_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC
    Int num_reads_raw = fastq_scan_raw_reads.read1_seq
    Int num_reads_clean = fastq_scan_clean_reads.read1_seq
    String fastq_scan_version = fastq_scan_clean_reads.version
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
    # Float s_gene_mean_coverage = stats_n_coverage_primtrim.s_gene_depth
    String samtools_version = stats_n_coverage.samtools_version
    # Lineage Assignment
    ## Put Quasitools stuff here
    String quasitools_version = quasitools.quasitools_version
    String quasitools_date = quasitools.quasitools_date
    File quasitools_coverage_file = quasitools.coverage_file
    File quasitools_dr_report = quasitools.dr_report
    File quasitools_hydra_vcf = quasitools.hydra_vcf
    File quasitools_mutations_report = quasitools.mutations_report
    # Clade Assigment
    ## Put HIVMMER stuff here
  }
}
