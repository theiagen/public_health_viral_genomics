version 1.0

import "../tasks/task_taxonID.wdl" as taxon_id
import "wf_read_QC_trim.wdl" as read_qc
import "../tasks/task_alignment.wdl" as align
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_versioning.wdl" as versioning

workflow freyja_fastq {
  input {
    File read1_raw
    File read2_raw
    File primer_bed
    File reference_genome
    Int trimmomatic_minlen = 25
    String samplename
    Int min_depth = 10
  }
  call read_qc.read_QC_trim {
    input:
      samplename = samplename,
      read1_raw  = read1_raw,
      read2_raw  = read2_raw,
      trimmomatic_minlen = trimmomatic_minlen
  }
  call align.bwa {
    input:
      samplename = samplename,
      reference_genome=reference_genome,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean
  }
  call consensus_call.primer_trim {
    input:
      samplename = samplename,
      primer_bed = primer_bed,
      bamfile = bwa.sorted_bam
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam,
      min_depth = min_depth
  }
  call taxon_id.freyja_one_sample as freyja {
    input:
      primer_trimmed_bam = primer_trim.trim_sorted_bam,
      samplename = samplename,
      reference_genome = reference_genome
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String freyja_fastq_wf_version = version_capture.phvg_version
    String freyja_fastq_wf_analysis_date = version_capture.date
    # Raw Read QC
    File read1_dehosted = read_QC_trim.read1_dehosted
    File read2_dehosted = read_QC_trim.read2_dehosted
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    Int num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    # Read Trim
    Int num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String trimmomatic_version = read_QC_trim.trimmomatic_version
    String bbduk_docker = read_QC_trim.bbduk_docker
    # Contaminent Check
    String kraken_version = read_QC_trim.kraken_version
    Float kraken_human = read_QC_trim.kraken_human
    Float kraken_sc2 = read_QC_trim.kraken_sc2
    String kraken_report = read_QC_trim.kraken_report
    Float kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Mapping and Alignment
    String bwa_version = bwa.bwa_version
    String samtools_version = bwa.sam_version
    String alignment_method = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
    File aligned_bam = primer_trim.trim_sorted_bam
    File aligned_bai = primer_trim.trim_sorted_bai
    Float primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String ivar_version_primtrim = primer_trim.ivar_version
    String samtools_version_primtrim = primer_trim.samtools_version
    String primer_bed_name = primer_trim.primer_bed_name
    # Alignment QC
    File alignment_stats = stats_n_coverage.stats
    File alignment_flagstat = stats_n_coverage.flagstat
    Float meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    Float alignment_mean_coverage = stats_n_coverage_primtrim.depth
    Float s_gene_mean_coverage = stats_n_coverage_primtrim.s_gene_depth
    Float s_gene_percent_coverage = stats_n_coverage_primtrim.s_gene_percent_coverage
    File percent_gene_coverage = stats_n_coverage_primtrim.percent_gene_coverage
    String samtools_version_stats = stats_n_coverage.samtools_version
    # Freyja Analysis
    File freyja_variants = freyja.freyja_variants
    File freyja_depths = freyja.freyja_depths
    File freyja_demixed = freyja.freyja_demixed
    String freyja_barcode_version = freyja.freyja_barcode_version
    String freyja_metadata_version = freyja.freyja_metadata_version
    }
}
