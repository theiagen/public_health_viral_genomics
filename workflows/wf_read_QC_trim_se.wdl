version 1.0

import "../tasks/task_qc_utils.wdl" as qc_utils
import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_taxonID.wdl" as taxonID

workflow read_QC_trim {
  meta {
    description: "Runs basic QC (fastq-scan), trimming (SeqyClean), and taxonomic ID (Kraken2) on illumina PE reads"
  }

  input {
    String samplename
    File read1_raw
    Int? trimmomatic_minlen = 25
    Int? trimmomatic_quality_trim_score = 30
    Int? trimmomatic_window_size = 4
    Int  bbduk_mem = 8
  }
#  call read_clean.ncbi_scrub_se {
#    input:
#      samplename = samplename,
#      read1 = read1_raw
#  }
  call read_clean.trimmomatic_se {
    input:
      samplename = samplename,
      read1 = read1_raw,
      trimmomatic_minlen = trimmomatic_minlen,
      trimmomatic_quality_trim_score = trimmomatic_quality_trim_score,
      trimmomatic_window_size = trimmomatic_window_size
  }
  call read_clean.bbduk_se {
    input:
      samplename = samplename,
      read1_trimmed = trimmomatic_se.read1_trimmed,
      memory = bbduk_mem
  }
  call qc_utils.fastq_scan_se as fastq_scan_raw {
    input:
      read1 = read1_raw
  }
  call qc_utils.fastq_scan_se as fastq_scan_clean {
    input:
      read1 = bbduk_se.read1_clean
  }
  call taxonID.kraken2 as kraken2_raw {
    input:
      samplename = samplename,
      read1 = bbduk_se.read1_clean
  }
#  call taxonID.kraken2 as kraken2_dehosted {
#    input:
#      samplename = samplename,
#      read1 = ncbi_scrub_se.read1_dehosted
#  }

  output {
    File    read1_clean               = bbduk_se.read1_clean

    Int     fastq_scan_number_reads       = fastq_scan_raw.read1_seq
    Int     fastq_scan_clean_number_reads = fastq_scan_clean.read1_seq

    String  kraken_version            = kraken2_raw.version
    Float   kraken_human              = kraken2_raw.percent_human
    Float   kraken_sc2                = kraken2_raw.percent_sc2
    String  kraken_report             = kraken2_raw.kraken_report
#    Float    kraken_human_dehosted    =    kraken2_dehosted.percent_human
#    Float    kraken_sc2_dehosted    =    kraken2_dehosted.percent_sc2
#    String    kraken_report_dehosted    =    kraken2_dehosted.kraken_report

    String  fastq_scan_version          = fastq_scan_raw.version
    String  bbduk_docker            = bbduk_se.bbduk_docker
    String  trimmomatic_version     = trimmomatic_se.version
  }
}
