version 1.0

import "../tasks/task_qc_utils.wdl" as qc_utils
import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_taxonID.wdl" as taxonID

workflow read_QC_trim {
  meta {
    description: "Runs basic QC (FastQC), trimming (SeqyClean), and taxonomic ID (Kraken2) on illumina PE reads"
  }

  input {
    String  samplename
    File    read1_raw
    File    read2_raw
    Int?    trimmomatic_minlen = 75
    Int?    trimmomatic_quality_trim_score = 30
    Int?    trimmomatic_window_size = 4
  }
  call read_clean.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw
  }  
  call read_clean.trimmomatic {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_pe.read1_dehosted,
      read2 = ncbi_scrub_pe.read2_dehosted,
      trimmomatic_minlen = trimmomatic_minlen,
      trimmomatic_quality_trim_score = trimmomatic_quality_trim_score,
      trimmomatic_window_size = trimmomatic_window_size
  }
  call read_clean.bbduk {
    input:
      samplename = samplename,
      read1_trimmed = trimmomatic.read1_trimmed,
      read2_trimmed = trimmomatic.read2_trimmed
  }
  call qc_utils.fastqc as fastqc_raw {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
  }
  call qc_utils.fastqc as fastqc_clean {
    input:
      read1 = bbduk.read1_clean,
      read2 = bbduk.read2_clean
  }
  call taxonID.kraken2 as kraken2_raw {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw
  }
  call taxonID.kraken2 as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_pe.read1_dehosted,
      read2 = ncbi_scrub_pe.read2_dehosted
  }
  output {
    File   read1_dehosted            = ncbi_scrub_pe.read1_dehosted
    File   read2_dehosted            = ncbi_scrub_pe.read2_dehosted
    Int    read1_human_spots_removed = ncbi_scrub_pe.read1_human_spots_removed
    Int    read2_human_spots_removed = ncbi_scrub_pe.read2_human_spots_removed
    
    File   read1_clean               = bbduk.read1_clean
    File   read2_clean               = bbduk.read2_clean

    Int    fastqc_raw1               = fastqc_raw.read1_seq
    Int    fastqc_raw2               = fastqc_raw.read2_seq
    String fastqc_raw_pairs          = fastqc_raw.read_pairs

    Int    fastqc_clean1             = fastqc_clean.read1_seq
    Int    fastqc_clean2             = fastqc_clean.read2_seq
    String fastqc_clean_pairs        = fastqc_clean.read_pairs

    String kraken_version            = kraken2_raw.version
    Float  kraken_human              = kraken2_raw.percent_human
    Float  kraken_sc2                = kraken2_raw.percent_sc2
    String kraken_report             = kraken2_raw.kraken_report
    Float  kraken_human_dehosted     = kraken2_dehosted.percent_human
    Float  kraken_sc2_dehosted       = kraken2_dehosted.percent_sc2
    String kraken_report_dehosted    = kraken2_dehosted.kraken_report

    String fastqc_version            = fastqc_raw.version
    String bbduk_docker              = bbduk.bbduk_docker
    String trimmomatic_version       = trimmomatic.version
  }
}
