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
  }

  call read_clean.seqyclean {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw,
  }
  call qc_utils.fastqc as fastqc_raw {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
  }
  call qc_utils.fastqc as fastqc_clean {
    input:
      read1 = seqyclean.read1_clean,
      read2 = seqyclean.read2_clean
  }
  call taxonID.kraken2 {
    input:
      samplename = samplename,
      read1 = seqyclean.read1_clean, 
      read2 = seqyclean.read2_clean
  }

  output {
  	File 	   read1_clean = seqyclean.read1_clean
  	File 	   read2_clean = seqyclean.read2_clean

    Int      fastqc_raw1      = fastqc_raw.read1_seq
    Int      fastqc_raw2      = fastqc_raw.read2_seq
    Int      fastqc_raw_pairs = fastqc_raw.read_pairs

    Int      seqy_pairs   = seqyclean.seqy_pairs
    Float    seqy_percent = seqyclean.seqy_percent

    Int      fastqc_clean1      = fastqc_clean.read1_seq
    Int      fastqc_clean2      = fastqc_clean.read2_seq
    Int      fastqc_clean_pairs = fastqc_clean.read_pairs

    Float    kraken_human  = kraken2.percent_human
    Float    kraken_sc2    = kraken2.percent_sc2
    File     kraken_report = kraken2.kraken_report

    String   fastqc_version    = fastqc_raw.version 
    String   seqyclean_version = seqyclean.version
    String   kraken_version    = kraken2.version
  }
}
