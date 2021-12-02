version 1.0

import "../tasks/task_taxonID.wdl" as taxon_id
import "wf_read_QC_trim.wdl" as read_qc
import "../tasks/task_alignment.wdl" as align
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/task_versioning.wdl" as versioning

workflow freyja_fastq {
  input {
    File read1_raw
    File read2_raw
    File primer_bed
    Int trimmomatic_minlen = 25
    String samplename
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
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean
  }
  call consensus_call.primer_trim {
    input:
      samplename = samplename,
      primer_bed = primer_bed,
      bamfile = bwa.sorted_bam
  }
  call taxon_id.freyja_one_sample as freyja {
    input:
      primer_trimmed_bam = primer_trim.trim_sorted_bam,
      samplename = samplename
  }
  call versioning.version_capture{
    input:
  }
  output {
    String titan_freyja_run_wf_version = version_capture.phvg_version
    String titan_freyja_run_wf_analysis_date = version_capture.date
    
    File read1_dehosted = read_QC_trim.read1_dehosted
    File read2_dehosted = read_QC_trim.read2_dehosted
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    Int fastqc_raw1 = read_QC_trim.fastqc_raw1
    Int fastqc_raw2 = read_QC_trim.fastqc_raw2
    String fastqc_raw_pairs = read_QC_trim.fastqc_raw_pairs
    String fastqc_version = read_QC_trim.fastqc_version

    Int fastqc_clean1 = read_QC_trim.fastqc_clean1
    Int fastqc_clean2 = read_QC_trim.fastqc_clean2
    String fastqc_clean_pairs = read_QC_trim.fastqc_clean_pairs
    String trimmomatic_version = read_QC_trim.trimmomatic_version
    String bbduk_docker = read_QC_trim.bbduk_docker

    String kraken_version = read_QC_trim.kraken_version
    Float kraken_human = read_QC_trim.kraken_human
    Float kraken_sc2 = read_QC_trim.kraken_sc2
    String kraken_report = read_QC_trim.kraken_report
    Float kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    
    String bwa_version = bwa.bwa_version
    String samtools_version = bwa.sam_version
    String alignment_method = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"

    File aligned_bam = primer_trim.trim_sorted_bam
    File aligned_bai = primer_trim.trim_sorted_bai
    Float primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String ivar_version_primtrim = primer_trim.ivar_version
    String samtools_version_primtrim = primer_trim.samtools_version
    String primer_bed_name = primer_trim.primer_bed_name

    File freyja_variants = freyja.freyja_variants
    File freyja_depths = freyja.freyja_depths
    File freyja_demixed = freyja.freyja_demixed
    }
}