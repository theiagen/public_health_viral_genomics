version 1.0

import "wf_read_QC_trim.wdl" as read_qc
import "../tasks/task_alignment.wdl" as align
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/task_versioning.wdl" as versioning
import "../workflows/wf_WasteWaterVariantCalling_modified.wdl" as wastewater

workflow titan_illumina_wwvc {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }

  input {
    Array[String]  samplename
    Array[File]    read1_raw
    Array[File]    read2_raw
    File           primer_bed
    File           reference_genome
    File           spike_bed
    File           spike_annotations
    Int            trimmomatic_minlen = 25

  }
  scatter (r1_r2 in zip(read1_raw, read2_raw)) {
    call read_qc.read_QC_trim {
      input:
        samplename = "wastewater_sample",
        read1_raw  = r1_r2.left,
        read2_raw  = r1_r2.right,
        trimmomatic_minlen = trimmomatic_minlen
    }
    call align.bwa {
      input:
        samplename = "wastewater_sample",
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean
    }
    call consensus_call.primer_trim {
      input:
        samplename = "wastewater_sample",
        primer_bed = primer_bed,
        bamfile = bwa.sorted_bam
    }
  }
   call wastewater.WasteWaterVariantCalling{
     input:
       sorted_bam = primer_trim.trim_sorted_bam,
       covid_genome = reference_genome,
       spike_bed = spike_bed,
       spike_annotations = spike_annotations,
       sample_id = samplename
   }
   call versioning.version_capture{
     input:
   }
   output {
     String titan_wwvc_version = version_capture.phvg_version
     String titan_wwcv_date = version_capture.date
   
     Array[File] addrg_bam = WasteWaterVariantCalling.addrg_bam
     Array[File] variants = WasteWaterVariantCalling.variants
     Array[File] sorted_vcf = WasteWaterVariantCalling.sorted_vcf
     Array[File] sample_spike_vcf = WasteWaterVariantCalling.sample_spike_vcf
     Array[File] sample_spike_tsv = WasteWaterVariantCalling.sample_spike_tsv
     Array[File] sample_spike_tsv_summary = WasteWaterVariantCalling.sample_spike_tsv_summary
     Array[File] sample_spike_tsv_dash = WasteWaterVariantCalling.sample_spike_tsv_dash
     Array[File] fill_NA_tsv = WasteWaterVariantCalling.fill_NA_tsv
     Array[File] allele_freq_tsv = WasteWaterVariantCalling.allele_freq_tsv
     Array[File] reformat_tsv_tsv = WasteWaterVariantCalling.reformat_tsv_tsv
     Array[File] sample_spike_tsv_counts = WasteWaterVariantCalling.sample_spike_tsv_counts
     Array[File] alignment_files = primer_trim.trim_sorted_bam
     File spike_summary_temp = WasteWaterVariantCalling.spike_summary_temp
     File spike_summary = WasteWaterVariantCalling.spike_summary
     File spike_dashboard = WasteWaterVariantCalling.spike_dashboard
     File spike_counts = WasteWaterVariantCalling.spike_counts
       
   }
}
  
  