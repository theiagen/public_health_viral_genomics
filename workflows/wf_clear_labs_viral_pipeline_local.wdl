version 1.0

import "wf_clearlabs_viral_refbased_assembly.wdl" as assembly
import "../tasks/task_amplicon_metrics.wdl" as assembly_metrics
import "../tasks/task_clear_labs_sample_metrics.wdl" as summary
import "wf_ont_sc2_pubRepo_submission.wdl" as submission


workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], File]] inputSamples
    Array[Array[String]] inputConfig
  }

  scatter (sample in inputSamples) {
    call assembly.viral_refbased_assembly {
      input:
        samplename = sample.left[0],
        clear_lab_fastq  = sample.right
    }

    call summary.sample_metrics {
      input:
        samplename        = sample.left[0],
        submission_id     = sample.left[1],
        collection_date   = sample.left[2],
        pangolin_lineage  = viral_refbased_assembly.pangolin_lineage,
        pangolin_aLRT     = viral_refbased_assembly.pangolin_aLRT,
        nextclade_clade   = viral_refbased_assembly.nextclade_clade,
        nextclade_aa_subs = viral_refbased_assembly.nextclade_aa_subs,
        nextclade_aa_dels = viral_refbased_assembly.nextclade_aa_dels,
        kraken_human      = viral_refbased_assembly.kraken_human,
        kraken_sc2        = viral_refbased_assembly.kraken_sc2,
        variant_num       = viral_refbased_assembly.variant_num,
        number_N          = viral_refbased_assembly.number_N,
        number_ATCG       = viral_refbased_assembly.number_ATCG,
        number_Degenerate = viral_refbased_assembly.number_Degenerate,
        number_Total      = viral_refbased_assembly.number_Total,
        coverage          = viral_refbased_assembly.coverage,
        depth             = viral_refbased_assembly.depth,
        meanbaseq_trim    = viral_refbased_assembly.meanbaseq_trim,
        meanmapq_trim     = viral_refbased_assembly.meanmapq_trim,
        coverage_trim     = viral_refbased_assembly.coverage_trim,
        depth_trim        = viral_refbased_assembly.depth_trim,
        amp_fail          = viral_refbased_assembly.amp_fail
    }

  }

  call assembly_metrics.bedtools_multicov {
  	input:
  	  bamfiles          = viral_refbased_assembly.trim_sorted_bam,
  	  baifiles          = viral_refbased_assembly.trim_sorted_bai,
  	  primtrim_bamfiles = viral_refbased_assembly.trim_sorted_bam,
  	  primtrim_baifiles = viral_refbased_assembly.trim_sorted_bai
  }

  call summary.merge_metrics {
    input:
      all_metrics = sample_metrics.single_metrics
  }

  output {
    Array[File]    kraken_report        = viral_refbased_assembly.kraken_report
    Array[File]    trim_sorted_bam      = viral_refbased_assembly.trim_sorted_bam
    Array[File]    trim_sorted_bai      = viral_refbased_assembly.trim_sorted_bai
    Array[File]    consensus_seq        = viral_refbased_assembly.consensus_seq
    Array[File]    samtools_stats       = viral_refbased_assembly.consensus_stats
    Array[File]    cov_hist             = viral_refbased_assembly.cov_hist
    Array[File]    cov_stats            = viral_refbased_assembly.cov_stats
    Array[File]    samtools_flagstat    = viral_refbased_assembly.consensus_flagstat
    Array[File]    pango_lineage_report = viral_refbased_assembly.pango_lineage_report
    Array[File]    amp_coverage         = viral_refbased_assembly.amp_coverage
    File           amp_multicov         = bedtools_multicov.amp_coverage
    File           merged_metrics       = merge_metrics.run_results

    Array[File?]   reads_submission      = viral_refbased_assembly.reads_submission
    Array[File]    deID_assembly        = viral_refbased_assembly.deID_assembly
    Array[File?]   genbank_assembly     = viral_refbased_assembly.genbank_assembly
    Array[File?]   gisaid_assembly      = viral_refbased_assembly.gisaid_assembly

  }
}
