version 1.0

import "wf_titan_clearlabs.wdl" as assembly
import "../tasks/task_amplicon_metrics.wdl" as assembly_metrics
import "../tasks/task_clear_labs_sample_metrics.wdl" as summary

workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], File]] inputSamples
    Array[Array[String]] inputConfig
  }

  scatter (sample in inputSamples) {
    call assembly.titan_clearlabs {
      input:
        samplename = sample.left[0],
        clear_lab_fastq  = sample.right
    }

    call summary.sample_metrics {
      input:
        samplename        = sample.left[0],
        submission_id     = sample.left[1],
        collection_date   = sample.left[2],
        pangolin_lineage  = titan_clearlabs.pangolin_lineage,
        pangolin_aLRT     = titan_clearlabs.pangolin_aLRT,
        nextclade_clade   = titan_clearlabs.nextclade_clade,
        nextclade_aa_subs = titan_clearlabs.nextclade_aa_subs,
        nextclade_aa_dels = titan_clearlabs.nextclade_aa_dels,
        kraken_human      = titan_clearlabs.kraken_human,
        kraken_sc2        = titan_clearlabs.kraken_sc2,
        number_N          = titan_clearlabs.number_N,
        number_ATCG       = titan_clearlabs.number_ATCG,
        number_Degenerate = titan_clearlabs.number_Degenerate,
        number_Total      = titan_clearlabs.number_Total,
        coverage          = titan_clearlabs.coverage,
        depth             = titan_clearlabs.depth,
        meanbaseq_trim    = titan_clearlabs.meanbaseq_trim,
        meanmapq_trim     = titan_clearlabs.meanmapq_trim,
        coverage_trim     = titan_clearlabs.coverage_trim,
        depth_trim        = titan_clearlabs.depth_trim,
        amp_fail          = titan_clearlabs.amp_fail
    }

  }

  call assembly_metrics.bedtools_multicov {
  	input:
  	  bamfiles          = titan_clearlabs.trim_sorted_bam,
  	  baifiles          = titan_clearlabs.trim_sorted_bai,
  	  primtrim_bamfiles = titan_clearlabs.trim_sorted_bam,
  	  primtrim_baifiles = titan_clearlabs.trim_sorted_bai
  }

  call summary.merge_metrics {
    input:
      all_metrics = sample_metrics.single_metrics
  }

  output {
    Array[File]    kraken_report        = titan_clearlabs.kraken_report
    Array[File]    trim_sorted_bam      = titan_clearlabs.trim_sorted_bam
    Array[File]    trim_sorted_bai      = titan_clearlabs.trim_sorted_bai
    Array[File]    consensus_seq        = titan_clearlabs.consensus_seq
    Array[File]    samtools_stats       = titan_clearlabs.consensus_stats
    Array[File]    cov_hist             = titan_clearlabs.cov_hist
    Array[File]    cov_stats            = titan_clearlabs.cov_stats
    Array[File]    samtools_flagstat    = titan_clearlabs.consensus_flagstat
    Array[File]    pango_lineage_report = titan_clearlabs.pango_lineage_report
    Array[File]    amp_coverage         = titan_clearlabs.amp_coverage
    File           amp_multicov         = bedtools_multicov.amp_coverage
    File           merged_metrics       = merge_metrics.run_results


  }
}
