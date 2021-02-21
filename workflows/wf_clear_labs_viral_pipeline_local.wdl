version 1.0

import "wf_titan_ont.wdl" as assembly
import "../tasks/task_amplicon_metrics.wdl" as assembly_metrics
import "../tasks/task_clear_labs_sample_metrics.wdl" as summary
import "wf_ont_sc2_pubRepo_submission.wdl" as submission


workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], File]] inputSamples
    Array[Array[String]] inputConfig
  }

  scatter (sample in inputSamples) {
    call assembly.titan_ont {
      input:
        samplename = sample.left[0],
        demultiplexed_reads  = sample.right
    }

    call summary.sample_metrics {
      input:
        samplename        = sample.left[0],
        submission_id     = sample.left[1],
        collection_date   = sample.left[2],
        pangolin_lineage  = titan_ont.pangolin_lineage,
        pangolin_aLRT     = titan_ont.pangolin_aLRT,
        nextclade_clade   = titan_ont.nextclade_clade,
        nextclade_aa_subs = titan_ont.nextclade_aa_subs,
        nextclade_aa_dels = titan_ont.nextclade_aa_dels,
        kraken_human      = titan_ont.kraken_human,
        kraken_sc2        = titan_ont.kraken_sc2,
        variant_num       = titan_ont.variant_num,
        number_N          = titan_ont.number_N,
        number_ATCG       = titan_ont.number_ATCG,
        number_Degenerate = titan_ont.number_Degenerate,
        number_Total      = titan_ont.number_Total,
        coverage          = titan_ont.coverage,
        depth             = titan_ont.depth,
        meanbaseq_trim    = titan_ont.meanbaseq_trim,
        meanmapq_trim     = titan_ont.meanmapq_trim,
        coverage_trim     = titan_ont.coverage_trim,
        depth_trim        = titan_ont.depth_trim,
        amp_fail          = titan_ont.amp_fail
    }

  }

  call assembly_metrics.bedtools_multicov {
  	input:
  	  bamfiles          = titan_ont.trim_sorted_bam,
  	  baifiles          = titan_ont.trim_sorted_bai,
  	  primtrim_bamfiles = titan_ont.trim_sorted_bam,
  	  primtrim_baifiles = titan_ont.trim_sorted_bai
  }

  call summary.merge_metrics {
    input:
      all_metrics = sample_metrics.single_metrics
  }

  output {
    Array[File]    kraken_report        = titan_ont.kraken_report
    Array[File]    trim_sorted_bam      = titan_ont.trim_sorted_bam
    Array[File]    trim_sorted_bai      = titan_ont.trim_sorted_bai
    Array[File]    consensus_seq        = titan_ont.consensus_seq
    Array[File]    samtools_stats       = titan_ont.consensus_stats
    Array[File]    cov_hist             = titan_ont.cov_hist
    Array[File]    cov_stats            = titan_ont.cov_stats
    Array[File]    samtools_flagstat    = titan_ont.consensus_flagstat
    Array[File]    pango_lineage_report = titan_ont.pango_lineage_report
    Array[File]    amp_coverage         = titan_ont.amp_coverage
    File           amp_multicov         = bedtools_multicov.amp_coverage
    File           merged_metrics       = merge_metrics.run_results


  }
}
