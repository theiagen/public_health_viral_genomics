version 1.0

import "wf_titan_illumina_pe.wdl" as assembly
import "../tasks/task_amplicon_metrics.wdl" as assembly_metrics
import "../tasks/task_sample_metrics.wdl" as summary


workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], Pair[File,File]]] inputSamples
    Array[Array[String]] inputConfig
  }

  scatter (sample in inputSamples) {
    call assembly.titan_illumina_pe {
      input:
        samplename = sample.left[0],
        read1_raw  = sample.right.left,
        read2_raw  = sample.right.right
    }

    call summary.sample_metrics {
      input:
        samplename        = sample.left[0],
        submission_id     = sample.left[1],
        collection_date   = sample.left[2],
        pangolin_lineage  = titan_illumina_pe.pangolin_lineage,
        pangolin_aLRT     = titan_illumina_pe.pangolin_aLRT,
        nextclade_clade   = titan_illumina_pe.nextclade_clade,
        nextclade_aa_subs = titan_illumina_pe.nextclade_aa_subs,
        nextclade_aa_dels = titan_illumina_pe.nextclade_aa_dels,
        fastqc_raw_pairs  = titan_illumina_pe.fastqc_raw_pairs,
        seqy_pairs        = titan_illumina_pe.seqy_pairs,
        seqy_percent      = titan_illumina_pe.seqy_percent,
        kraken_human      = titan_illumina_pe.kraken_human,
        kraken_sc2        = titan_illumina_pe.kraken_sc2,
        variant_num       = titan_illumina_pe.variant_num,
        number_N          = titan_illumina_pe.number_N,
        number_ATCG       = titan_illumina_pe.number_ATCG,
        number_Degenerate = titan_illumina_pe.number_Degenerate,
        number_Total      = titan_illumina_pe.number_Total,
        coverage          = titan_illumina_pe.coverage,
        depth             = titan_illumina_pe.depth,
        meanbaseq_trim    = titan_illumina_pe.meanbaseq_trim,
        meanmapq_trim     = titan_illumina_pe.meanmapq_trim,
        coverage_trim     = titan_illumina_pe.coverage_trim,
        depth_trim        = titan_illumina_pe.depth_trim,
        amp_fail          = titan_illumina_pe.amp_fail
    }

  }

  call assembly_metrics.bedtools_multicov {
  	input:
  	  bamfiles          = titan_illumina_pe.sorted_bam,
  	  baifiles          = titan_illumina_pe.sorted_bai,
  	  primtrim_bamfiles = titan_illumina_pe.trim_sorted_bam,
  	  primtrim_baifiles = titan_illumina_pe.trim_sorted_bai
  }

  call summary.merge_metrics {
    input:
      all_metrics = sample_metrics.single_metrics
  }

  output {
    Array[File]    read1_clean          = titan_illumina_pe.read1_clean
    Array[File]    read2_clean          = titan_illumina_pe.read2_clean
    Array[File]    kraken_report        = titan_illumina_pe.kraken_report
    Array[File]    sorted_bam           = titan_illumina_pe.sorted_bam
    Array[File]    sorted_bai           = titan_illumina_pe.sorted_bai
    Array[File]    trim_sorted_bam      = titan_illumina_pe.trim_sorted_bam
    Array[File]    trim_sorted_bai      = titan_illumina_pe.trim_sorted_bai
    Array[File]    consensus_seq        = titan_illumina_pe.consensus_seq
    Array[File]    samtools_stats       = titan_illumina_pe.consensus_stats
    Array[File]    cov_hist             = titan_illumina_pe.cov_hist
    Array[File]    cov_stats            = titan_illumina_pe.cov_stats
    Array[File]    samtools_flagstat    = titan_illumina_pe.consensus_flagstat
    Array[File]    pango_lineage_report = titan_illumina_pe.pango_lineage_report
    Array[File]    amp_coverage         = titan_illumina_pe.amp_coverage
    File           amp_multicov         = bedtools_multicov.amp_coverage
    File           merged_metrics       = merge_metrics.run_results

  }
}
