version 1.0

import "wf_titan_illumina_se.wdl" as assembly
import "../tasks/task_amplicon_metrics.wdl" as assembly_metrics
import "../tasks/task_sample_metrics.wdl" as summary


workflow nCoV19_pipeline {
input {
  Array[Pair[Array[String], File]] inputSamples
  File primer_bed
}

  scatter (sample in inputSamples) {
    call assembly.titan_illumina_se {
      input:
        samplename = sample.left[0],
        primer_bed = primer_bed,
        read1_raw  = sample.right
    }

    call summary.sample_metrics {
      input:
        samplename        = sample.left[0],
        submission_id     = sample.left[1],
        collection_date   = sample.left[2],
        pangolin_lineage  = titan_illumina_se.pango_lineage,
        pangolin_aLRT     = titan_illumina_se.pangolin_aLRT,
        nextclade_clade   = titan_illumina_se.nextclade_clade,
        nextclade_aa_subs = titan_illumina_se.nextclade_aa_subs,
        nextclade_aa_dels = titan_illumina_se.nextclade_aa_dels,
        fastqc_raw_pairs  = titan_illumina_se.fastqc_number_reads,
        kraken_human      = titan_illumina_se.kraken_human,
        kraken_sc2        = titan_illumina_se.kraken_sc2,
        number_N          = titan_illumina_se.number_N,
        number_ATCG       = titan_illumina_se.assembly_length_unambiguous,
        number_Degenerate = titan_illumina_se.number_Degenerate,
        number_Total      = titan_illumina_se.number_Total,
        meanbaseq_trim    = titan_illumina_se.meanbaseq_trim,
        meanmapq_trim     = titan_illumina_se.meanmapq_trim,
        coverage_trim     = titan_illumina_se.percent_reference_coverage,
        depth_trim        = titan_illumina_se.assembly_mean_coverage,
    }

  }

  call summary.merge_metrics {
    input:
      all_metrics = sample_metrics.single_metrics
  }

  output {
    Array[File]    read1_clean          = titan_illumina_se.read1_clean
    Array[File]    kraken_report        = titan_illumina_se.kraken_report
    Array[File]    trim_sorted_bam      = titan_illumina_se.aligned_bam
    Array[File]    trim_sorted_bai      = titan_illumina_se.aligned_bai
    Array[File]    consensus_seq        = titan_illumina_se.assembly_fasta
    Array[File]    samtools_stats       = titan_illumina_se.consensus_stats
    Array[File]    samtools_flagstat    = titan_illumina_se.consensus_flagstat
    File           merged_metrics       = merge_metrics.run_results

  }
}
