version 1.0

import "wf_read_QC_trim.wdl" as read_qc
import "../tasks/task_alignment.wdl" as align
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_amplicon_metrics.wdl" as amplicon_metrics
import "../tasks/task_ncbi.wdl" as ncbi


workflow titan_illumina_pe {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }

  input {
    String  samplename
    String  seq_method="Illumina paired-end"
    File    read1_raw
    File    read2_raw
    File    primer_bed
    String  pangolin_docker_image = "staphb/pangolin:2.3.2-pangolearn-2021-02-21"

  }

  call read_qc.read_QC_trim {
    input:
      samplename = samplename,
      read1_raw  = read1_raw,
      read2_raw  = read2_raw
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
  call consensus_call.variant_call {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam
  }
  call consensus_call.consensus {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = bwa.sorted_bam
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = primer_trim.trim_sorted_bam
  }
  call taxon_ID.pangolin2 {
    input:
      samplename = samplename,
      fasta = consensus.consensus_seq,
      docker = pangolin_docker_image
  }
  call taxon_ID.nextclade_one_sample {
    input:
      genome_fasta = consensus.consensus_seq
  }
  call amplicon_metrics.bedtools_cov {
    input:
      bamfile = bwa.sorted_bam,
      baifile = bwa.sorted_bai
  }
  call ncbi.vadr {
    input:
      genome_fasta = consensus.consensus_seq,
      samplename = samplename
  }
  output {

    String  seq_platform = seq_method

    File    read1_clean        = read_QC_trim.read1_clean
    File    read2_clean        = read_QC_trim.read2_clean
    Int     fastqc_raw1        = read_QC_trim.fastqc_raw1
    Int     fastqc_raw2        = read_QC_trim.fastqc_raw2
    Int     fastqc_raw_pairs   = read_QC_trim.fastqc_raw_pairs
    String  fastqc_version     = read_QC_trim.fastqc_version

    Int     seqy_pairs         = read_QC_trim.seqy_pairs
    Float   seqy_percent       = read_QC_trim.seqy_percent
    Int     fastqc_clean1      = read_QC_trim.fastqc_clean1
    Int     fastqc_clean2      = read_QC_trim.fastqc_clean2
    Int     fastqc_clean_pairs = read_QC_trim.fastqc_clean_pairs
    String  seqyclean_version  = read_QC_trim.seqyclean_version

    Float   kraken_human       = read_QC_trim.kraken_human
    Float   kraken_sc2         = read_QC_trim.kraken_sc2
    String  kraken_version     = read_QC_trim.kraken_version
    String  kraken_report      = read_QC_trim.kraken_report

    File    sorted_bam         = bwa.sorted_bam
    File    sorted_bai         = bwa.sorted_bai
    String  bwa_version        = bwa.bwa_version
    String  sam_version        = bwa.sam_version
    String assembly_method     = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"

    File    trim_sorted_bam            = primer_trim.trim_sorted_bam
    File    trim_sorted_bai            = primer_trim.trim_sorted_bai
    String  ivar_version_primtrim      = primer_trim.ivar_version
    String  samtools_version_primtrim  = primer_trim.samtools_version

    File    consensus_seq              = consensus.consensus_seq
    Int     number_N                   = consensus.number_N
    Int     number_ATCG                = consensus.number_ATCG
    Int     number_Degenerate          = consensus.number_Degenerate
    Int     number_Total               = consensus.number_Total
    String  ivar_version_consensus     = consensus.ivar_version
    String  samtools_version_consensus = consensus.samtools_version

    File    consensus_stats        = stats_n_coverage.stats
    File    cov_hist               = stats_n_coverage.cov_hist
    File    cov_stats              = stats_n_coverage.cov_stats
    File    consensus_flagstat     = stats_n_coverage.flagstat
    Float   coverage               = stats_n_coverage.coverage
    Float   depth                  = stats_n_coverage.depth
    Float   meanbaseq_trim         = stats_n_coverage_primtrim.meanbaseq
    Float   meanmapq_trim          = stats_n_coverage_primtrim.meanmapq
    Float   coverage_trim          = stats_n_coverage_primtrim.coverage
    Float   depth_trim             = stats_n_coverage_primtrim.depth
    String  samtools_version_stats = stats_n_coverage.samtools_version

    String  pangolin_lineage       = pangolin2.pangolin_lineage
    Float   pangolin_aLRT          = pangolin2.pangolin_aLRT
    File    pango_lineage_report   = pangolin2.pango_lineage_report
    String  pangolin_version       = pangolin2.version
    String  pangolin_docker       = pangolin2.pangolin_docker


    File    nextclade_json         = nextclade_one_sample.nextclade_json
    File    auspice_json           = nextclade_one_sample.auspice_json
    File    nextclade_tsv          = nextclade_one_sample.nextclade_tsv
    String  nextclade_clade        = nextclade_one_sample.nextclade_clade
    String  nextclade_aa_subs      = nextclade_one_sample.nextclade_aa_subs
    String  nextclade_aa_dels      = nextclade_one_sample.nextclade_aa_dels
    String  nextclade_version      = nextclade_one_sample.nextclade_version

    Int     amp_fail               = bedtools_cov.amp_fail
    File    amp_coverage           = bedtools_cov.amp_coverage
    String  bedtools_version       = bedtools_cov.version

    File vadr_alterts_list = vadr.alerts_list
    Int vadr_num_alerts = vadr.num_alerts
  }
}
