version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_amplicon_metrics.wdl" as amplicon_metrics
import "../tasks/task_ncbi.wdl" as ncbi

workflow titan_clearlabs {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }

  input {
    String  samplename
    String  seq_method="ONT via Clear Labs WGS"
    String? artic_primer_version="V3"
    String  pangolin_docker_image = "staphb/pangolin:2.3.2-pangolearn-2021-02-21"
    File  clear_lab_fastq
    Int?  normalise=20000

  }

  call medaka.consensus {
    input:
      samplename = samplename,
      filtered_reads = clear_lab_fastq,
      artic_primer_version = artic_primer_version,
      normalise = normalise
  }

  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = consensus.sorted_bam
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = consensus.trim_sorted_bam
  }
  call taxon_ID.pangolin2 {
    input:
      samplename = samplename,
      fasta = consensus.consensus_seq,
      docker = pangolin_docker_image
  }
  call taxon_ID.kraken2 {
    input:
      samplename = samplename,
      read1 = clear_lab_fastq
  }
  call taxon_ID.nextclade_one_sample {
    input:
      genome_fasta = consensus.consensus_seq
  }
  call amplicon_metrics.bedtools_cov {
    input:
      bamfile = consensus.trim_sorted_bam,
      baifile = consensus.trim_sorted_bai
  }
  call ncbi.vadr {
    input:
      genome_fasta = consensus.consensus_seq,
      samplename = samplename
  }

  output {

    String  seq_platform = seq_method

    Float   kraken_human       = kraken2.percent_human
    Float   kraken_sc2         = kraken2.percent_sc2
    String  kraken_version     = kraken2.version
    String  kraken_report      = kraken2.kraken_report

    File    trim_sorted_bam         = consensus.trim_sorted_bam
    File    trim_sorted_bai         = consensus.trim_sorted_bai
    String  artic_version           = consensus.artic_pipeline_version
    File    consensus_seq              = consensus.consensus_seq
    Int     number_N                   = consensus.number_N
    Int     number_ATCG                = consensus.number_ATCG
    Int     number_Degenerate          = consensus.number_Degenerate
    Int     number_Total               = consensus.number_Total
    Float   pool1_percent              = consensus.pool1_percent
    Float   pool2_percent              = consensus.pool2_percent
    String  assembly_method     = consensus.artic_pipeline_version

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
