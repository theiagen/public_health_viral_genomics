version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/quality_control/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../tasks/task_versioning.wdl" as versioning
import "../tasks/task_sc2_gene_coverage.wdl" as sc2_calculation

workflow theiacov_clearlabs {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }
  input {
    String samplename
    File clear_lab_fastq
    String seq_method = "OXFORD_NANOPORE"
    File primer_bed
    Int normalise = 20000
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2022-09-27T12:00:00Z"
    String medaka_docker = "quay.io/staphb/artic-ncov2019:1.3.0-medaka-1.4.3"
    String? nextclade_dataset_name
    File? reference_genome
    String organism = "sars-cov-2"
    String? target_org
  }
  call fastq_scan.fastq_scan_se as fastq_scan_raw_reads {
    input:
      read1 = clear_lab_fastq
  }
  call read_clean.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = clear_lab_fastq
  }
  call fastq_scan.fastq_scan_se as fastq_scan_clean_reads {
    input:
      read1 = ncbi_scrub_se.read1_dehosted
  }
  call taxon_ID.kraken2 as kraken2_raw {
    input:
      samplename = samplename,
      read1 = clear_lab_fastq,
      target_org = target_org
  }  
  call taxon_ID.kraken2 as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted,
      target_org = target_org
  }
  call medaka.consensus {
    input:
      samplename = samplename,
      filtered_reads = ncbi_scrub_se.read1_dehosted,
      primer_bed = primer_bed,
      normalise = normalise,
      docker = medaka_docker,
      reference_genome = reference_genome,
      organism = organism
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = consensus.sorted_bam
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = consensus.consensus_seq,
      reference_genome = reference_genome
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = consensus.trim_sorted_bam
  }
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call taxon_ID.pangolin4 {
      input:
        samplename = samplename,
        fasta = consensus.consensus_seq
    }
    call sc2_calculation.sc2_gene_coverage {
      input: 
        samplename = samplename,
        bamfile = consensus.trim_sorted_bam,
        min_depth = 20
    }
  }
  if (organism == "MPXV") {
    # MPXV specific tasks
  }
  if (organism == "MPXV" || organism == "sars-cov-2"){
    # tasks specific to either MPXV or sars-cov-2
    call taxon_ID.nextclade_one_sample {
      input:
      genome_fasta = consensus.consensus_seq,
      dataset_name = select_first([nextclade_dataset_name, organism]),
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
    }
    call taxon_ID.nextclade_output_parser_one_sample {
      input:
      nextclade_tsv = nextclade_one_sample.nextclade_tsv
    }
    call ncbi.vadr {
      input:
        genome_fasta = consensus.consensus_seq,
        assembly_length_unambiguous = consensus_qc.number_ATCG
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_clearlabs_version = version_capture.phvg_version
    String theiacov_clearlabs_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC
    File reads_dehosted = ncbi_scrub_se.read1_dehosted
    Int num_reads_raw = fastq_scan_raw_reads.read1_seq
    Int num_reads_clean = fastq_scan_clean_reads.read1_seq
    String fastq_scan_version = fastq_scan_raw_reads.version
    String kraken_version = kraken2_raw.version
    Float kraken_human = kraken2_raw.percent_human
    Float kraken_sc2 = kraken2_raw.percent_sc2
    String? kraken_target_org = kraken2_raw.percent_target_org
    String? kraken_target_org_name = kraken2_raw.kraken_target_org
    File kraken_report = kraken2_raw.kraken_report
    Float kraken_human_dehosted = kraken2_dehosted.percent_human
    Float kraken_sc2_dehosted = kraken2_dehosted.percent_sc2
    String? kraken_target_org_dehosted = kraken2_dehosted.percent_target_org
    File kraken_report_dehosted = kraken2_dehosted.kraken_report
    # Read Alignment
    File aligned_bam = consensus.trim_sorted_bam
    File aligned_bai = consensus.trim_sorted_bai
    File variants_from_ref_vcf = consensus.medaka_pass_vcf
    String artic_version = consensus.artic_pipeline_version
    String artic_docker = consensus.artic_pipeline_docker
    String medaka_reference = consensus.medaka_reference
    String primer_bed_name = consensus.primer_bed_name
    File assembly_fasta = consensus.consensus_seq
    String assembly_method = "TheiaCoV (~{version_capture.phvg_version}): ~{consensus.artic_pipeline_version}"
    File reads_aligned = consensus.reads_aligned
    # Assembly QC
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Lineage Assignment
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report= pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Alignment QC
    File consensus_stats = stats_n_coverage.stats
    File consensus_flagstat = stats_n_coverage.flagstat
    Float meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    Float assembly_mean_coverage = stats_n_coverage_primtrim.depth
    String samtools_version_stats = stats_n_coverage.samtools_version
    # SC2 specific
    Float? sc2_s_gene_mean_coverage = sc2_gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = sc2_gene_coverage.sc2_s_gene_percent_coverage
    File? sc2_all_genes_percent_coverage = sc2_gene_coverage.sc2_all_genes_percent_coverage
    # Clade Assigment
    File? nextclade_json = nextclade_one_sample.nextclade_json
    File? auspice_json = nextclade_one_sample.auspice_json
    File? nextclade_tsv = nextclade_one_sample.nextclade_tsv
    String? nextclade_version = nextclade_one_sample.nextclade_version
    String? nextclade_docker = nextclade_one_sample.nextclade_docker
    String nextclade_ds_tag = nextclade_dataset_tag
    String? nextclade_aa_subs = nextclade_output_parser_one_sample.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser_one_sample.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser_one_sample.nextclade_clade
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
  }
}
