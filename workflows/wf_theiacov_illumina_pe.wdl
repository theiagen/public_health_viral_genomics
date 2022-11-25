version 1.0

import "wf_read_QC_trim.wdl" as read_qc
import "../tasks/task_alignment.wdl" as align
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/quality_control/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_versioning.wdl" as versioning
import "../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../tasks/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../tasks/task_irma.wdl" as irma
import "../tasks/task_abricate_flu.wdl" as abricate_flu

workflow theiacov_illumina_pe {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File read2_raw
    File? primer_bed
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2022-09-27T12:00:00Z"
    String? nextclade_dataset_name
    File? reference_genome
    Int min_depth = 100
    String organism = "sars-cov-2"
    Boolean trim_primers = true
    File? adapters
    File? phix
    String nextclade_flu_h1n1_tag = "2022-06-08T12:00:00Z"
    String nextclade_flu_h3n2_tag = "2022-06-08T12:00:00Z"
    String nextclade_flu_vic_tag = "2022-06-08T12:00:00Z"
    String nextclade_flu_yam_tag = "2022-07-27T12:00:00Z"
    Int? genome_length
  }
  call read_qc.read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = read1_raw,
      read2_raw = read2_raw,
      adapters = adapters,
      phix = phix
  }
  if (organism != "flu"){
    call align.bwa {
      input:
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        reference_genome = reference_genome
    }
    if (trim_primers){
      call consensus_call.primer_trim {
        input:
          samplename = samplename,
          primer_bed = select_first([primer_bed]),
          bamfile = bwa.sorted_bam
      }
      call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
      input:
        samplename = samplename,
        bamfile = primer_trim.trim_sorted_bam,
      }
    }
    call consensus_call.variant_call {
      input:
        samplename = samplename,
        bamfile = select_first([primer_trim.trim_sorted_bam,bwa.sorted_bam]),
        reference_genome = reference_genome,
        variant_min_depth = min_depth
    }
    call consensus_call.consensus {
      input:
        samplename = samplename,
        bamfile = select_first([primer_trim.trim_sorted_bam,bwa.sorted_bam]),
        reference_genome = reference_genome,
        consensus_min_depth = min_depth
    }
    call consensus_qc_task.consensus_qc as consensus_qc {
        input:
          assembly_fasta =  consensus.consensus_seq,
          reference_genome = reference_genome,
          genome_length = genome_length
      }
    call assembly_metrics.stats_n_coverage {
      input:
        samplename = samplename,
        bamfile = bwa.sorted_bam
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
          bamfile = bwa.sorted_bam,
          min_depth = min_depth
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
    String assembly_method_nonflu = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
  }
  if (organism == "flu"){
    # flu-specific tasks
    call irma.irma {
      input:
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        samplename = samplename,
    }
    call abricate_flu.abricate_flu {
      input:
        assembly = irma.irma_assembly_fasta,
        samplename = samplename,
        nextclade_flu_h1n1_tag = nextclade_flu_h1n1_tag,
        nextclade_flu_h3n2_tag = nextclade_flu_h3n2_tag,
        nextclade_flu_vic_tag = nextclade_flu_vic_tag,
        nextclade_flu_yam_tag = nextclade_flu_yam_tag,
    }
    if (abricate_flu.run_nextclade && irma.ha_seg_exists) {
      call taxon_ID.nextclade_one_sample as nextclade_one_sample_flu {
        input:
          genome_fasta = select_first([irma.seg4_ha_assembly]),
          dataset_name = abricate_flu.nextclade_name,
          dataset_reference = abricate_flu.nextclade_ref,
          dataset_tag = abricate_flu.nextclade_ds_tag
      }
      call taxon_ID.nextclade_output_parser_one_sample as  nextclade_output_parser_one_sample_flu {
        input:
        nextclade_tsv = nextclade_one_sample_flu.nextclade_tsv
      }
      call consensus_qc_task.consensus_qc as consensus_qc_flu {
        input:
          assembly_fasta =  irma.irma_assembly_fasta,
          reference_genome = reference_genome,
          genome_length = genome_length
      }
    }
  }

  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_illumina_pe_version = version_capture.phvg_version
    String theiacov_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String  seq_platform = seq_method
    # Read QC
    File read1_dehosted = read_QC_trim.read1_dehosted
    File read2_dehosted = read_QC_trim.read2_dehosted
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    Int num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    Int num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String trimmomatic_version = read_QC_trim.trimmomatic_version
    String bbduk_docker = read_QC_trim.bbduk_docker
    String kraken_version = read_QC_trim.kraken_version
    Float kraken_human = read_QC_trim.kraken_human
    Float kraken_sc2 = read_QC_trim.kraken_sc2
    String? kraken_target_org = read_QC_trim.kraken_target_org
    String? kraken_target_org_name = read_QC_trim.kraken_target_org_name
    File kraken_report = read_QC_trim.kraken_report
    Float kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String? kraken_target_org_dehosted =read_QC_trim.kraken_target_org_dehosted
    File kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment
    String? bwa_version = bwa.bwa_version
    String? samtools_version = bwa.sam_version
    File? read1_aligned = bwa.read1_aligned
    File? read2_aligned = bwa.read2_aligned
    String assembly_method = "TheiaCoV (~{version_capture.phvg_version}): " + select_first([assembly_method_nonflu,irma.irma_version])
    String aligned_bam =  select_first([primer_trim.trim_sorted_bam,bwa.sorted_bam,""])
    String aligned_bai = select_first([primer_trim.trim_sorted_bai,bwa.sorted_bai,""])
    Float? primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String? ivar_version_primtrim = primer_trim.ivar_version
    String? samtools_version_primtrim = primer_trim.samtools_version
    String? primer_bed_name = primer_trim.primer_bed_name
    File? ivar_tsv = variant_call.sample_variants_tsv
    File? ivar_vcf = variant_call.sample_variants_vcf
    String? ivar_variant_version = variant_call.ivar_version
    # Assembly QC
    File assembly_fasta = select_first([consensus.consensus_seq,irma.irma_assembly_fasta])
    String? ivar_version_consensus = consensus.ivar_version
    String? samtools_version_consensus = consensus.samtools_version
    Int number_N = select_first([consensus_qc.number_N,consensus_qc_flu.number_N])
    Int assembly_length_unambiguous =  select_first([consensus_qc.number_ATCG,consensus_qc_flu.number_ATCG])
    Int number_Degenerate =  select_first([consensus_qc.number_Degenerate,consensus_qc_flu.number_Degenerate])
    Int number_Total =  select_first([consensus_qc.number_Total,consensus_qc_flu.number_Total])
    Float percent_reference_coverage =  select_first([consensus_qc.percent_reference_coverage,consensus_qc_flu.percent_reference_coverage])
    Int consensus_n_variant_min_depth = min_depth
    # Alignment QC
    File? consensus_stats = stats_n_coverage.stats
    File? consensus_flagstat = stats_n_coverage.flagstat
    String meanbaseq_trim = select_first([stats_n_coverage_primtrim.meanbaseq, stats_n_coverage.meanbaseq,""])
    String meanmapq_trim = select_first([stats_n_coverage_primtrim.meanmapq, stats_n_coverage.meanmapq,""])
    String assembly_mean_coverage = select_first([stats_n_coverage_primtrim.depth, stats_n_coverage.depth,""])
    String? samtools_version_stats = stats_n_coverage.samtools_version
    # SC2 specific
    Float? sc2_s_gene_mean_coverage = sc2_gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = sc2_gene_coverage.sc2_s_gene_percent_coverage
    File? sc2_all_genes_percent_coverage = sc2_gene_coverage.sc2_all_genes_percent_coverage
    # SC2 Lineage Assignment
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Clade Assigment
    File nextclade_json = select_first([nextclade_one_sample.nextclade_json, nextclade_one_sample_flu.nextclade_json])
    File auspice_json = select_first([ nextclade_one_sample.auspice_json, nextclade_one_sample_flu.auspice_json])
    File nextclade_tsv = select_first([nextclade_one_sample.nextclade_tsv, nextclade_one_sample_flu.nextclade_tsv])
    String nextclade_version = select_first([nextclade_one_sample.nextclade_version, nextclade_one_sample_flu.nextclade_version])
    String nextclade_docker = select_first([nextclade_one_sample.nextclade_docker, nextclade_one_sample_flu.nextclade_docker])
    String nextclade_ds_tag = select_first([abricate_flu.nextclade_ds_tag, nextclade_dataset_tag])
    String nextclade_aa_subs = select_first([nextclade_output_parser_one_sample.nextclade_aa_subs, nextclade_output_parser_one_sample_flu.nextclade_aa_subs])
    String nextclade_aa_dels = select_first([nextclade_output_parser_one_sample.nextclade_aa_dels, nextclade_output_parser_one_sample_flu.nextclade_aa_dels])
    String nextclade_clade = select_first([nextclade_output_parser_one_sample.nextclade_clade, nextclade_output_parser_one_sample_flu.nextclade_clade])
    String? nextclade_lineage = nextclade_output_parser_one_sample.nextclade_lineage
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # Flu Outputs
    String? irma_version = irma.irma_version
    String? irma_subtype = irma.irma_subtype
    String? abricate_flu_type = abricate_flu.abricate_flu_type
    String? abricate_flu_subtype =  abricate_flu.abricate_flu_subtype
    File? abricate_flu_results = abricate_flu.abricate_flu_results
    String? abricate_flu_database =  abricate_flu.abricate_flu_database
    String? abricate_flu_version = abricate_flu.abricate_flu_version
  }
}
