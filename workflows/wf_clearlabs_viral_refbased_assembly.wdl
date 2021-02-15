version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/task_consensus_call.wdl" as consensus_call
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_amplicon_metrics.wdl" as amplicon_metrics
import "../tasks/task_ncbi.wdl" as ncbi
import "wf_ont_sc2_pubRepo_submission.wdl" as submission

workflow viral_refbased_assembly {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }

  input {
    String  samplename
    String? artic_primer_version="V3"
    File    clear_lab_fastq
    Int?      normalise=20000

    String	   	submission_id
    String 		  collection_date
    String    	gisaid_submitter
    String    	iso_state
    String    	iso_continent
    String      iso_country
    String    	originating_lab
    String    	origLab_address
    String      BioProject
    String    	submitting_lab
    String    	subLab_address
    String    	Authors
    String    	organism = "Severe acute respiratory syndrome coronavirus 2"
    String    	iso_org = "SARS-CoV-2"
    String    	iso_host = "Human"
    String    	assembly_or_consensus = "consensus"
    String    	seq_platform = "Nanopore via Clear Labs Dx WGS SARS-CoV-2"

    String    passage_details="Original"
    String    gender="unknown"
    String    patient_age="unknown"
    String    patient_status="unknown"
    String    specimen_source=""
    String    outbreak=""
    String    last_vaccinated=""
    String    treatment=""

  }

  call medaka.consensus {
    input:
      samplename = samplename,
      filtered_reads = clear_lab_fastq,
      artic_primer_version = artic_primer_version,
      normalise = normalise


  }
  call consensus_call.variant_call {
    input:
      samplename = samplename,
      bamfile = consensus.trim_sorted_bam
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
      fasta = consensus.consensus_seq
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
  if (vadr.vadr_result) {
  call submission.SC2_submission_files as vadr_passed_submissions{
  	input:
      samplename = samplename,
		 	submission_id = submission_id,
		  collection_date = collection_date,
		  sequence = consensus.consensus_seq,
		  reads = clear_lab_fastq,

	  	organism = organism,
	  	iso_org = iso_org,
	  	iso_host = iso_host,
	  	iso_country = iso_country,
	  	assembly_or_consensus = assembly_or_consensus,

		 	gisaid_submitter = gisaid_submitter,
     	iso_state = iso_state,
     	iso_continent = iso_continent,
     	seq_platform = seq_platform,
     	artic_pipeline_version = consensus.artic_pipeline_version,
    	originating_lab = originating_lab,
    	origLab_address = origLab_address,
      BioProject = BioProject,
    	submitting_lab = submitting_lab,
    	subLab_address = subLab_address,
    	Authors = Authors,

      passage_details = passage_details,
      gender = gender,
      patient_age = patient_age,
      patient_status = patient_status,
      specimen_source = specimen_source,
      outbreak = outbreak,
      last_vaccinated = last_vaccinated,
      treatment = treatment
  }
  }
  if (! vadr.vadr_result) {
  call submission.SC2_submission_files as vadr_warning_submissions{
  	input:
      samplename = samplename,
		 	submission_id = submission_id,
		  collection_date = collection_date,
      sequence = consensus.consensus_seq,
 		  reads = clear_lab_fastq,

	  	organism = organism,
	  	iso_org = iso_org,
	  	iso_host = iso_host,
	  	iso_country = iso_country,
	  	assembly_or_consensus = assembly_or_consensus,

		 	gisaid_submitter = gisaid_submitter,
     	iso_state = iso_state,
     	iso_continent = iso_continent,
     	seq_platform = seq_platform,
     	artic_pipeline_version = consensus.artic_pipeline_version,
    	originating_lab = originating_lab,
    	origLab_address = origLab_address,
      BioProject = BioProject,
    	submitting_lab = submitting_lab,
    	subLab_address = subLab_address,
    	Authors = Authors,

      passage_details = passage_details,
      gender = gender,
      patient_age = patient_age,
      patient_status = patient_status,
      specimen_source = specimen_source,
      outbreak = outbreak,
      last_vaccinated = last_vaccinated,
      treatment = treatment
  }
  }

  output {

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
    String  artic_pipeline_version     = consensus.artic_pipeline_version

    Int     variant_num                = variant_call.variant_num

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

    File?      vadr_passed_deID_assembly      = vadr_passed_submissions.deID_assembly
    File?     vadr_passed_genbank_assembly   = vadr_passed_submissions.genbank_assembly
    File?     vadr_passed_genbank_metadata   = vadr_passed_submissions.genbank_metadata
    File?     vadr_passed_gisaid_assembly    = vadr_passed_submissions.gisaid_assembly
    File?     vadr_passed_gisaid_metadata    = vadr_passed_submissions.gisaid_metadata

    File?      vadr_warning_deID_assembly      = vadr_warning_submissions.deID_assembly
    File?     vadr_warning_genbank_assembly   = vadr_warning_submissions.genbank_assembly
    File?     vadr_warning_genbank_metadata   = vadr_warning_submissions.genbank_metadata
    File?     vadr_warning_gisaid_assembly    = vadr_warning_submissions.gisaid_assembly
    File?     vadr_warning_gisaid_metadata    = vadr_warning_submissions.gisaid_metadata

  }
}
