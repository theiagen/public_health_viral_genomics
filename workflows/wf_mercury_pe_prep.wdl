version 1.0

import "../tasks/task_pe_pub_repo_submission.wdl" as submission

workflow mercury_pe_prep {
	input {
		String 		samplename
		String		submission_id
		String 		collection_date
		File 		sequence
		File 		read1
		File 		read2

		String    	organism = "Severe acute respiratory syndrome coronavirus 2"
	  String    	iso_org = "SARS-CoV-2"
	  String    	iso_host = "Human"
	  String    	assembly_or_consensus = "consensus"

		String    	gisaid_submitter
    String    	iso_state
	  String    	iso_country
    String    	iso_continent
    String    	seq_platform
    String    	assembly_method
    String    	originating_lab
    String    	origLab_address
    String      BioProject
    String    	submitting_lab
    String    	subLab_address
    String    	Authors

		String    passage_details="Original"
		String    gender="unknown"
		String    patient_age="unknown"
		String    patient_status="unknown"
		String    specimen_source=""
		String    outbreak=""
		String    last_vaccinated=""
		String    treatment=""
		String 		iso_county=""

	  # Optional inputs/user-defined thresholds for generating submission files
		Float		coverage = 100.00
		Float		coverage_threshold = 85.00
		Int 		number_N_threshold = 15000
		Int 		number_Total_threshold = 25000
		Int 		number_ATCG_gisaid = 25000
		Int 		number_ATCG_genbank = 25000
	}

	call submission.sra {
		input:
			submission_id = submission_id,
			read1 = read1,
			read2 = read2
	}

	call submission.deidentify {
		input:
			samplename    = samplename,
			submission_id = submission_id,
			sequence      = sequence
	}

	if (coverage >= coverage_threshold) {
		if (deidentify.number_N <= number_N_threshold) {
			if (deidentify.number_Total >= number_Total_threshold) {
				if (deidentify.number_ATCG >= number_ATCG_gisaid) {
					call submission.gisaid {
						input:
							samplename      				= samplename,
							submission_id   				= submission_id,
							collection_date 				= collection_date,
							sequence        				= sequence,
							iso_host         				= iso_host,
							iso_country 				    = iso_country,
							gisaid_submitter 				= gisaid_submitter,
							iso_state        				= iso_state,
							iso_continent 					= iso_continent,
							iso_county							= iso_county,
							seq_platform   					= seq_platform,
							assembly_method	= assembly_method,
							originating_lab 				= originating_lab,
							origLab_address  				= origLab_address,
							submitting_lab  				= submitting_lab,
							subLab_address 					= subLab_address,
							Authors          				= Authors,

							passage_details = passage_details,
							gender = gender,
							patient_age = patient_age,
							patient_status = patient_status,
							specimen_source = specimen_source,
							outbreak = outbreak,
							last_vaccinated = last_vaccinated,
							treatment = treatment

					}
					call submission.genbank {
						input:
							samplename      = samplename,
							submission_id   = submission_id,
							collection_date = collection_date,
							sequence        = sequence,
							organism        = organism,
							iso_org	        = iso_org,
							iso_host        = iso_host,
							iso_country     = iso_country,
							specimen_source = specimen_source,
							BioProject      = BioProject
					}
				}
			}
		}
	}

	output {
	    File?     read1_submission   = sra.read1_submission
	    File?     read2_submission   = sra.read2_submission
	    File?     SE_read_submission = sra.SE_read_submission
	    File      deID_assembly      = deidentify.deID_assembly
	    File?     genbank_assembly   = genbank.genbank_assembly
	    File?     genbank_metadata   = genbank.genbank_metadata
	    File?     gisaid_assembly    = gisaid.gisaid_assembly
	    File?     gisaid_metadata    = gisaid.gisaid_metadata
	}
}


#coverage >= coverage_gisaid && number_N <= number_N_gisaid &&
