version 1.0

import "../tasks/task_se_pub_repo_submission.wdl" as submission

workflow mercury_batch {
	input {
		Array[File] genbank_single_submission_fasta
		Array[File] genbank_single_submission_meta
		Array[File] gisaid_single_submission_fasta
		Array[File] gisaid_single_submission_meta
		Array[String] samplename
		Array[Int] vadr_num_alerts

	}

	call submission.compile as genbank_compile{
		input:
			single_submission_fasta=genbank_single_submission_fasta,
	    		single_submission_meta=genbank_single_submission_meta,
					samplename=samplename,
					vadr_num_alerts=vadr_num_alerts,
	    		repository="GenBank"
	}
	call submission.compile as gisaid_compile{
		input:
			single_submission_fasta=gisaid_single_submission_fasta,
	  	single_submission_meta=gisaid_single_submission_meta,
			samplename=samplename,
			vadr_num_alerts=vadr_num_alerts,
	    repository="GISAID"
	}


	output {
	    File?      GenBank_upload_meta  = genbank_compile.upload_meta
	    File?      GenBank_upload_fasta = genbank_compile.upload_fasta
			File			 GenBank_batched_samples = genbank_compile.batched_samples
			File			 GenBank_excluded_samples = genbank_compile.excluded_samples
	    File?      GISAID_upload_meta  = gisaid_compile.upload_meta
	    File?      GISAID_upload_fasta = gisaid_compile.upload_fasta
			File			 GISAID_batched_samples = gisaid_compile.batched_samples
			File			 GISAID_excluded_samples = gisaid_compile.excluded_samples

	}
}
