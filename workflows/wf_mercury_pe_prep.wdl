version 1.0

import "../tasks/task_versioning.wdl" as versioning
import "../tasks/task_pub_repo_prep.wdl" as submission_prep

workflow mercury_pe_prep {
  input {
    #required files
    File assembly_fasta
    File read1_dehosted
    File read2_dehosted
    
    #required metadata
    String authors
    String bioproject_accession
    String biosample_accession
    String collecting_lab
    String collecting_lab_address
    String collection_date
    String continent
    String country
    String gisaid_submitter
    String host_disease
    String host_sci_name
    String isolate
    String organism
    Int number_N
    String state
    String submission_id
    String submitting_lab
    String submitting_lab_address
    
    #optional metadata
    String? county
    String? gender
    String? isolation_source
    String? patient_age
    String? patient_age_bin
    String? patient_age_unit
    String? purpose_of_sampling
    String? purpose_of_sampling_details
    String? purpose_of_sequencing
    String? sequencing_protocol_name
    String? specimen_processing

    # Optional user-defined thresholds for generating submission files
    Int number_N_threshold = 5000
  }
  
  if (number_N <= number_N_threshold) {
    call submission_prep.prep_one_sample {
      input:
        assembly_fasta = assembly_fasta,
        read1_dehosted = read1_dehosted,
        read2_dehosted = read2_dehosted,
        authors = authors,
        bioproject_accession = bioproject_accession,
        biosample_accession = biosample_accession,
        collecting_lab = collecting_lab,
        collecting_lab_address = collecting_lab_address,
        collection_date = collection_date,
        continent = continent,
        country = country,
        gisaid_submitter = gisaid_submitter,
        host_disease = host_disease,
        host_sci_name = host_sci_name,
        isolate = isolate,
        organism = organism,
        state = state,
        submission_id = submission_id,
        submitting_lab = submitting_lab,
        submitting_lab_address = submitting_lab_address,
        gender = gender,
        patient_age = patient_age,
        county = county,
        specimen_processing = specimen_processing,
        patient_age_unit = patient_age_unit,
        patient_age_bin = patient_age_bin,
        purpose_of_sampling = purpose_of_sampling,
        purpose_of_sampling_details = purpose_of_sampling_details,
        purpose_of_sequencing = purpose_of_sequencing,
        sequencing_protocol_name = sequencing_protocol_name 
    }
  }

  call versioning.version_capture{
    input:
  }
  output {
    String mercury_pe_prep_version = version_capture.phvg_version
    String mercury_pe_prep_analysis_date = version_capture.date
    
    File? biosample_attributes = prep_one_sample.biosample_attributes
    File? sra_metadata = prep_one_sample.sra_metadata
    File? genbank_assembly = prep_one_sample.genbank_assembly
    File? genbank_modifier = prep_one_sample.genbank_modifier
    File? gisaid_assembly = prep_one_sample.gisaid_assembly
    File? gisaid_metadata = prep_one_sample.gisaid_metadata
  }
}


#coverage >= coverage_gisaid && number_N <= number_N_gisaid &&
