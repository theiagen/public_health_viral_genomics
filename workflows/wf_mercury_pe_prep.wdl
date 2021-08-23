version 1.0

import "../tasks/task_versioning.wdl" as versioning
import "../tasks/task_pub_repo_prep.wdl" as submission_prep

workflow mercury_pe_prep {
  input {
    #required files
    File assembly_fasta
    File read1_dehosted
    File read2_dehosted
    
    #required metadata (titan gc outputs)
    String assembly_method
    String assembly_mean_coverage
    
    #required metadata (user inputs)
    String authors
    String bioproject_accession
    String collecting_lab
    String collecting_lab_address
    String collection_date
    String continent
    String country
    String dehosting_method = "NCBI Human Scrubber"
    String gisaid_submitter
    String gisaid_organism = "hCoV-19"
    String filetype = "fastq"
    String host ="Human"
    String host_disease
    String host_sci_name = "Homo sapien"
    String instrument_model
    String isolation_source
    String library_id
    String library_layout = "paired"
    String library_selection
    String library_source
    String library_strategy
    String organism
    Int number_N
    String seq_platform
    String state
    String submission_id
    String submitting_lab
    String submitting_lab_address  
    
    #optional metadata
    String? amplicon_primer_scheme
    String? amplicon_size
    String? biosample_accession
    String? gisaid_accession
    String? county
    String? patient_gender
    String? patient_age
    String? purpose_of_sampling
    String? purpose_of_sequencing
    String? submitter_email
    String? treatment

    # Optional user-defined thresholds for generating submission files
    Int number_N_threshold = 5000
  }
  
  if (number_N <= number_N_threshold) {
    call submission_prep.ncbi_prep_one_sample {
      input: 
        amplicon_primer_scheme = amplicon_primer_scheme,
        amplicon_size = amplicon_size,
        assembly_fasta = assembly_fasta,
        assembly_method = assembly_method,
        bioproject_accession = bioproject_accession,
        biosample_accession = biosample_accession,
        collecting_lab = collecting_lab,
        collection_date = collection_date,
        country = country,
        dehosting_method = dehosting_method,
        design_description = "Whole genome sequencing of ~{organism}",
        filetype = filetype,
        gisaid_accession = gisaid_accession,
        gisaid_organism = gisaid_organism,
        host = host,
        host_disease = host_disease,
        host_sci_name = host_sci_name,
        instrument_model = instrument_model,
        isolation_source = isolation_source,
        library_id = library_id,
        library_layout = library_layout,
        library_selection = library_selection,
        library_source = library_source,
        library_strategy = library_strategy,
        organism = organism,
        patient_age = patient_age,
        patient_gender = patient_gender,
        purpose_of_sampling = purpose_of_sampling,
        purpose_of_sequencing = purpose_of_sequencing,
        read1_dehosted = read1_dehosted,
        read2_dehosted = read2_dehosted,
        seq_platform = seq_platform,
        state = state,
        submission_id = submission_id,
        submitter_email = submitter_email,
        treatment = treatment	 
    }
    call submission_prep.gisaid_prep_one_sample {
      input: 
        assembly_fasta = assembly_fasta,
        authors = authors,
        assembly_method = assembly_method,
        collecting_lab = collecting_lab,
        collecting_lab_address = collecting_lab_address,
        collection_date = collection_date,
        continent = continent,
        assembly_mean_coverage = assembly_mean_coverage,
        country = country,
        gisaid_submitter = gisaid_submitter,
        host = host,
        organism = gisaid_organism,
        seq_platform = seq_platform,
        state = state,
        submission_id = submission_id,
        submitting_lab = submitting_lab,
        submitting_lab_address = submitting_lab_address,
        county = county,
        patient_gender = patient_gender,
        patient_age = patient_age,
        purpose_of_sequencing = purpose_of_sequencing,
        treatment = treatment
    }
  }

  call versioning.version_capture{
    input:
  }
  output {
    String mercury_pe_prep_version = version_capture.phvg_version
    String mercury_pe_prep_analysis_date = version_capture.date
    
    File? biosample_attributes = ncbi_prep_one_sample.biosample_attributes
    File? sra_metadata = ncbi_prep_one_sample.sra_metadata
    File? genbank_assembly = ncbi_prep_one_sample.genbank_assembly
    File? genbank_modifier = ncbi_prep_one_sample.genbank_modifier
    File? sra_read1 = ncbi_prep_one_sample.sra_read1
    File? sra_read2 = ncbi_prep_one_sample.sra_read2
    Array[File]? sra_reads = ncbi_prep_one_sample.sra_reads
    
    File? gisaid_assembly = gisaid_prep_one_sample.gisaid_assembly
    File? gisaid_metadata = gisaid_prep_one_sample.gisaid_metadata
  }
}


#coverage >= coverage_gisaid && number_N <= number_N_gisaid &&
