version 1.0

task supermassive_file_wrangling {
  input {
    String table_name
    String workspace_name
    String project_name
    File? input_table
    Array[String] sample_names
    String organism = "SARS-CoV-2"

  }
  command <<<
    # when running on terra, comment out all input_table mentions
    python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --entity_type ~{table_name} --tsv_filename ~{table_name}-data.tsv
    
    # when running locally, use the input_table in place of downloading from Terra
    #cp ~{input_table} ~{table_name}-data.tsv

    python3 <<CODE 
    import pandas as pd 
    import numpy as np 
    import os 

    # read export table into pandas
    tablename = ~{table_name}-data.tsv 
    table = pd.read_csv(tablename, delimiter='\t', header=0)

    # extract the samples for upload from the entire table
    table = table[table["~{table_name}_id"].isin("~{sep='*' sample_names}".split("*"))]

    # set required and optional metadata fields based on the organism type

    ## work out what to do with the optional biosample_accession and gisaid_organism variable as those currently have default values.
    if ("~{organism}" == "SARS-CoV-2"):
      required_metadata = ["assembly_fasta", "read1_dehosted", "assembly_method", "bioproject_accession", "collecting_lab", "collection_date", "country", "design_description", "dehosting_method", "filetype", "host_disease", "host", "host_sci_name", "instrument_model", "isolation_source", "library_id", "library_layout", "library_selection", "library_source", "library_strategy", "organism", "seq_platform", "state", "submission_id"]
      optional_metadata = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "biosample_accession", "gisaid_accession",  "gisaid_organism", "patient_age", "patient_gender", "purpose_of_sampling", "purpose_of_sequencing", "submitter_email", "treatment"]

    elif ("~{organism}" == "MPXV"):
      ## TO DO: need to change these to include MPXV-required fields
      required_metadata = ["assembly_fasta", "read1_dehosted", "assembly_method", "bioproject_accession", "collecting_lab", "collection_date", "country", "design_description", "dehosting_method", "filetype", "host_disease", "host", "host_sci_name", "instrument_model", "isolation_source", "library_id", "library_layout", "library_selection", "library_source", "library_strategy", "organism", "seq_platform", "state", "submission_id"]
      optional_metadata = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "biosample_accession", "gisaid_accession",  "gisaid_organism", "patient_age", "patient_gender", "purpose_of_sampling", "purpose_of_sequencing", "submitter_email", "treatment"]

    else:
      raise Exception('Only "SARS-CoV-2" and "MPXV" are supported as acceptable input for the \'organism\' variable at this time. You entered "~{organism}".')
    


    # sra metadata fields:
    

    # remove rows with blank cells from table -- figure out how to specify which row was blank
    table.replace(r'^\s+$', np.nan, regex=True) # replace blank cells with NaNs 
    excluded_samples = table[table[required_metadata].isna().any(axis=1)] # write out all rows that are required with NaNs to a new table
    excluded_samples["~{table_name}_id"].to_csv("excluded_samples.tsv", sep='\t', index=False, header=False) # write the excluded names out to a file
    table.dropna(subset=required_metadata, axis=0, how='any', inplace=True) # remove all rows that are required with NaNs from table
    
    # extract the required metadata from the table




  >>>
  output {

  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}