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
      # BIOSAMPLE_ACCESSION = {populate with biosample accession}
      # ISOLATE = ORGANISM/HOST/COUNTRY/SUBMISSION_ID/YEAR
      # GISAID_VIRUS_NAME = hCoV-19/COUNTRY/SUBMISSION_ID/YEAR
      biosample_required = ["submission_id", "bioproject_accession", "organism", "collecting_lab", "collection_date", "country", "state", "host_sci_name", "host_disease", "isolate", "isolation_source"]
      biosample_optional = ["treatment", "gisaid_accession", "gisaid_virus_name", "patient_age", "patient_gender", "purpose_of_sampling", "purpose_of_sequencing"]

      sra_required = ["bioproject_accession", "submission_id", "library_id", "organism", "isolation_source", "library_strategy", "library_source", "library_selection", "library_layout", "seq_platform", "instrument_model", "design_description", "filetype", "read1_dehosted"]
      sra_optional = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "assembly_method", "dehosting_method", "submitter_email"]

      genbank_required = ["submission_id", "country", "host_sci_name", "isolate", "collection_date", "isolation_source", "biosample_accession", "bioproject_accession"]

      # TYPE IS BETACORONAVIRUS
      # ORGANISM = hCoV-19  
      gisaid_required = ["gisaid_submitter", "organism", "country", "submission_id", "year", "type", "passage_details", "collection_date", "continent", "country", "state", "host", "seq_platform", "assembly_method", "assembly_mean_coverage", "collecting_lab", "collecting_lab_address", "submitting_lab", "submitting_lab_address", "authors"]
      gisaid_optional = ["county", "purpose_of_sequencing", "patient_gender", "patient_age", "patient_status", "specimen_source", "outbreak", "last_vaccinated", "treatment"]


    elif ("~{organism}" == "MPXV"):
      ## TO DO: need to change these to include MPXV-required fields
      required_metadata = ["assembly_fasta", "read1_dehosted", "assembly_method", "bioproject_accession", "collecting_lab", "collection_date", "country", "design_description", "dehosting_method", "filetype", "host_disease", "host", "host_sci_name", "instrument_model", "isolation_source", "library_id", "library_layout", "library_selection", "library_source", "library_strategy", "organism", "seq_platform", "state", "submission_id"]
      optional_metadata = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "biosample_accession", "gisaid_accession",  "gisaid_organism", "patient_age", "patient_gender", "purpose_of_sampling", "purpose_of_sequencing", "submitter_email", "treatment"]

      biosample_required = ["submission_id", "organism", "collected_by", "collection_date", "geo_loc_name", "host", "host_disease", "isolation_source", "lat_lon", "isolation_type"]
      biosample_optional = ["sample_title", "bioproject_accession", "attribute_package", "strain", "isolate", "culture_collection", "genotype", "host_age", "host_description", "host_disease_outcome", "host_disease_stage", "host_health_state", "host_sex", "host_subject_id", "host_tissue_sampled", "passage_history", "pathotype", "serotype", "serovar", "specimen_voucher", "subgroup", "subtype", "description"] 

      sra_required = ["bioproject_accession", "submission_id", "library_id", "organism", "isolation_source", "library_strategy", "library_source", "library_selection", "library_layout", "seq_platform", "instrument_model", "design_description", "filetype", "read1_dehosted"]
      sra_optional = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "assembly_method", "dehosting_method", "submitter_email"]

      # ADD BIOPROJECT_ACCESSION TO FASTA HEADER
      bankit_required = ["submission_id", "isolate", "collection_date", "country", "host"]
      bankit_optional = ["isolation_source", "passage_details"]

      gisaid_required = ["gisaid_submitter", "organism", "country", "submission_id", "year", "passage_details", "collection_date", "continent", "country", "state", "host", "seq_platform", "assembly_method", "assembly_mean_coverage", "collecting_lab", "collecting_lab_address", "submitting_lab", "submitting_lab_address", "authors"]
      gisaid_optional = ["county", "purpose_of_sequencing", "patient_gender", "patient_age", "patient_status", "specimen_source", "outbreak", "last_vaccinated", "treatment"]

      


    else:
      raise Exception('Only "SARS-CoV-2" and "MPXV" are supported as acceptable input for the \'organism\' variable at this time. You entered "~{organism}".')
    


    # sra metadata fields:
    

    # remove rows with blank cells from table -- figure out how to specify which row was blank
    table.replace(r'^\s+$', np.nan, regex=True) # replace blank cells with NaNs 
    excluded_samples = table[table[required_metadata].isna().any(axis=1)] # write out all rows that are required with NaNs to a new table
    excluded_samples["~{table_name}_id"].to_csv("excluded_samples.tsv", sep='\t', index=False, header=False) # write the excluded names out to a file
    table.dropna(subset=required_metadata, axis=0, how='any', inplace=True) # remove all rows that are required with NaNs from table

    # remove rows that have > 0 vadr alerts
    table.drop(table.index[table["vadr_num_alerts"] > 0], inplace=True)
    
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