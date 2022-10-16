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



    #######Frank was here
    if ("~{organism}" == "SARS-CoV-2"):
      # remove rows that have > 0 vadr alerts
      table.drop(table.index[table["vadr_num_alerts"] > 0], inplace=True)
    else:
      print("Organism is MPXV, no VADR filtering performed")

    if ("~{organism}" == "SARS-CoV-2"):
      # extract the required metadata from the table
      biosample_required_df = table[[biosample_required]]
      biosample_optional_df = table[[biosample_optional]]
      sra_required_df = table[[sra_required]]
      sra_optional_df = table[[sra_optional]]
      genbank_required_df = table[[genbank_required]]
      gisaid_required_df = table[[gisaid_required]]
      gisaid_optional_df = table[[gisaid_optional]]
      # combine required and optional dfs
      biosample_metadata_df = pd.concat([biosample_required_df, biosample_optional_df], axis=1)
      sra_metadata_df = pd.concat([sra_required_df, sra_optional_df], axis=1)
      genbank_metadata_df = genbank_required_df
      gisaid_metadata_df = pd.concat([gisaid_required_df, gisaid_optional_df], axis=1)
      # remove empty columns
      biosample_metadata_df_clean = biosample_metadata_df.dropna(axis='columns', how='all')
      sra_metadata_df_clean = sra_metadata_df.dropna(axis='columns', how='all')
      bankit_metadata_df_clean = bankit_metadata_df.dropna(axis='columns', how='all')
      gisaid_metadata_df_clean = gisaid_metadata_df.dropna(axis='columns', how='all')
      # determine which columns if any were dropped
        # biosample
      biosample_metadata_df_col_headers = list(biosample_metadata_df.columns)
      biosample_metadata_df_clean_col_headers = list(biosample_metadata_df_clean.columns)
      biosample_blank_col_headers = []
      for i in biosample_metadata_df_col_headers:
          if element not in biosample_metadata_df_clean_col_headers:
              biosample_blank_col_headers.append(i)
        # sra
      sra_metadata_df_col_headers = list(sra_metadata_df.columns)
      sra_metadata_df_clean_col_headers = list(sra_metadata_df_clean.columns)
      sra_blank_col_headers = []
      for i in sra_metadata_df_col_headers:
          if element not in sra_metadata_df_clean_col_headers:
              sra_blank_col_headers.append(i)
        # genbank
      genbank_metadata_df_col_headers = list(genbank_metadata_df.columns)
      genbank_metadata_df_clean_col_headers = list(genbank_metadata_df_clean.columns)
      genbank_blank_col_headers = []
      for i in genbank_metadata_df_col_headers:
          if element not in genbank_metadata_df_clean_col_headers:
              genbank_blank_col_headers.append(i)
        # gisaid
      gisaid_metadata_df_col_headers = list(gisaid_metadata_df.columns)
      gisaid_metadata_df_clean_col_headers = list(gisaid_metadata_df_clean.columns)
      gisaid_blank_col_headers = []
      for i in gisaid_metadata_df_col_headers:
          if element not in gisaid_metadata_df_clean_col_headers:
              gisaid_blank_col_headers.append(i)
      # print metadata files to tsvs
      biosample_metadata_df_clean.to_csv('biosample_metadata.tsv', sep="\t")
      sra_metadata_df_clean.to_csv('sra_metadata.tsv', sep="\t")
      genbank_metadata_df_clean.to_csv('genbank_metadata.tsv', sep="\t")
      gisaid_metadata_df_clean.to_csv('gisaid_metadata.tsv', sep="\t")
      # print blank headers list to tsvs
      biosample_blank_col_headers.to_csv('biosample_blank_col_headers.tsv', sep="\t")
      sra_blank_col_headers.to_csv('sra_blank_col_headers.tsv', sep="\t")
      genbank_blank_col_headers.to_csv('genbank_blank_col_headers.tsv', sep="\t")
      gisaid_blank_col_headers.to_csv('gisaid_blank_col_headers.tsv', sep="\t")
    elif ("~{organism}" == "MPXV"):
       # extract the required metadata from the table
      biosample_required_df = table[[biosample_required]]
      biosample_optional_df = table[[biosample_optional]]
      sra_required_df = table[[sra_required]]
      sra_optional_df = table[[sra_optional]]
      bankit_required_df = table[[bankit_required]]
      bankit_optional_df = table[[bankit_optional]]
      gisaid_required_df = table[[gisaid_required]]
      gisaid_optional_df = table[[gisaid_optional]]
      # combine required and optional dfs
      biosample_metadata_df = pd.concat([biosample_required_df, biosample_optional_df], axis=1)
      sra_metadata_df = pd.concat([sra_required_df, sra_optional_df], axis=1)
      bankit_metadata_df = pd.concat([bankit_required_df, bankit_optional_df], axis=1)
      gisaid_metadata_df = pd.concat([gisaid_required_df, gisaid_optional_df], axis=1)
      # remove empty columns
      biosample_metadata_df_clean = biosample_metadata_df.dropna(axis='columns', how='all')
      sra_metadata_df_clean = sra_metadata_df.dropna(axis='columns', how='all')
      bankit_metadata_df_clean = bankit_metadata_df.dropna(axis='columns', how='all')
      gisaid_metadata_df_clean = gisaid_metadata_df.dropna(axis='columns', how='all')
      # determine which columns if any were dropped
        # biosample
      biosample_metadata_df_col_headers = list(biosample_metadata_df.columns)
      biosample_metadata_df_clean_col_headers = list(biosample_metadata_df_clean.columns)
      biosample_blank_col_headers = []
      for i in biosample_metadata_df_col_headers:
          if element not in biosample_metadata_df_clean_col_headers:
              biosample_blank_col_headers.append(i)
        # sra
      sra_metadata_df_col_headers = list(sra_metadata_df.columns)
      sra_metadata_df_clean_col_headers = list(sra_metadata_df_clean.columns)
      sra_blank_col_headers = []
      for i in sra_metadata_df_col_headers:
          if element not in sra_metadata_df_clean_col_headers:
              sra_blank_col_headers.append(i)
        # bankit
      bankit_metadata_df_col_headers = list(bankit_metadata_df.columns)
      bankit_metadata_df_clean_col_headers = list(bankit_metadata_df_clean.columns)
      bankit_blank_col_headers = []
      for i in bankit_metadata_df_col_headers:
          if element not in bankit_metadata_df_clean_col_headers:
              bankit_blank_col_headers.append(i)
        # gisaid
      gisaid_metadata_df_col_headers = list(gisaid_metadata_df.columns)
      gisaid_metadata_df_clean_col_headers = list(gisaid_metadata_df_clean.columns)
      gisaid_blank_col_headers = []
      for i in gisaid_metadata_df_col_headers:
          if element not in gisaid_metadata_df_clean_col_headers:
              gisaid_blank_col_headers.append(i)
      # print metadata filesto tsvs
      biosample_metadata_df_clean.to_csv('biosample_metadata.tsv', sep="\t")
      sra_metadata_df_clean.to_csv('sra_metadata.tsv', sep="\t")
      bankit_metadata_df_clean.to_csv('genbank_metadata.tsv', sep="\t")
      gisaid_metadata_df_clean.to_csv('gisaid_metadata.tsv', sep="\t")
      # print blank headers list to tsvs
      biosample_blank_col_headers.to_csv('biosample_blank_col_headers.tsv', sep="\t")
      sra_blank_col_headers.to_csv('sra_blank_col_headers.tsv', sep="\t")
      bankit_blank_col_headers.to_csv('bankit_blank_col_headers.tsv', sep="\t")
      gisaid_blank_col_headers.to_csv('gisaid_blank_col_headers.tsv', sep="\t")
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