version 1.0

task sm_metadata_wrangling {
  input {
    String table_name
    String workspace_name
    String project_name
    File? input_table
    Array[String] sample_names
    String organism = "SARS-CoV-2"
    String output_name

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

      sra_required = ["bioproject_accession", "submission_id", "library_ID", "organism", "isolation_source", "library_strategy", "library_source", "library_selection", "library_layout", "seq_platform", "instrument_model", "design_description", "filetype", "read1_dehosted"]
      sra_optional = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "assembly_method", "dehosting_method", "submitter_email"]

      genbank_required = ["submission_id", "country", "host_sci_name", "isolate", "collection_date", "isolation_source", "biosample_accession", "bioproject_accession"]

      # TYPE IS BETACORONAVIRUS
      # ORGANISM = hCoV-19  
      gisaid_required = ["gisaid_submitter", "organism", "country", "submission_id", "year", "type", "passage_details", "collection_date", "continent", "country", "state", "host", "seq_platform", "assembly_method", "assembly_mean_coverage", "collecting_lab", "collecting_lab_address", "submitting_lab", "submitting_lab_address", "authors"]
      gisaid_optional = ["county", "purpose_of_sequencing", "patient_gender", "patient_age", "patient_status", "specimen_source", "outbreak", "last_vaccinated", "treatment"]


      required_metadata = biosample_required + sra_required + genbank_required + gisaid_required

    elif ("~{organism}" == "MPXV"):
      ## TO DO: need to change these to include MPXV-required fields
      #  required_metadata = ["assembly_fasta", "read1_dehosted", "assembly_method", "bioproject_accession", "collecting_lab", "collection_date", "country", "design_description", "dehosting_method", "filetype", "host_disease", "host", "host_sci_name", "instrument_model", "isolation_source", "library_id", "library_layout", "library_selection", "library_source", "library_strategy", "organism", "seq_platform", "state", "submission_id"]
      #  optional_metadata = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "biosample_accession", "gisaid_accession",  "gisaid_organism", "patient_age", "patient_gender", "purpose_of_sampling", "purpose_of_sequencing", "submitter_email", "treatment"]

      biosample_required = ["submission_id", "organism", "collected_by", "collection_date", "geo_loc_name", "host", "host_disease", "isolation_source", "lat_lon", "isolation_type"]
      biosample_optional = ["sample_title", "bioproject_accession", "attribute_package", "strain", "isolate", "culture_collection", "genotype", "host_age", "host_description", "host_disease_outcome", "host_disease_stage", "host_health_state", "host_sex", "host_subject_id", "host_tissue_sampled", "passage_history", "pathotype", "serotype", "serovar", "specimen_voucher", "subgroup", "subtype", "description"] 

      sra_required = ["bioproject_accession", "submission_id", "library_id", "organism", "isolation_source", "library_strategy", "library_source", "library_selection", "library_layout", "seq_platform", "instrument_model", "design_description", "filetype", "read1_dehosted"]
      sra_optional = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "assembly_method", "dehosting_method", "submitter_email"]

      # ADD BIOPROJECT_ACCESSION TO FASTA HEADER
      bankit_required = ["submission_id", "isolate", "collection_date", "country", "host"]
      bankit_optional = ["isolation_source", "passage_details"]

      gisaid_required = ["gisaid_submitter", "organism", "country", "submission_id", "year", "passage_details", "collection_date", "continent", "country", "state", "host", "seq_platform", "assembly_method", "assembly_mean_coverage", "collecting_lab", "collecting_lab_address", "submitting_lab", "submitting_lab_address", "authors"]
      gisaid_optional = ["county", "purpose_of_sequencing", "patient_gender", "patient_age", "patient_status", "specimen_source", "outbreak", "last_vaccinated", "treatment"]


      required_metadata = biosample_required + sra_required + bankit_requried + gisaid_required



    else:
      raise Exception('Only "SARS-CoV-2" and "MPXV" are supported as acceptable input for the \'organism\' variable at this time. You entered "~{organism}".')
    

    # remove rows with blank cells from table -- figure out how to specify which row was blank
    table.replace(r'^\s+$', np.nan, regex=True) # replace blank cells with NaNs 
    excluded_samples = table[table[required_metadata].isna().any(axis=1)] # write out all rows that are required with NaNs to a new table
    table.dropna(subset=required_metadata, axis=0, how='any', inplace=True) # remove all rows that are required with NaNs from table

    # for row in excluded_samples, compare column names to required_metadata
    # append missing column names to list
    # print sample_name and then list of missing column names
   
    # excluded_samples["~{table_name}_id"].to_csv("excluded_samples.tsv", sep='\t', index=False, header=False) # write the excluded names out to a file
   
    if ("~{organism}" == "SARS-CoV-2"):
      # remove rows that have > 0 vadr alerts
      table.drop(table.index[table["vadr_num_alerts"] > 0], inplace=True)
    else:
      print("Organism is MPXV, no VADR filtering performed")

    if ("~{organism}" == "SARS-CoV-2"):

      
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
     biosample_metadata = table[[biosample_required]].copy()
      for column in biosample_optional:
        if column in table.columns:
          biosample_metadata[column] = table[column]
      biosample_metadata.rename(columns={"submission_id" : "sample_name", "" : ""}, inplace=True)

      sra_metadata = table[[sra_required]].copy()
      for column in sra_optional:
        if column in table.columns:
          sra_metadata[column] = table[column]
        else: # add the column
      sra_metadata.rename(columns={"submission_id" : "sample_name"}, inplace=True)

      bankit_metadata = table[[bankit_required]].copy()
      for column in bankit_optional:
        if column in table.columns:
          bankit_metadata[column] = table[column]
      
      
      # gisaid combined metadata
      gisaid_metadata = table[[gisaid_required]].copy() # all required with old header
      for column in gisaid_optional:                    # do any optional things exist in original table?
        if column in table.columns:                     # if so,
          gisaid_metadata[column] = table[column]       # add to gisaid table with old header
        else:
          gisaid_metadata[column] = ""                 # add the optional column even though it's blank for easy renaming



      # now have full table with all possible columns
      # then we can rename any old header
      # using fancy dictionary

      # gisaid output file 
      gisaid_out = open("~{output_name}_gisaid_metadata.csv", "w")
      gisaid_out.write("submitter,fn,pox_virus_name,pox_passage,pox_collection_date,pox_location,pox_add_location,pox_host,pox_add_host_info,pox_sampling_strategy,pox_gender,pox_patient_age,pox_patient_status,pox_specimen,pox_outbreak,pox_last_vaccinated,pox_treatment,pox_seq_technology,pox_assembly_method,pox_coverage,pox_orig_lab,pox_orig_lab_addr,pox_provider_sample_id,pox_subm_lab,pox_subm_lab_addr,pox_subm_sample_id,pox_authors")
      gisaid_out.write("Submitter,FASTA filename,Virus name,Passage details/history,Collection date,Location,Additional location information,Host,Additional host information,Sampling Strategy,Gender,Patient age,Patient status,Specimen source,Outbreak,Last vaccinated,Treatment,Sequencing technology,Assembly method,Depth of coverage,Originating lab,Address,Sample ID given by the sample provider,Submitting lab,Address,Sample ID given by the submitting laboratory,Authors")
      # how to ensure each column gets written in the correct order??

      # parse date for year
      def year_getter(x):
        for name, values in x.iteritems():
          for item in values:
            if item != [0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]:
              return  
              # the date is formatted incorrectly
              # deal with this later
            else:
              return item.split("-")[0]

      # create gisaid-specific variables
      gisaid_metadata["year"] = df.apply(year_getter(gisaid_metadata["collection_date"]))
      gisaid_metadata["pox_virus_name"] = (gisaid_metadata["organism"] + "/" + gisaid_metadata["country"] + "/" + gisaid_metadata["submission_id"] + "/" + gisaid_metadata["year"])
      gisaid_metadata["pox_location"] = (gisaid_metadata["continent"] + " / " + gisaid_metadata["country"] + " / " + gisaid_metadata["state"]) 
      gisaid_metadata.drop("year", "continent", "country", "state"], axis=1)
      if gisaid_metadata["county"]:
        gisaid_metadata["pox_location"] = (gisaid_metadata["pox_location"] + " / " + gisaid_optional["county"])
        gisaid_metadata.drop(["county"], axis=1)


      
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

      >> # END PYTHON CODE

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