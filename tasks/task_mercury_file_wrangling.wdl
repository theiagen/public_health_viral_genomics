version 1.0

task sm_metadata_wrangling { # the sm stands for supermassive
  input {
    String table_name
    String workspace_name
    String project_name
    File? input_table
    Array[String] sample_names
    String organism = "sars-cov-2"
    String output_name
    String gcp_bucket_uri
    Int vadr_alert_limit = 0 # only for SC2
  }
  command <<<
    # when running on terra, comment out all input_table mentions
    python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --entity_type ~{table_name} --tsv_filename ~{table_name}-data.tsv
    
    # when running locally, use the input_table in place of downloading from Terra
    #cp ~{input_table} ~{table_name}-data.tsv

    echo "DEBUG: Now entering Python block to perform parsing of metadata"

    python3 <<CODE 
    import pandas as pd 
    import numpy as np 
    import os 
    import sys

    # set a function to grab only the year from the date
    def year_getter(x):
      for name, values in x.iteritems():
        for item in values:
          if item != [0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]:
            print("Incorrect collection date format; collection date must be in YYYY-MM-DD format.")
            sys.exit(1)  
          else:
            return item.split("-")[0]

    # set a function to remove NA values and return the cleaned table and a table of excluded samples
    def remove_nas(table, required_metadata):
      # remove rows with blank cells from table -- figure out how to specify which row was blank
      table.replace(r'^\s+$', np.nan, regex=True) # replace blank cells with NaNs 
      excluded_samples = table[table[required_metadata].isna().any(axis=1)] # write out all rows that are required with NaNs to a new table
      excluded_samples = excluded_samples[excluded_samples.columns.intersection(required_metadata)] # remove all optional rows so only required rows are shown
      excluded_samples = excluded_samples.loc[:, excluded_samples.isna().any()] # remove all NON-NA columns so only columns with NAs remain; Shelly is a wizard and I love her 
      table.dropna(subset=required_metadata, axis=0, how='any', inplace=True) # remove all rows that are required with NaNs from table

      return table, excluded_samples

    # read exported Terra table into pandas
    tablename = ~{table_name}-data.tsv 
    table = pd.read_csv(tablename, delimiter='\t', header=0)

    # extract the samples for upload from the entire table
    table = table[table["~{table_name}_id"].isin("~{sep='*' sample_names}".split("*"))]


    # make some standard variables that are used multiple times
    table["year"] = df.apply(year_getter(table["collection_date"]))
    table["isolate"] = (table["organism"] + "/" + table["host"] + "/" + table["submission_id"] + "/" + table["year"])
    table["biosample_accession"] = "{populate_with_BioSample_accession}"

    # set required and optional metadata fields based on the organism type
    if ("~{organism}" == "sars-cov-2"):
      print("Organism is SARS-CoV-2; performing VADR check")
      # perform vadr alert check
      table.drop(table.index[table["vadr_num_alerts"] > ~{vadr_alert_limit}], inplace=True)

      # set default values
      table["gisaid_organism"] = "hCoV-19"
      table["gisaid_virus_name"] = (table["organism"] + "/" + table["country"] + "/" + table["submission_id"] + "/" + table["year"])

      # set required and optional metadata fields
      biosample_required = ["submission_id", "bioproject_accession", "organism", "collecting_lab", "collection_date", "country", "state", "host_sci_name", "host_disease", "isolate", "isolation_source"]
      biosample_optional = ["treatment", "gisaid_accession", "gisaid_virus_name", "patient_age", "patient_gender", "purpose_of_sampling", "purpose_of_sequencing"]
  
      sra_required = ["bioproject_accession", "submission_id", "library_ID", "organism", "isolation_source", "library_strategy", "library_source", "library_selection", "library_layout", "seq_platform", "instrument_model", "design_description", "filetype", "read1_dehosted"]
      sra_optional = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "assembly_method", "dehosting_method", "submitter_email"]

      genbank_required = ["submission_id", "country", "host_sci_name", "isolate", "collection_date", "isolation_source", "biosample_accession", "bioproject_accession", "assembly_fasta"]

      gisaid_required = ["gisaid_submitter", "gisaid_virus_name", "submission_id", "type", "passage_details", "collection_date", "continent", "country", "state", "host", "seq_platform", "assembly_fasta", "assembly_method", "assembly_mean_coverage", "collecting_lab", "collecting_lab_address", "submitting_lab", "submitting_lab_address", "authors"]
      gisaid_optional = ["county", "purpose_of_sequencing", "patient_gender", "patient_age", "patient_status", "specimen_source", "outbreak", "last_vaccinated", "treatment"]

      required_metadata = biosample_required + sra_required + genbank_required + gisaid_required
      table, excluded_samples = remove_nas(table, required_metadata)
      excluded_samples.to_csv("~{output_name}_excluded_samples.tsv", sep='\t', index=False)

      # SC2 BIOSAMPLE
      print("DEBUG: creating biosample metadata table...")
      biosample_metadata = table[[biosample_required]].copy()
      for column in biosample_optional:
        if column in table.columns:
          biosample_metadata[column] = table[column]
        else:
          biosample_metadata[column] = ""
      
      biosample_metadata["geo_loc_name"] = biosample_metadata["country"] + ": " + biosample_metadata["state"]
      biosample_metadata.drop(["country", "state"], axis=1)
      biosample_metadata.rename(columns={"submission_id" : "sample_name", "host_sci_name" : "host", "treatment" : "antiviral_treatment_agent", "patient_gender" : "host_sex", "patient_age" : "host_age", "collecting_lab" : "collected_by"}, inplace=True)

      biosample_metadata.to_csv("~{output_name}_biosample_metadata.tsv", sep='\t', index=False)

      # SC2 SRA
      print("DEBUG: creating sra metadata table...")
      sra_metadata = table[[sra_required]].copy()
      for column in sra_optional:
        if column in table.columns:
          sra_metadata[column] = table[column]
        else: # add the column
          sra_metadata[column] = ""
      sra_metadata.rename(columns={"submission_id" : "sample_name", "seq_platform" : "platform", "amplicon_primer_scheme" : "amplicon_PCR_primer_scheme", "assembly_method" : "raw_sequence_data_processing_method", "submitter_email" : "sequence_submitter_contact_email"}, inplace=True)
      sra_metadata["biosample_accession"] = "{populate with BioSample accession}"
      sra_metadata["title"] = "Genomic sequencing of " + sra_metadata["organism"] + ": " + sra_metadata["isolation_source"]
      sra_metadata.drop(["organism", "isolation_source"], axis=1)

      # prettify the filenames and rename them to be sra compatible; write out copy commands to a file to rename and move later
      sra_metadata["filename"] = sra_metadata["sample_name"] + "_R1.fastq.gz"
      sra_metadata["copy_command_r1"] = "gsutil -m cp " + sra_metadata["read1"] + " " + ~{gcp_bucket_uri} + "/" + sra_metadata["filename"]
      sra_metadata["copy_command_r1"].to_csv("sra-file-transfer.sh", index=False, header=False)
      sra_metadata.drop(["copy_command_r1", "read1"], axis=1)
      if "read2_dehosted" in table.columns: # enable optional single end submission
        sra_metadata["filename2"] = sra_metadata["sample_name"] + "_R2.fastq.gz"
        sra_metadata["copy_command_r2"] = "gsutil -m cp " + sra_metadata["read2"] + " " + ~{gcp_bucket_uri} + "/" + sra_metadata["filename2"]
        sra_metadata["copy_command_r2"].to_csv("sra-file-transfer.sh", mode='a', index=False, header=False)
        sra_metadata.drop(["copy_command_r2", "read2"], axis=1)

      sra_metadata.to_csv("~{output_name}_sra_metadata.tsv", sep='\t', index=False)

      # GENBANK
      print("DEBUG: creating genbank metadata table...")
      genbank_metadata = table[[genbank_required]].copy()
      genbank_metadata.rename(columns={"submission_id" : "Sequence_ID", "host_sci_name" : "host", "collection_date" : "collection-date", "isolation_source" : "isolation-source", "biosample_accession" : "BioSample", "bioproject_accession" : "BioProject"}, inplace=True)

      # prep for file manipulation and manuevering 
      genbank_metadata["cp"] = "gsutil cp"
      genbank_metadata["fn"] = genbank_metadata["Sequence_ID"] + "_genbank_untrimmed.fasta"
      genbank_metadata.to_csv("genbank-file-transfer.sh", sep=' ', header=False, index=False, columns = ["cp", "assembly_fasta", "fn"])
      genbank_metadata.drop(["cp", "assembly_fasta"], axis=1)

      # replace the first line of every fasta file (>Sample_ID) with the gisaid virus name instead (>covv_virus_name)
      # since gisaid virus name includes '/', use '|' in sed command instead
      genbank_metadata["rename_fasta_header"] = "sed -i '1s|.*|>" + genbank_metadata["Sequence_ID"] + "|' " + genbank_metadata["fn"]
      genbank_metadata.to_csv("genbank-fasta-manipulation.sh", header=False, index=False, columns = ["rename_fasta_header"])
      genbank_metadata.drop("rename_fasta_header", "fn")

      genbank_metadata.to_csv("~{output_name}_genbank_metadata.tsv", sep='\t', index=False)

      # SC2 GISAID
      print("DEBUG: creating gisaid metadata file...")
      gisaid_metadata = table[[gisaid_required]].copy() 
      for column in gisaid_optional:                  
        if column in table.columns:                     
          gisaid_metadata[column] = table[column]  
        else:
          gisaid_metadata[column] = "" 

      # create gisaid-specific variables; drop variables that are not included in gisaid metadata
      gisaid_metadata["covv_location"] = (gisaid_metadata["continent"] + " / " + gisaid_metadata["country"] + " / " + gisaid_metadata["state"]) 
      gisaid_metadata.drop(["continent", "country", "state"], axis=1)
      if gisaid_metadata["county"]:
        gisaid_metadata["covv_location"] = (gisaid_metadata["covv_location"] + " / " + gisaid_optional["county"])
        gisaid_metadata.drop("county", axis=1)
      gisaid_metadata["covv_type"] = "betacoronavirus"

      # make new column for filename
      gisaid_metadata["fn"] = gisaid_metadata[submission_id] + "_gisaid.fasta"
      gisaid_metadata.drop("submission_id", axis=1)

      # write out the command to rename the assembly files to a file for bash to move about
      gisaid_metadata["cp"] = "gsutil cp"
      gisaid_metadata.to_csv("gisaid-file-transfer.sh", sep=' ', header=False, index=False, columns = ["cp", "assembly_fasta", "fn"])
      gisaid_metadata.drop("cp", axis=1)

      # replace the first line of every fasta file (>Sample_ID) with the gisaid virus name instead (>covv_virus_name)
      # since gisaid virus name includes '/', use '|' in sed command instead
      gisaid_metadata["rename_fasta_header"] = "sed -i '1s|.*|>" + gisaid_metadata["gisaid_virus_name"] + "|' " + gisaid_metadata["fn"]
      gisaid_metadata.to_csv("gisaid-fasta-manipulation.sh", header=False, index=False, columns = ["rename_fasta_header"])
      gisaid_metadata.drop("rename_fasta_header")

      # make dictionary for renaming headers
      # format: {original : new} or {metadata_formatter : gisaid_format}
      gisaid_rename_headers = {"gisaid_virus_name" : "covv_virus_name", "gisaid_submitter" : "submitter", "passage_details" : "covv_passage", "collection_date" : "covv_collection_date", "seq_platform" : "covv_seq_technology", "host" : "covv_host", "assembly_method" : "covv_assembly_method", "assembly_mean_coverage" : "covv_coverage", "collecting_lab" : "covv_orig_lab", "collecting_lab_address" : "covv_orig_lab_addr", "submitting_lab" : "pos_subm_lab", "submitting_lab_address" : "covv_subm_lab_addr", "authors" : "covv_authors", "purpose_of_sequencing" : "covv_sampling_strategy", "patient_gender" : "covv_gender", "patient_age" : "covv_patient_age", "patient_status" : "covv_patient_status", "specimen_source" : "covv_specimen", "outbreak" : "covv_outbreak", "last_vaccinated" : "covv_last_vaccinated", "treatment" : "covv_treatment"}
      
      # rename columns
      gisaid_metadata.rename(columns=gisaid_rename_headers, inplace=True)

      # write out to gisaid output metadata file
      ###### THIS DOES NOT INCLUDE THE SECONDARY HEADER LINE
      gisaid_metadata.to_csv("~{output_name}_gisaid_metadata.csv", sep=',', index=False)

    elif ("~{organism}" == "mpox"):
      print("Organism is mpox, no VADR filtering performed")
      table["gisaid_organism"] = "mpx/A"
      table["gisaid_virus_name"] = (table["organism"] + "/" + table["country"] + "/" + table["submission_id"] + "/" + table["year"])

      biosample_required = ["submission_id", "organism", "collected_by", "collection_date",  "country", "state", "host", "host_disease", "isolation_source", "lat_lon", "isolation_type"]
      biosample_optional = ["sample_title", "bioproject_accession", "strain", "isolate", "culture_collection", "genotype", "host_age", "host_description", "host_disease_outcome", "host_disease_stage", "host_health_state", "host_sex", "host_subject_id", "host_tissue_sampled", "passage_history", "pathotype", "serotype", "serovar", "specimen_voucher", "subgroup", "subtype", "description"] 

      sra_required = ["bioproject_accession", "submission_id", "library_id", "organism", "isolation_source", "library_strategy", "library_source", "library_selection", "library_layout", "seq_platform", "instrument_model", "design_description", "filetype", "read1_dehosted"]
      sra_optional = ["read2_dehosted", "amplicon_primer_scheme", "amplicon_size", "assembly_method", "dehosting_method", "submitter_email"]

      # ADD BIOPROJECT_ACCESSION TO FASTA HEADER (?)
      bankit_required = ["submission_id", "isolate", "collection_date", "country", "host", "assembly_fasta"]
      bankit_optional = ["isolation_source", "bioproject_accession"]

      gisaid_required = ["gisaid_submitter", "gisaid_virus_name", "submission_id", "passage_details", "collection_date", "continent", "country", "state", "host", "seq_platform", "assembly_fasta", "assembly_method", "assembly_mean_coverage", "collecting_lab", "collecting_lab_address", "submitting_lab", "submitting_lab_address", "authors"]
      gisaid_optional = ["county", "purpose_of_sequencing", "patient_gender", "patient_age", "patient_status", "specimen_source", "outbreak", "last_vaccinated", "treatment"]

      print("DEBUG: removing rows with NAs in required columns...")
      # remove all rows with NAs in required columns and capture which rows and columns have those NAs
      required_metadata = biosample_required + sra_required + bankit_requried + gisaid_required
      table, excluded_samples = remove_nas(table, required_metadata)
      excluded_samples.to_csv("~{output_name}_excluded_samples.tsv", sep='\t', index=False)
     
      # BIOSAMPLE
      print("DEBUG: creating biosample metadata table...")
      biosample_metadata = table[[biosample_required]].copy()
      for column in biosample_optional:
        if column in table.columns:
          biosample_metadata[column] = table[column]
        else:
          biosample_metadata[column] = ""
      biosample_metadata.rename(columns={"submission_id" : "sample_name"}, inplace=True)
      biosample_metadata["geo_loc_name"] = biosample_metadata["country"] + ": " + biosample_metadata["state"]
      biosample_metadata.drop(["country", "state"], axis=1)

      biosample_metadata.to_csv("~{output_name}_biosample_metadata.tsv", sep='\t', index=False)

      # SRA
      print("DEBUG: creating sra metadata table...")
      sra_metadata = table[[sra_required]].copy()
      for column in sra_optional:
        if column in table.columns:
          sra_metadata[column] = table[column]
        else: # add the column
          sra_metadata[column] = ""
      sra_metadata.rename(columns={"submission_id" : "sample_name", "seq_platform" : "platform", "amplicon_primer_scheme" : "amplicon_PCR_primer_scheme", "assembly_method" : "raw_sequence_data_processing_method", "submitter_email" : "sequence_submitter_contact_email"}, inplace=True)
      sra_metadata["biosample_accession"] = "{populate with BioSample accession}"
      sra_metadata["title"] = "Genomic sequencing of " + sra_metadata["organism"] + ": " + sra_metadata["isolation_source"]
      sra_metadata.drop(["organism", "isolation_source"], axis=1)
       
      # prettify the filenames and rename them to be sra compatible
      sra_metadata["filename"] = sra_metadata["sample_name"] + "_R1.fastq.gz"
      sra_metadata["copy_command_r1"] = "gsutil -m cp " + sra_metadata["read1"] + " " + ~{gcp_bucket_uri} + "/" + sra_metadata["filename"]
      sra_metadata["copy_command_r1"].to_csv("sra-file-transfer.sh", index=False, header=False)
      sra_metadata.drop(["copy_command_r1", "read1"], axis=1)
      if "read2_dehosted" in table.columns: # enable optional single end submission
        sra_metadata["filename2"] = sra_metadata["sample_name"] + "_R2.fastq.gz"
        sra_metadata["copy_command_r2"] = "gsutil -m cp " + sra_metadata["read2"] + " " + ~{gcp_bucket_uri} + "/" + sra_metadata["filename2"]
        sra_metadata["copy_command_r2"].to_csv("sra-file-transfer.sh", mode='a', index=False, header=False)
        sra_metadata.drop(["copy_command_r2", "read2"], axis=1)

      sra_metadata.to_csv("~{output_name}_sra_metadata.tsv", sep='\t', index=False)

      # BANKIT
      print("DEBUG: creating bankit sqn file...")
      bankit_metadata = table[[bankit_required]].copy()
      for column in bankit_optional:
        if column in table.columns:
          bankit_metadata[column] = table[column]
        else:
          bankit_metadata[column] = ""
      
      bankit_metadata.rename({"submission_id" : "Sequence_ID", "isolate" : "Isolate", "collection_date" : "Collection_date", "country" : "Country", "host" : "Host", "isolation_source" : "Isolation_source"}, inplace=True)
      
      bankit_metadata["cp"] = "gsutil cp"
      bankit_metadata["fn"] = bankit_metadata["Sequence_ID"] + "_bankit.fasta"
      bankit_metadata.to_csv("bankit-file-transfer.sh", sep=' ', header=False, index=False, columns = ["cp", "assembly_fasta", "fn"])
      bankit_metadata.drop(["cp", "assembly_fasta"], axis=1)

      # replace the first line of every fasta file (>Sample_ID) with the gisaid virus name instead (>covv_virus_name)
      # since gisaid virus name includes '/', use '|' in sed command instead
      bankit_metadata["rename_fasta_header"] = "sed -i '1s|.*|>" + bankit_metadata["Sequence_ID"] + "|' " + bankit_metadata["fn"]
      bankit_metadata.to_csv("bankit-fasta-manipulation.sh", header=False, index=False, columns = ["rename_fasta_header"])
      bankit_metadata.drop("rename_fasta_header", "fn")

      bankit_metadata.to_csv("~{output_name}.src", sep='\t', index=False)


      # GISAID
      print("DEBUG: creating gisaid metadata file...")
      gisaid_metadata = table[[gisaid_required]].copy() 
      for column in gisaid_optional:                  
        if column in table.columns:                     
          gisaid_metadata[column] = table[column]  
        else:
          gisaid_metadata[column] = "" 

      # create gisaid-specific variables; drop variables that are not included in gisaid metadata
      gisaid_metadata["pox_location"] = (gisaid_metadata["continent"] + " / " + gisaid_metadata["country"] + " / " + gisaid_metadata["state"]) 
      gisaid_metadata.drop(["continent", "country", "state"], axis=1)
      if gisaid_metadata["county"]:
        gisaid_metadata["pox_location"] = (gisaid_metadata["pox_location"] + " / " + gisaid_optional["county"])
        gisaid_metadata.drop("county", axis=1)

      # make new column for filename
      gisaid_metadata["fn"] = gisaid_metadata[submission_id] + "_gisaid.fasta"
      gisaid_metadata.drop("submission_id", axis=1)

      # write out the command to rename the assembly files to a file for bash to move about
      gisaid_metadata["cp"] = "gsutil cp"
      gisaid_metadata.to_csv("gisaid-file-transfer.sh", sep=' ', header=False, index=False, columns = ["cp", "assembly_fasta", "fn"])
      gisaid_metadata.drop(["cp", "assembly_fasta"], axis=1)

      # replace the first line of every fasta file (>Sample_ID) with the gisaid virus name instead (>pox_virus_name)
      # since gisaid virus name includes '/', use '|' in sed command instead
      gisaid_metadata["rename_fasta_header"] = "sed -i '1s|.*|>" + gisaid_metadata["gisaid_virus_name"] + "|' " + gisaid_metadata["fn"]
      gisaid_metadata.to_csv("gisaid-fasta-manipulation.sh", header=False, index=False, columns = ["rename_fasta_header"])
      gisaid_metadata.drop("rename_fasta_header")

      # make dictionary for renaming headers
      # format: {original : new} or {metadata_formatter : gisaid_format}
      gisaid_rename_headers = {"gisaid_virus_name" : "pox_virus_name", "gisaid_submitter" : "submitter", "passage_details" : "pox_passage", "collection_date" : "pox_collection_date", "seq_platform" : "pox_seq_technology", "host" : "pox_host", "assembly_method" : "pox_assembly_method", "assembly_mean_coverage" : "pox_coverage", "collecting_lab" : "pox_orig_lab", "collecting_lab_address" : "pox_orig_lab_addr", "submitting_lab" : "pos_subm_lab", "submitting_lab_address" : "pox_subm_lab_addr", "authors" : "pox_authors", "purpose_of_sequencing" : "pox_sampling_strategy", "patient_gender" : "pox_gender", "patient_age" : "pox_patient_age", "patient_status" : "pox_patient_status", "specimen_source" : "pox_specimen_source", "outbreak" : "pox_outbreak", "last_vaccinated" : "pox_last_vaccinated", "treatment" : "pox_treatment"}
      
      # rename columns
      gisaid_metadata.rename(columns=gisaid_rename_headers, inplace=True)

      # write out to gisaid output metadata file
      ###### THIS DOES NOT INCLUDE THE SECONDARY HEADER LINE
      gisaid_metadata.to_csv("~{output_name}_gisaid_metadata.csv", sep=',', index=False)

      #gisaid_out.write("submitter,fn,pox_virus_name,pox_passage,pox_collection_date,pox_location,pox_host,pox_sampling_strategy,pox_gender,pox_patient_age,pox_patient_status,pox_specimen,pox_outbreak,pox_last_vaccinated,pox_treatment,pox_seq_technology,pox_assembly_method,pox_coverage,pox_orig_lab,pox_orig_lab_addr,pox_subm_lab,pox_subm_lab_addr,pox_authors")
      
    else:
      raise Exception('Only "SARS-CoV-2" and "MPXV" are supported as acceptable input for the \'organism\' variable at this time. You entered "~{organism}".')
  
    CODE

    echo "DEBUG: performing file transfers and manipulations"
    # this version of gsutil only works on python2.7
    export CLOUDSDK_PYTHON=python2.7

    # transfer gisaid files, alter header lines, then concatenate all gisaid fasta files
    bash gisaid-file-transfer.sh
    bash gisaid-fasta-manipulation.sh
    cat *_gisaid.fasta > ~{output_name}_gisaid.fasta

    if [[ ~{organism} == "sars-cov-2" ]]; then
      bash genbank-file-transfer.sh
      bash genbank-fasta-manipulation.sh
      cat *_genbank_untrimmed.fasta > ~{output_name}_genbank_untrimmed.fasta
    fi

    if [[ ~{organism} == "mpox" ]] ; then
      bash bankit-file-transfer.sh
      bash bankit-fasta-manipulation.sh
      cat *_bankit.fasta > ~{output_name}.fsa
    fi

    # transfer sra files to gcp bucket
    bash sra-file-transfer.sh    

    unset CLOUDSDK_PYTHON   # reset env var


  >>>
  output {
    File excluded_samples = "~{output_name}_excluded_samples.tsv"
    File biosample_metadata = "~{output_name}_biosample_metadata.tsv"
    File sra_metadata = "~{output_name}_sra_metadata.tsv"
    File? genbank_metadata = "~{output_name}_genbank_metadata.tsv"
    File? genbank_untrimmed_fasta = "~{output_name}_genbank_untrimmed.fasta"
    File? bankit_metadata = "~{output_name}.src"
    File? bankit_fasta = "~{output_name}.fsa"
    File gisaid_metadata = "~{output_name}_gisaid_metadata.csv"
    File gisaid_fasta = "~{output_name}_gisaid.fasta"
  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task trim_genbank_fastas {
  input {
    File genbank_untrimmed_fasta
    String output_name
    Int minlen = 50
    Int maxlen = 30000
  }
  command <<<
    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      ~{genbank_untrimmed_fasta} \
      --minlen ~{minlen} \
      --maxlen ~{maxlen} \
      > ~{output_name}_genbank.fasta
  >>>
  output {
    File genbank_fasta = "~{output_name}_genbank.fasta"
  }
  runtime {
    docker: "quay.io/staphb/vadr:1.3"
    memory: "1 GB"
    cpu: 1    
    disks: "local-disk 25 SSD"
    preemptible: 0
    maxRetries: 3
  }
}


## I think this works, but honestly not sure.
task table2asn {
  input {
    File authors_sbt # have users provide the .sbt file for MPXV submission-- it can be created here: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
    File bankit_fasta
    File bankit_metadata
    String output_name
  }
  command <<<
    # using this echo statement so the fasta file doesn't have a wiggly line
    echo "~{bankit_fasta} file needs to be localized for the program to access"

    # rename authors_sbt to contain output_name
    mv ~{authors_sbt} ~{output_name}.sbt

    # convert the data into a sqn file so it can be emailed to NCBI
    table2asn -t ~{output_name}.sbt \
      -src-file ~{bankit_metadata} \
      -indir . \
      -a s # inputting a set of fasta data

  >>>
  output {
    File sqn_file = "~{output_name}.sqn"
  }
  runtime {
    docker: "staphb/ncbi-table2asn:1.26.678"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 25 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
