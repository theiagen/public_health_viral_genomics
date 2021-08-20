version 1.0

task ncbi_prep_one_sample {
  input {
    #required files
    File assembly_fasta
    File read1_dehosted
    File read2_dehosted
    
    #required metadata
    String authors
    String bioproject_accession
    String biosample_accession = "{populate_with_bioSample_accession}"
    String collecting_lab
    String collecting_lab_address
    String collection_date
    String continent
    String country
    String gisaid_submitter
    String host_disease
    String host
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
    
    #runtime
    String docker_image = "theiagen/utility:1.1"
    Int  mem_size_gb = 1
    Int CPUs = 1
    Int disk_size = 25
    Int preemptible_tries = 0
  }
  command <<<
    #Check date format
    if [[ ~{collection_date} != [0-3][0-9]/[0-1][0-9]/[0-9][0-9][0-9][0-9] ]]
    then 
      echo "Incorrect collection date format; collection date must be in  YYYY-MM-DD format." 1>&2
      exit 1
    else
      year=$(echo ~{collection_date} | cut -f 1 -d '-')
    fi
       
    #Format GenBank metadata and assembly
    isolate="~{organism}/~{host}/~{country}/~{submission_id}/${year}"
    
    ##GenBank assembly
    ###removing leading Ns, folding sequencing to 75 bp wide, and adding metadata for genbank submissions
    echo ">~{submission_id}" > ~{submission_id}.genbank.fa
    grep -v ">" ~{assembly_fasta} | sed 's/^N*N//g' | fold -w 75 >> ~{submission_id}_genbank.fasta
    
    ##GenBank modifier
    echo -e "Sequence_ID\tcountry\thost\tisolate\tcollection-date\tisolation-source\tBioSample\tBioProject\tnote" > ~{submission_id}_genbank_modifier.tsv
    echo -e "~{submission_id}\t~{country}\t~{host_sci_name}\t${isolate}\t~{collection_date}\t~{isolation_source}\t~{biosample_accession}\t{bioproject_accession} >> ~{submission_id}_genbank_modifier.tsv"
    
    #Format GISAID metadata and assembly

    ##GISAID assembly
    echo ">~{submission_id}" > ~{submission_id}.genbank.fa
    grep -v ">" ~{assembly_fasta} | sed 's/^N*N//g' | fold -w 75 >> ~{submission_id}_genbank.fasta
    
    ##GISAID modifier
    echo -e "Sequence_ID\tcountry\thost\tisolate\tcollection-date\tisolation-source\tBioSample\tBioProject\tnote" > ~{submission_id}_genbank_modifier.tsv
    echo -e "~{submission_id}\t~{country}\t~{host_sci_name}\t${isolate}\t~{collection_date}\t~{isolation_source}\t~{biosample_accession}\t{bioproject_accession} >> ~{submission_id}_genbank_modifier.tsv"
  >>>

  output {
    File biosample_attributes = "~{submission_id}_biosample_attributes.tsv"
    File sra_metadata = "~{submission_id}_sra_metadata.tsv"
    File genbank_assembly = "~{submission_id}_genbank.fasta"
    File genbank_modifier = "~{submission_id}_genbank_modifier.tsv"
  }

  runtime {
      docker:       "theiagen/utility:1.1"
      memory:       "~{mem_size_gb} GB"
      cpu:          CPUs
      disks:        "local-disk ~{disk_size} SSD"
      preemptible:  preemptible_tries
      maxRetries:   3
  }
}
task gisaid_prep_one_sample {
  input {
    #required files
    File assembly_fasta
   
    #required metadata
    String authors
    String assembly_method
    String collecting_lab
    String collecting_lab_address
    String collection_date
    String continent
    String assembly_mean_coverage
    String country
    String gisaid_submitter
    String host
    String organism = "hCoV-19"
    String seq_platform
    String state
    String submission_id
    String submitting_lab
    String submitting_lab_address
    String type="betacoronavirus"
    
    #optional metadata
    String? county
    String? gender
    String? last_vaccinated
    String? passage_details
    String? patient_age
    String? patient_status
    String? purpose_of_sequencing
    String? outbreak
    String? specimen_source
    String? treatment
    
    #runtime
    String docker_image = "theiagen/utility:1.1"
    Int  mem_size_gb = 1
    Int CPUs = 1
    Int disk_size = 25
    Int preemptible_tries = 0
  }
  command <<<
    #Check date format
    if [[ ~{collection_date} != [0-3][0-9]/[0-1][0-9]/[0-9][0-9][0-9][0-9] ]]
    then 
      echo "Incorrect collection date format; collection date must be in  YYYY-MM-DD format." 1>&2
      exit 1
    else
      year=$(echo ~{collection_date} | cut -f 1 -d '-')
    fi 
    
    #Format GISAID metadata and assembly
    ##GISAID assembly
    echo ">~{organism}/${country}/${submission_id}/$year" > ${submission_id}.gisaid.fa
    grep -v ">" ~{assembly_fasta} >> ${submission_id}_gisaid.fasta
    
    ##GISAID tMetadata
    echo "submitter,fn,covv_virus_name,covv_type,covv_passage,covv_collection_date,covv_location,covv_add_location,covv_host,covv_add_host_info,covv_sampling_strategy,covv_gender,covv_patient_age,covv_patient_status,covv_specimen,covv_outbreak,covv_last_vaccinated,covv_treatment,covv_seq_technology,covv_assembly_method,covv_coverage,covv_orig_lab,covv_orig_lab_addr,covv_provider_sample_id,covv_subm_lab,covv_subm_lab_addr,covv_subm_sample_id,covv_authors,covv_comment,comment_type" >  ${submission_id}.gisaidMeta.csv

    echo "Submitter,FASTA filename,Virus name,Type,Passage details/history,Collection date,Location,Additional location information,Host,Additional host information,Sampling Strategy,Gender,Patient age,Patient status,Specimen source,Outbreak,Last vaccinated,Treatment,Sequencing technology,Assembly method,Coverage,Originating lab,Address,Sample ID given by the sample provider,Submitting lab,Address,Sample ID given by the submitting laboratory,Authors,Comment,Comment Icon" >> ${submission_id}.gisaidMeta.csv

    echo "\"~{gisaid_submitter}\",\"~{submission_id}.gisaid.fa\",\"~{organism}/~{country}/~{submission_id}/$year\",\"~{type}\",\"~{passage_details}\",\"~{collection_date}\",\"~{continent}/~{country}/~{state}/~{county}\",,\"~{host}\",,\"~{purpose_of_sequencing}\",\"~{gender}\",\"~{patient_age}\",\"~{patient_status}\",\"~{specimen_source}\",\"~{outbreak}\",\"~{last_vaccinated}\",\"~{treatment}\",\"~{seq_platform}\",\"~{assembly_method}\",\"~{assembly_mean_coverage}\",\"~{collecting_lab}\",\"~{collecting_lab_address}\",,\"~{submitting_lab}\",\"~{submitting_lab_address}\",,\"~{authors}\",," >> ${submission_id}.gisaidMeta.csv 

  >>>

  output {
    File gisaid_assembly = "~{submission_id}_gisaid.fasta"
    File gisaid_metadata = "~{submission_id}_gisaid_metadata.tsv"
  }

  runtime {
      docker:       "~{docker_image}"
      memory:       "~{mem_size_gb} GB"
      cpu:          CPUs
      disks:        "local-disk ~{disk_size} SSD"
      preemptible:  preemptible_tries
      maxRetries:   3
  }
}


task compile {
input {
  Array[File] single_submission_fasta
  Array[File] single_submission_meta
  Array[Int]  vadr_num_alerts
  Int         vadr_threshold=0
  String      repository
  String      docker_image = "theiagen/utility:1.1"
  Int         mem_size_gb = 1
  Int         CPUs = 1
  Int         disk_size = 25
  Int         preemptible_tries = 0
}

  command <<<

  assembly_array=(~{sep=' ' single_submission_fasta})
  meta_array=(~{sep=' ' single_submission_meta})
  vadr_array=(~{sep=' ' vadr_num_alerts})

  # remove samples that excede vadr threshold
  for index in ${!assembly_array[@]}; do
    assembly=${assembly_array[$index]}
    meta=${meta_array[$index]}
    vadr=${vadr_array[$index]}

    if [ "${vadr}" -gt "~{vadr_threshold}" ]; then
      assembly_array=( "${assembly_array[@]/$assembly}" )
      meta_array=( "${meta_array[@]/$meta}" )
    fi
    done


  head -n -1 ${meta_array[1]} > ~{repository}_upload_meta.csv
  for i in ${meta_array[*]}; do
      echo $i
      tail -n1 $i >> ~{repository}_upload_meta.csv
  done

  cat ${assembly_array[*]} > ~{repository}_upload.fasta

  >>>

  output {
    File    upload_meta   = "${repository}_upload_meta.csv"
    File    upload_fasta  = "${repository}_upload.fasta"

  }

  runtime {
      docker:       docker_image
      memory:       "~{mem_size_gb} GB"
      cpu:          CPUs
      disks:        "local-disk ~{disk_size} SSD"
      preemptible:  preemptible_tries
      maxRetries:   3
  }
}
