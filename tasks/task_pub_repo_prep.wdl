version 1.0

task prep_one_sample {
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
    String host_sci_name
    String isolate
    String organism
    Int number_N
    String state
    String submission_id
    String submitting_lab
    String submitting_lab_address
    
    #optional metadata
    String? gender
    String? patient_age
    String? county
    String? specimen_processing
    String? patient_age_unit
    String? patient_age_bin
    String? purpose_of_sampling
    String? purpose_of_sampling_details
    String? purpose_of_sequencing
    String? sequencing_protocol_name
    
    #runtime
    String docker_image = "theiagen/utility:1.1"
    Int  mem_size_gb = 1
    Int CPUs = 1
    Int disk_size = 25
    Int preemptible_tries = 0
  }
  command {
    # de-identified consensus/assembly sequence
    
  }

  output {
    File biosample_attributes = "~{submission_id}_biosample_attributes.tsv"
    File sra_metadata = "~{submission_id}_sra_metadata.tsv"
    File genbank_assembly = "~{submission_id}_genbank_assembly.tsv"
    File genbank_modifier = "~{submission_id}_genbank_modifier.tsv"
    File gisaid_assembly = "~{submission_id}_gisaid_assembly.tsv"
    File gisaid_metadata = "~{submission_id}_gisaid_metadata.tsv"
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
