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
  >>>

  output {
    File biosample_attributes = "~{submission_id}_biosample_attributes.tsv"
    File sra_metadata = "~{submission_id}_sra_metadata.tsv"
    File genbank_assembly = "~{submission_id}_genbank.fasta"
    File genbank_modifier = "~{submission_id}_genbank_modifier.tsv"
    File gisaid_assembly = "~{submission_id}_gisaid.fasta"
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
