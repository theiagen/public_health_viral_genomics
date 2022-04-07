version 1.0

workflow freyja_update {
  input {
    String gcp_uri
  }
  call freyja_update_refs {
    input:
      gcp_uri = gcp_uri
  }
  call transfer_files {
    input:
      updated_barcodes = freyja_update_refs.updated_barcodes,
      updated_lineages = freyja_update_refs.updated_lineages,
      update_log = freyja_update_refs.update_log,
      gcp_uri = gcp_uri
  }
  output {
  }
}
task freyja_update_refs {
  input {
    String gcp_uri
    String docker = "staphb/freyja:1.3.4"
  }
  command <<<
  # Create updated refrence files
  mkdir freyja_update_refs && cd freyja_update_refs
  echo "PWD: $PWD"
  freyja update --outdir $PWD
  
  find / -name "lineagePaths.txt"
  echo "Freyja reference files created using the freyja update command; Freyja Docker Image utilized: ~{docker}. More information can be found at https://github.com/andersen-lab/Freyja" > $PWD/update_log.txt
  >>>
  runtime {
    memory: "4 GB"
    cpu: 2
    docker: "~{docker}"
    disks: "local-disk 100 HDD"
  }
  output {
    File updated_barcodes = "freyja_update_refs/usher_barcodes.csv"
    File updated_lineages = "freyja_update_refs/curated_lineages.json"
    File update_log = "freyja_update_refs/update_log.txt"
  }
}
task transfer_files {
  input {
    String gcp_uri
    File updated_barcodes
    File updated_lineages
    File update_log
  }
  command <<<
  # transfer_files to specified gcp_uri
  date_tag=$(date +"%Y-%m-%d")
  
  gsutil -m cp ~{updated_barcodes} ~{updated_lineages} ~{update_log} ~{gcp_uri}/${date_tag}
  
  >>>
  runtime {
    memory: "4 GB"
    cpu: 2
    docker: "theiagen/utility:1.1"
    disks: "local-disk 100 HDD"
  }
  output {    
  }
}