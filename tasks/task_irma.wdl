version 1.0

task irma {
  input {
    File read1
    File read2
    String samplename
    Boolean keep_ref_deletions = true
    String irma_module = "FLU-utr"
    Boolean from_sra = false
    String read_basename = basename(read1)
    String docker = "quay.io/biocontainers/irma:1.0.2--pl5321hdfd78af_2"
  }
  command <<<
    date | tee DATE
    # capture irma vesion
    IRMA | head -n1 | awk -F' ' '{ print $5 }' | tee VERSION
    # set config if needed
    if [ ~{keep_ref_deletions} ]; then 
      touch irma_config.sh
      echo 'DEL_TYPE="NNN"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi
    # format reads from sra
    if ~{from_sra} ; then
      sra_id=$(echo "~{read_basename}" | awk -F "_" '{ print $1 }')
      zcat ~{read1} | awk '{print (NR%4 == 1) ? "@${sra_id}-" ++i " 1:1" : $0}' | gzip -c > "${sra_id}-irmafix_R1.fastq.gz"
      zcat ~{read2} | awk '{print (NR%4 == 1) ? "@${sra_id}-" ++i " 2:2" : $0}' | gzip -c > "${sra_id}-irmafix_R2.fastq.gz"
      #if from sra, i don't know if IRMA run below will read the original read1 or the fixed reads
      #whcih of these would work?
      #read1=zcat ~{read1} | awk '{print (NR%4 == 1) ? "@${sra_id}-" ++i " 1:1" : $0}' | gzip -c > "${sra_id}-irmafix_R1.fastq.gz"
      #read1="${sra_id}-irmafix_R1.fastq.gz"            
    fi
    # run IRMA 
    IRMA ~{irma_module} ~{read1} ~{read2} ~{samplename}

    # cat consensus assemblies
    if [ -d "~{samplename}/amended_consensus/" ]; then
      cat ~{samplename}/amended_consensus/*.fa > ~{samplename}.irma.consensus.fasta
    fi
  >>>
  output {
    File irma_consensus = "~{samplename}.irma.consensus.fasta"
    String irma_version = read_string("VERSION")
    String irma_pipeline_date = read_string("DATE")
    File? read1_irmafix = "~{samplename}-irmafix_R1.fastq.gz"
    File? read2_irmafix = "~{samplename}-irmafix_R2.fastq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible:  0
  }
}