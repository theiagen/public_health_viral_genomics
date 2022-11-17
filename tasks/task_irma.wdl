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
	Int memory = 8
	Int cpu = 4
  }
  command <<<
    date | tee DATE
    #capture reads as bash variables
    read1=~{read1}
    read2=~{read2}
    # capture irma vesion
    IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION
    # set config if needed
    if ~{keep_ref_deletions}; then 
      touch irma_config.sh
      echo 'DEL_TYPE="NNN"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi
    # format reads from sra
    if ~{from_sra} ; then
      echo "SRA reads will be formatted to meet IRMA input requirements"
      sra_id=$(echo "~{read_basename}" | awk -F "_" '{ print $1 }')
      zcat ~{read1} | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 1:1" : $0}' | gzip -c > "${sra_id}-irmafix_R1.fastq.gz"
      zcat ~{read2} | awk '{print (NR%4 == 1) ? "@'${sra_id}'-" ++i " 2:2" : $0}' | gzip -c > "${sra_id}-irmafix_R2.fastq.gz"
      #modify read variables
      read1="${sra_id}-irmafix_R1.fastq.gz"
      read2="${sra_id}-irmafix_R2.fastq.gz"  
    fi
    # run IRMA 
    IRMA "~{irma_module}" "${read1}" "${read2}" ~{samplename}
    # cat consensus assemblies
    if [ -d "~{samplename}/amended_consensus/" ]; then
      cat ~{samplename}/amended_consensus/*.fa > ~{samplename}.irma.consensus.fasta
    fi
  >>>
  output {
    File irma_assembly_fasta = "~{samplename}.irma.consensus.fasta"
    File seg4_ha_assembly = "~{samplename}/amended_consensus/~{samplename}_4.fa"
    String irma_version = read_string("VERSION")
    String irma_pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible:  0
  }
}