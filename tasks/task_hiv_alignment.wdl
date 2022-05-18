version 1.0

task bwa {
  input {
    File read1
    File read2
    String samplename
    File reference_genome = "gs://theiagen-public-files/terra/hivgc-files/NC_001802.1.fasta"
    Int cpu = 6
  }
  command <<<
    # date and version control
    date | tee DATE
    echo "BWA $(bwa 2>&1 | grep Version )" | tee BWA_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # set reference genome
    ref_genome="~{reference_genome}"
    bwa index "~{reference_genome}"

    # Map with BWA MEM
    echo "Running bwa mem -t ~{cpu} ${ref_genome} ~{read1} ~{read2} | samtools sort | samtools view -F 4 -o ~{samplename}.sorted.bam "
    bwa mem \
    -t ~{cpu} \
    "${ref_genome}" \
    ~{read1} ~{read2} |\
    samtools sort | samtools view -F 4 -o ~{samplename}.sorted.bam

    # index BAMs
    samtools index ~{samplename}.sorted.bam
  >>>
  output {
    String bwa_version = read_string("BWA_VERSION")
    String sam_version = read_string("SAMTOOLS_VERSION")
    File sorted_bam = "${samplename}.sorted.bam"
    File sorted_bai = "${samplename}.sorted.bam.bai"
  }
  runtime {
    docker: "quay.io/staphb/ivar:1.3.1-titan"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task mafft {
  input {
    Array[File] genomes
    Int cpu = 16
  }
  command <<<
    # date and version control
    date | tee DATE
    mafft_vers=$(mafft --version)
    echo Mafft $(mafft_vers) | tee VERSION

    cat ~{sep=" " genomes} | sed 's/Consensus_//;s/.consensus_threshold.*//' > assemblies.fasta
    mafft --thread -~{cpu} assemblies.fasta > msa.fasta
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File msa = "msa.fasta"
  }
  runtime {
    docker: "quay.io/staphb/mafft:7.450"
    memory: "32 GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
