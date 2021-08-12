version 1.0


task demultiplexing {

  input {
    Array[File] basecalled_reads
    String?     run_prefix="artic_ncov2019"
    Int?        normalise=200
    Int?        cpu=8
  }

  command{
    guppy_barcoder -t \$cpus --require_barcodes_both_ends -i . -s . --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg  barcode_arrs_nb96.cfg" -q 0 -r

  }

  output {
    Array[File]       demultiplexed_reads = glob("*.fq.gz")
  }

  runtime {
    docker:       "genomicpariscentre/guppy"
    memory:       "16 GB"
    cpu:          8
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task read_filtering {

  input {
    File        demultiplexed_reads
    String      samplename
    String?     run_prefix="artic_ncov2019"
    Int?        min_length=400
    Int?        max_length=700
    Int?        cpu=8
  }

  command{
    # date and version control
    mkdir ~{samplename}
    cp ~{demultiplexed_reads} ~{samplename}/
    echo "DIRNAME: $(dirname)"
    artic guppyplex --min-length ${min_length} --max-length ${max_length} --directory ~{samplename} --prefix ${run_prefix}

  }

  output {
    File       filtered_reads = "${run_prefix}_~{samplename}.fastq"
  }

  runtime {

    docker:       "staphb/artic-ncov2019:1.3.0"
    memory:       "16 GB"
    cpu:          8
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task consensus {
  ## Need to output multiple directories
  input {
    String  samplename
    File    filtered_reads
    File    primer_bed
    Int?    normalise=20000
    Int?    cpu=8
    String  medaka_model="r941_min_high_g360"
    String  docker="staphb/artic-ncov2019:1.3.0"
  }
  String primer_name = basename(primer_bed)
  
  command{
    # setup custom primer scheme (/V is required by Artic)
    mkdir -p ./primer-schemes/nCoV-2019/Vuser
    cp /primer-schemes/nCoV-2019/V3/nCoV-2019.reference.fasta ./primer-schemes/nCoV-2019/Vuser/nCoV-2019.reference.fasta
    cp ${primer_bed} ./primer-schemes/nCoV-2019/Vuser/nCoV-2019.scheme.bed

    # version control
    echo "Medaka via $(artic -v)" | tee VERSION
    echo "${primer_name}" | tee PRIMER_NAME
    artic minion --medaka --medaka-mode ~{medaka_model} --normalise ~{normalise} --threads ~{cpu} --scheme-directory ./primer-schemes --read-file ~{filtered_reads} nCoV-2019/Vuser ~{samplename}

    # clean up fasta header
    echo ">~{samplename}" > ~{samplename}.medaka.consensus.fasta
    grep -v ">" ~{samplename}.consensus.fasta >> ~{samplename}.medaka.consensus.fasta
  }

  output {
    File    consensus_seq = "~{samplename}.medaka.consensus.fasta"
    File    sorted_bam = "~{samplename}.trimmed.rg.sorted.bam"
    File    trim_sorted_bam = "~{samplename}.primertrimmed.rg.sorted.bam"
    File    trim_sorted_bai = "~{samplename}.primertrimmed.rg.sorted.bam.bai"
    File    medaka_pass_vcf = "~{samplename}.pass.vcf.gz" 
    String  artic_pipeline_version = read_string("VERSION")
    String  primer_bed_name = read_string("PRIMER_NAME")
  }

  runtime {
    docker:       "~{docker}"
    memory:       "16 GB"
    cpu:          8
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
