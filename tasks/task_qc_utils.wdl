version 1.0

task fastqc {

  input {
    File        read1
    File        read2
    String      read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    String      read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
    Int?        cpus = 2
  }

  command {
    # capture date and version
    date | tee DATE
    fastqc --version | grep FastQC | tee VERSION

    fastqc --outdir $PWD --threads ${cpus} ${read1} ${read2}

    unzip -p ${read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
    unzip -p ${read2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

    READ1_SEQS=$(unzip -p ${read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
    READ2_SEQS=$(unzip -p ${read2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

    if [ $READ1_SEQS == $READ2_SEQS ]; then
      read_pairs=$READ1_SEQS
    else
      read_pairs="Uneven pairs: R1=$READ1_SEQS, R2=$READ2_SEQS"
    fi
    echo $read_pairs | tee READ_PAIRS
  }

  output {
    File       fastqc1_html = "${read1_name}_fastqc.html"
    File       fastqc1_zip = "${read1_name}_fastqc.zip"
    File       fastqc2_html = "${read2_name}_fastqc.html"
    File       fastqc2_zip = "${read2_name}_fastqc.zip"
    Int        read1_seq = read_string("READ1_SEQS")
    Int        read2_seq = read_string("READ2_SEQS")
    String        read_pairs = read_string("READ_PAIRS")
    String     version = read_string("VERSION")
    String     pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "quay.io/staphb/fastqc:0.11.9"
    memory:       "4 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
task fastqc_se {

  input {
    File        read1
    String      read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    Int?        cpus = 2
  }

  command {
    # capture date and version
    date | tee DATE
    fastqc --version | grep FastQC | tee VERSION

    fastqc --outdir $PWD --threads ${cpus} ${read1}

    unzip -p ${read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS

    READ_SEQS=$(unzip -p ${read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

    echo $read_pairs | tee READ_PAIRS
  }

  output {
    File       fastqc_html = "${read1_name}_fastqc.html"
    File       fastqc_zip = "${read1_name}_fastqc.zip"
    Int        number_reads = read_string("READ1_SEQS")
    String     version = read_string("VERSION")
    String     pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "quay.io/staphb/fastqc:0.11.8"
    memory:       "4 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
task fastq_scan {

  input {
    File        read1
    File        read2
    String      read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    String      read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
  }

  command <<<
    # capture date and version
    date | tee DATE
    fastq-scan -v | tee VERSION

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then 
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi
    
    # capture forward read stats
    eval "${cat_reads} ~{read1}" | fastq-scan | tee ~{read1_name}_fastq-scan.json >(jq .qc_stats.read_total > READ1_SEQS)
    read1_seqs=$(cat READ1_SEQS)
    eval "${cat_reads} ~{read2}" | fastq-scan | tee ~{read2_name}_fastq-scan.json >(jq .qc_stats.read_total > READ2_SEQS)
    read2_seqs=$(cat READ2_SEQS) 
      
    # capture number of read pairs
    if [ ${read1_seqs} == $read2_seqs]; then
      read_pairs=${read1_seqs}
    else
      read_pairs="Uneven pairs: R1=${read1_seqs}, R2=${read2_seqs}"
    fi

    echo $read_pairs | tee READ_PAIRS
  >>>

  output {
    File       read1_fastq_scan_report = "~{read1_name}_fastq-scan.json"
    File       read2_fastq_scan_report = "~{read2_name}_fastq-scan.json"
    Int        read1_seq = read_string("READ1_SEQS")
    Int        read2_seq = read_string("READ2_SEQS")
    String     read_pairs = read_string("READ_PAIRS")
    String     version = read_string("VERSION")
    String     pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "quay.io/staphb/fastq-scan:0.4.3"
    memory:       "2 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
task fastq_scan_se {

  input {
    File        read1
    String      read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
  }

  command <<<
    # capture date and version
    date | tee DATE
    fastq-scan -v | tee VERSION

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then 
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi
    
    # capture forward read stats
    eval "${cat_reads} ~{read1}" | fastq-scan | tee ~{read1_name}_fastq-scan.json >(jq .qc_stats.read_total > READ1_SEQS)
  

  >>>

  output {
    File       fastq_scan_report = "~{read1_name}_fastq-scan.json"
    Int        read1_seq = read_string("READ1_SEQS")
    String     version = read_string("VERSION")
    String     pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "quay.io/staphb/fastq-scan:0.4.3"
    memory:       "2 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
task consensus_qc {

  input {
    File        assembly_fasta

  }

  command <<<
    # capture date and version
    date | tee DATE
    
    num_N=$( grep -v ">" ~{assembly_fasta} | grep -o 'N' | wc -l )
    if [ -z "$num_N" ] ; then num_N="0" ; fi
    echo $num_N | tee NUM_N

    num_ACTG=$( grep -v ">" ~{assembly_fasta} | grep -o -E "C|A|T|G" | wc -l )
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    echo $num_ACTG | tee NUM_ACTG

    # calculate percent coverage (Wu Han-1 genome length: 29903bp)
    python3 -c "print ( round( ($num_ACTG / 29903 ) * 100, 2 ) )" | tee PERCENT_REF_COVERAGE

    num_degenerate=$( grep -v ">" ~{assembly_fasta} | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
    echo $num_degenerate | tee NUM_DEGENERATE

    num_total=$( grep -v ">" ~{assembly_fasta} | grep -o -E '[A-Z]' | wc -l )
    if [ -z "$num_total" ] ; then num_total="0" ; fi
    echo $num_total | tee NUM_TOTAL
  >>>

  output {
    Int       number_N = read_string("NUM_N")
    Int       number_ATCG = read_string("NUM_ACTG")
    Int       number_Degenerate = read_string("NUM_DEGENERATE")
    Int       number_Total = read_string("NUM_TOTAL")
    Float     percent_reference_coverage = read_string("PERCENT_REF_COVERAGE")
  }

  runtime {
    docker:       "quay.io/theiagen/utility:1.1"    
    memory:       "2 GB"
    cpu:          1
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
