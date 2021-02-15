version 1.0

task seqyclean {
  input {
    File        read1
    File        read2
    String      samplename
    String?     adapters = "/Adapters_plus_PhiX_174.fasta"
    Int?        seqyclean_minlen = 25
    String?     seqyclean_qual = "20 20"
    Boolean?    compress = true
    Boolean?    seqyclean_dup = false
    Boolean?    seqyclean_no_adapter_trim = false
    Int?        cpus = 16
  }
  
  command {
    # date and version control
    date | tee DATE
    echo "Seqyclean $(seqyclean -h | grep Version)" | tee VERSION

    seqyclean \
    -minlen ${seqyclean_minlen} \
    -qual ${seqyclean_qual} \
    -c ${adapters} \
    ${true="-dup" false="" seqyclean_dup} \
    ${true="-no_adapter_trim" false="" seqyclean_no_adapter_trim} \
    ${true="-gz" false="" compress} \
    -t ${cpus} \
    -1 ${read1} \
    -2 ${read2} \
    -o ${samplename}

    # Capture metrics for summary file
    cut -f 58 ${samplename}_SummaryStatistics.tsv | grep -v "PairsKept" | head -n 1 | tee PAIRS_KEPT
    cut -f 59 ${samplename}_SummaryStatistics.tsv | grep -v "Perc_Kept" | head -n 1 | tee PERCENT_KEPT
  }

  output {
    File       read1_clean = "${samplename}_PE1.fastq.gz"
    File       read2_clean = "${samplename}_PE2.fastq.gz"
    String     version = read_string("VERSION")
    String     pipeline_date = read_string("DATE")
    Int        seqy_pairs = read_string("PAIRS_KEPT")
    Float      seqy_percent = read_string("PERCENT_KEPT")
  }

  runtime {
      docker:       "staphb/seqyclean:1.10.09"
      memory:       "8 GB"
      cpu:          2
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}
