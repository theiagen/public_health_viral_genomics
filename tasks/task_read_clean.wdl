version 1.0

task ncbi_scrub_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "gcr.io/ncbi-sys-gcr-public-research/sra-human-scrubber@sha256:b7dba71079344daea4ea3363e1a67fa54edb7ec65459d039669c68a66d38b140"
    Int disk_size = 100
  }
  String r1_filename = basename(read1)
  String r2_filename = basename(read2)
  command <<<
    # date and version control
    date | tee DATE

    # unzip fwd file as scrub tool does not take in .gz fastq files
    if [[ "~{read1}" == *.gz ]]
    then
      gunzip -c ~{read1} > r1.fastq
      read1_unzip=r1.fastq
    else
      read1_unzip=~{read1}
    fi

    # dehost reads
    /opt/scrubber/scripts/scrub.sh -n ${read1_unzip} |& tail -n1 | awk -F" " '{print $1}' > FWD_SPOTS_REMOVED

    # gzip dehosted reads
    gzip ${read1_unzip}.clean -c > ~{samplename}_R1_dehosted.fastq.gz

    # do the same on read
    # unzip file if necessary
    if [[ "~{read2}" == *.gz ]]
    then
      gunzip -c ~{read2} > r2.fastq
      read2_unzip=r2.fastq
    else
      read2_unzip=~{read2}
    fi

    # dehost reads
    /opt/scrubber/scripts/scrub.sh -n ${read2_unzip} |& tail -n1 | awk -F" " '{print $1}' > REV_SPOTS_REMOVED

    # gzip dehosted reads
    gzip ${read2_unzip}.clean -c > ~{samplename}_R2_dehosted.fastq.gz
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    File read2_dehosted = "~{samplename}_R2_dehosted.fastq.gz"
    Int read1_human_spots_removed = read_int("FWD_SPOTS_REMOVED")
    Int read2_human_spots_removed = read_int("REV_SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
  }
  runtime {
      docker: "~{docker}"
      memory: "8 GB"
      cpu: 4
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 3
  }
}

task ncbi_scrub_se {
  input {
    File read1
    String samplename
    String docker = "gcr.io/ncbi-sys-gcr-public-research/sra-human-scrubber@sha256:b7dba71079344daea4ea3363e1a67fa54edb7ec65459d039669c68a66d38b140"
    Int disk_size = 100
  }
  String r1_filename = basename(read1)
  command <<<
    # date and version control
    date | tee DATE

    # unzip fwd file as scrub tool does not take in .gz fastq files
    if [[ "~{read1}" == *.gz ]]
    then
      gunzip -c ~{read1} > r1.fastq
      read1_unzip=r1.fastq
    else
      read1_unzip=~{read1}
    fi

    # dehost reads
    /opt/scrubber/scripts/scrub.sh -n ${read1_unzip} |& tail -n1 | awk -F" " '{print $1}' > FWD_SPOTS_REMOVED

    # gzip dehosted reads
    gzip ${read1_unzip}.clean -c > ~{samplename}_R1_dehosted.fastq.gz
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    Int read1_human_spots_removed = read_int("FWD_SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task seqyclean {
  input {
    File read1
    File read2
    String samplename
    String? adapters = "/Adapters_plus_PhiX_174.fasta"
    Int? seqyclean_minlen = 15
    String? seqyclean_qual = "20 20"
    Boolean? compress = true
    Boolean? seqyclean_dup = false
    Boolean? seqyclean_no_adapter_trim = false
    Int? cpu = 16
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    echo "Seqyclean $(seqyclean -h | grep Version)" | tee VERSION

    seqyclean \
    -minlen ~{seqyclean_minlen} \
    -qual ~{seqyclean_qual} \
    -c ~{adapters} \
    ~{true="-dup" false="" seqyclean_dup} \
    ~{true="-no_adapter_trim" false="" seqyclean_no_adapter_trim} \
    ~{true="-gz" false="" compress} \
    -t ~{cpu} \
    -1 ~{read1} \
    -2 ~{read2} \
    -o ~{samplename}

    # Capture metrics for summary file
    cut -f 58 ~{samplename}_SummaryStatistics.tsv | grep -v "PairsKept" | head -n 1 | tee PAIRS_KEPT
    cut -f 59 ~{samplename}_SummaryStatistics.tsv | grep -v "Perc_Kept" | head -n 1 | tee PERCENT_KEPT
  >>>
  output {
    File read1_clean = "~{samplename}_PE1.fastq.gz"
    File read2_clean = "~{samplename}_PE2.fastq.gz"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    Int seqy_pairs = read_string("PAIRS_KEPT")
    Float seqy_percent = read_string("PERCENT_KEPT")
  }
  runtime {
      docker: "quay.io/staphb/seqyclean:1.10.09"
      memory: "8 GB"
      cpu: 2
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 3
  }
}

task trimmomatic {
  input {
    File read1
    File read2
    String samplename
    String docker = "quay.io/staphb/trimmomatic:0.39"
    Int? trimmomatic_minlen = 75
    Int? trimmomatic_window_size=4
    Int? trimmomatic_quality_trim_score=30
    Int? threads = 4
    String? trimmomatic_args
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    trimmomatic PE \
    ~{trimmomatic_args} \
    -threads ~{threads} \
    ~{read1} ~{read2} \
    -baseout ~{samplename}.fastq.gz \
    SLIDINGWINDOW:~{trimmomatic_window_size}:~{trimmomatic_quality_trim_score} \
    MINLEN:~{trimmomatic_minlen} &> ~{samplename}.trim.stats.txt

  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File trimmomatic_stats = "~{samplename}.trim.stats.txt"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task trimmomatic_se {
  input {
    File read1
    String samplename
    String docker = "quay.io/staphb/trimmomatic:0.39"
    Int? trimmomatic_minlen = 25
    Int? trimmomatic_window_size = 4
    Int? trimmomatic_quality_trim_score = 30
    Int? threads = 4
    String? trimmomatic_args
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    trimmomatic SE \
    ~{trimmomatic_args} \
    -threads ~{threads} \
    ~{read1} \
    ~{samplename}_trimmed.fastq.gz \
    SLIDINGWINDOW:~{trimmomatic_window_size}:~{trimmomatic_quality_trim_score} \
    MINLEN:~{trimmomatic_minlen} > ~{samplename}.trim.stats.txt
  >>>
  output {
    File read1_trimmed = "~{samplename}_trimmed.fastq.gz"
    File trimmomatic_stats = "~{samplename}.trim.stats.txt"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task fastp {
  input {
    File read1
    File read2
    String samplename
    String docker = "quay.io/staphb/fastp:0.23.2"
    Int? fastp_minlen = 75
    Int? fastp_window_size=4
    Int? fastp_quality_trim_score=30
    # adapter trimming and quality filtering are enabled by default, the flags below disable these functions to match trimmomatic
    String? fastp_args = "--disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --dont_eval_duplication"
    Int? threads = 4
    Int disk_size = 100
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} --in2 ~{read2} \
    --out1 ~{samplename}_1P.fastq.gz --out2 ~{samplename}_2P.fastq.gz \
    --cut_right --cut_window_size ~{fastp_window_size} --cut_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_minlen} \
    --thread ~{threads} \
    ~{fastp_args} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task fastp_se {
  input {
    File read1
    String samplename
    String docker = "quay.io/staphb/fastp:0.23.2"
    Int? fastp_minlen = 75
    Int? fastp_window_size=4
    Int? fastp_quality_trim_score=30
    Int? threads = 4
    # adapter trimming and quality filtering are enabled by default, the flags below disable these functions
    String? fastp_params = "--disable_adapter_trimming --disable_quality_filtering"
    Int disk_size = 100
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} \
    --out1 ~{samplename}_1P.fastq.gz \
    --cut_right --cut_window_size ~{fastp_window_size} --cut_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_minlen} \
    --thread ~{threads} \
    ~{fastp_params} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task bbduk {
  input {
    File read1_trimmed
    File read2_trimmed
    String samplename
    Int memory = 8
    String docker = "quay.io/staphb/bbtools:38.76"
    Int disk_size = 100
    File? adapters
    File? phix
  }
  command <<<
    # date and version control
    date | tee DATE

    # set adapter fasta
    if [[ ! -z "~{adapters}" ]]; then
      echo "Using user supplied FASTA file for adapters..."
      adapter_fasta="~{adapters}"
    else
      echo "User did not supply adapters FASTA file, using default adapters.fa file..."
      adapter_fasta="/bbmap/resources/adapters.fa" 
    fi

    # set phix fasta
    if [[ ! -z "~{phix}" ]]; then
      echo "Using user supplied FASTA file for phiX..."
      phix_fasta="~{phix}"
    else
      echo "User did not supply phiX FASTA file, using default phix174_ill.ref.fa.gz file..."
      phix_fasta="/bbmap/resources/phix174_ill.ref.fa.gz"
    fi

    repair.sh in1=~{read1_trimmed} in2=~{read2_trimmed} out1=~{samplename}.paired_1.fastq.gz out2=~{samplename}.paired_2.fastq.gz

    bbduk.sh in1=~{samplename}.paired_1.fastq.gz in2=~{samplename}.paired_2.fastq.gz out1=~{samplename}.rmadpt_1.fastq.gz out2=~{samplename}.rmadpt_2.fastq.gz ref=${adapter_fasta} stats=~{samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo ordered=t

    bbduk.sh in1=~{samplename}.rmadpt_1.fastq.gz in2=~{samplename}.rmadpt_2.fastq.gz out1=~{samplename}_1.clean.fastq.gz out2=~{samplename}_2.clean.fastq.gz outm=~{samplename}.matched_phix.fq ref=${phix_fasta} k=31 hdist=1 stats=~{samplename}.phix.stats.txt ordered=t
  >>>
  output {
    File read1_clean = "${samplename}_1.clean.fastq.gz"
    File read2_clean = "${samplename}_2.clean.fastq.gz"
    File adapter_stats = "${samplename}.adapters.stats.txt"
    File phiX_stats = "${samplename}.phix.stats.txt"
    String bbduk_docker = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task bbduk_se {
  input {
    File read1_trimmed
    String samplename
    String docker = "quay.io/staphb/bbtools:38.76"
    Int memory = 8
    Int disk_size = 100
    File? adapters
    File? phix
  }
  command <<<
    # date and version control
    date | tee DATE

    # set adapter fasta
    if [[ ! -z "~{adapters}" ]]; then
      echo "Using user supplied FASTA file for adapters..."
      adapter_fasta="~{adapters}"
    else
      echo "User did not supply adapters FASTA file, using default adapters.fa file..."
      adapter_fasta="/bbmap/resources/adapters.fa" 
    fi

    # set phix fasta
    if [[ ! -z "~{phix}" ]]; then
      echo "Using user supplied FASTA file for phiX..."
      phix_fasta="~{phix}"
    else
      echo "User did not supply phiX FASTA file, using default phix174_ill.ref.fa.gz file..."
      phix_fasta="/bbmap/resources/phix174_ill.ref.fa.gz"
    fi

    bbduk.sh in1=~{read1_trimmed} out1=~{samplename}.rmadpt_1.fastq.gz ref=${adapter_fasta} stats=~{samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo ordered=t

    bbduk.sh in1=~{read1_trimmed} out1=~{samplename}_1.clean.fastq.gz outm=~{samplename}.matched_phix.fq ref=${phix_fasta} k=31 hdist=1 stats=~{samplename}.phix.stats.txt ordered=t
  >>>
  output {
    File read1_clean = "${samplename}_1.clean.fastq.gz"
    File adapter_stats = "${samplename}.adapters.stats.txt"
    File phiX_stats = "${samplename}.phix.stats.txt"
    String bbduk_docker   = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
