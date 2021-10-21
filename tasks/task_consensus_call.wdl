version 1.0

task primer_trim {

  input {
    File     bamfile
    String   samplename
    File     primer_bed
    Boolean? keep_noprimer_reads=true
  }
  String primer_name = basename(primer_bed)

  command {
    # date and version control
    echo "${primer_name}" | tee PRIMER_NAME
    date | tee DATE
    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # trimming primers
    ivar trim \
    ${true="-e" false="" keep_noprimer_reads} \
    -i ${bamfile} \
    -b ${primer_bed} \
    -p ${samplename}.primertrim | tee IVAR_OUT

    # sorting and indexing the trimmed bams
    samtools sort \
    ${samplename}.primertrim.bam \
    -o ${samplename}.primertrim.sorted.bam

    samtools index ${samplename}.primertrim.sorted.bam

    PCT=$(grep "Trimmed primers from" IVAR_OUT | perl -lape 's/Trimmed primers from (\S+)%.*/$1/')
    echo $PCT
    if [[ $PCT = -* ]]; then echo 0; else echo $PCT; fi > IVAR_TRIM_PCT
  }

  output {
    File   trimmed_bam = "${samplename}.primertrim.bam"
    File   trim_sorted_bam = "${samplename}.primertrim.sorted.bam"
    File   trim_sorted_bai = "${samplename}.primertrim.sorted.bam.bai"
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String pipeline_date = read_string("DATE")
    Float  primer_trimmed_read_percent = read_float("IVAR_TRIM_PCT")
    String primer_bed_name = read_string("PRIMER_NAME")
  }

  runtime {
    docker:       "quay.io/staphb/ivar:1.3.1-titan"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task variant_call {

  input {
    File     bamfile
    String   samplename
    String?  ref_genome = "/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"
    String?  ref_gff = "/reference/GCF_009858895.2_ASM985889v3_genomic.gff"
    Boolean? count_orphans = true
    Int?     max_depth = "600000"
    Boolean? disable_baq = true
    Int?     min_bq = "0"
    Int?     min_qual = "20"
    Float?   min_freq = "0.6"
    Int?     min_depth = "10"
  }

  command {
    # date and version control
    date | tee DATE
    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # call variants
    samtools mpileup \
    ${true = "-A" false = "" count_orphans} \
    -d ${max_depth} \
    ${true = "-B" false = "" disable_baq} \
    -Q ${min_bq} \
    --reference ${ref_genome} \
    ${bamfile} | \
    ivar variants \
    -p ${samplename}.variants \
    -q ${min_qual} \
    -t ${min_freq} \
    -m ${min_depth} \
    -r ${ref_genome} \
    -g ${ref_gff}

    # Convert TSV to VCF
    ivar_variants_to_vcf.py ${samplename}.variants.tsv ${samplename}.variants.vcf

    variants_num=$(grep "TRUE" ${samplename}.variants.tsv | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
    echo $variants_num | tee VARIANT_NUM
	}

  output {
    Int       variant_num = read_string("VARIANT_NUM")
    File      sample_variants_tsv = "${samplename}.variants.tsv"
    File      sample_variants_vcf = "${samplename}.variants.vcf"
    String    ivar_version = read_string("IVAR_VERSION")
    String    samtools_version = read_string("SAMTOOLS_VERSION")
    String    pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "quay.io/staphb/ivar:1.3.1-titan"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task consensus {

  input {
    File        bamfile
    String      samplename
    String?     ref_genome = "/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"
    String?     ref_gff = "/reference/GCF_009858895.2_ASM985889v3_genomic.gff"
    Boolean?    count_orphans = true
    Int?        max_depth = "600000"
    Boolean?    disable_baq = true
    Int?        min_bq = "0"
    Int?        min_qual = "20"
    Float?      min_freq = "0.6"
    Int?        min_depth = "10"
    String?     char_unknown = "N"
  }

  command {
    # date and version control
    date | tee DATE
    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # call consensus
    samtools mpileup \
    ${true = "--count-orphans" false = "" count_orphans} \
    -d ${max_depth} \
    ${true = "--no-BAQ" false = "" disable_baq} \
    -Q ${min_bq} \
    --reference ${ref_genome} \
    ${bamfile} | \
    ivar consensus \
    -p ${samplename}.consensus \
    -q ${min_qual} \
    -t ${min_freq} \
    -m ${min_depth} \
    -n ${char_unknown}

    # clean up fasta header
    echo ">${samplename}" > ${samplename}.ivar.consensus.fasta
    grep -v ">" ~{samplename}.consensus.fa >> ${samplename}.ivar.consensus.fasta
  }

  output {
    File      consensus_seq = "${samplename}.ivar.consensus.fasta"
    String    ivar_version = read_string("IVAR_VERSION")
    String    samtools_version = read_string("SAMTOOLS_VERSION")
    String    pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "quay.io/staphb/ivar:1.3.1-titan"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
