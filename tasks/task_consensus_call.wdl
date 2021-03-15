version 1.0

task primer_trim {

  input {
    File        bamfile
    String      samplename
    String?     primer_bed = "/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.bed"
    Boolean?    keep_noprimer_reads=true
  }

  command {
    # date and version control
    date | tee DATE
    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # trimming primers
    ivar trim \
    ${true="-e" false="" keep_noprimer_reads} \
    -i ${bamfile} \
    -b ${primer_bed} \
    -p ${samplename}.primertrim

    # sorting and indexing the trimmed bams
    samtools sort \
    ${samplename}.primertrim.bam \
    -o ${samplename}.primertrim.sorted.bam

    samtools index ${samplename}.primertrim.sorted.bam
  }

  output {
    File      trimmed_bam = "${samplename}.primertrim.bam"
    File 	  trim_sorted_bam = "${samplename}.primertrim.sorted.bam"
    File 	  trim_sorted_bai = "${samplename}.primertrim.sorted.bam.bai"
    String    ivar_version = read_string("IVAR_VERSION")
    String 	  samtools_version = read_string("SAMTOOLS_VERSION")
    String    pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "staphb/ivar:1.2.2_artic20200528"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task variant_call {

  input {
    File        bamfile
    String      samplename
    String? 	ref_genome = "/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"
    String? 	ref_gff = "/reference/GCF_009858895.2_ASM985889v3_genomic.gff"
    Boolean?    count_orphans = true
    Int?        max_depth = "600000"
    Boolean?    disable_baq = true
    Int?        min_bq = "0"
    Int?        min_qual = "20"
    Float?      min_freq = "0.6"
    Int?        min_depth = "10"
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

    variants_num=$(grep "TRUE" ${samplename}.variants.tsv | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
    echo $variants_num | tee VARIANT_NUM
	}

  output {
 	Int       variant_num = read_string("VARIANT_NUM")
 	File  	  sample_variants = "${samplename}.variants.tsv"
    String    ivar_version = read_string("IVAR_VERSION")
    String    samtools_version = read_string("SAMTOOLS_VERSION")
    String    pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "staphb/ivar:1.2.2_artic20200528"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
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

    num_N=$( grep -v ">" ${samplename}.consensus.fa | grep -o 'N' | wc -l )
    if [ -z "$num_N" ] ; then num_N="0" ; fi
    echo $num_N | tee NUM_N

    num_ACTG=$( grep -v ">" ${samplename}.consensus.fa | grep -o -E "C|A|T|G" | wc -l )
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    echo $num_ACTG | tee NUM_ACTG

    num_degenerate=$( grep -v ">" ${samplename}.consensus.fa | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
    echo $num_degenerate | tee NUM_DEGENERATE

    num_total=$( grep -v ">" ${samplename}.consensus.fa | grep -o -E '[A-Z]' | wc -l )
    if [ -z "$num_total" ] ; then num_total="0" ; fi
    echo $num_total | tee NUM_TOTAL

    # clean up fasta header
    echo ">${samplename}" > ${samplename}.ivar.consensus.fasta
    grep -v ">" ~{samplename}.consensus.fa >> ${samplename}.ivar.consensus.fasta
  }

  output {
    File      consensus_seq = "${samplename}.ivar.consensus.fasta"
    Int       number_N = read_string("NUM_N")
    Int       number_ATCG = read_string("NUM_ACTG")
    Int       number_Degenerate = read_string("NUM_DEGENERATE")
    Int       number_Total = read_string("NUM_TOTAL")
    String    ivar_version = read_string("IVAR_VERSION")
    String    samtools_version = read_string("SAMTOOLS_VERSION")
    String    pipeline_date = read_string("DATE")
  }

  runtime {
    docker:       "staphb/ivar:1.2.2_artic20200528"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
