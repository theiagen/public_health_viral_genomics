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
  }
}

task read_filtering {

  input {
    File demultiplexed_reads
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

    docker:       "theiagen/artic-ncov2019:1.1.3"
    memory:       "16 GB"
    cpu:          8
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task consensus {
  ## Need to output multiple directories
  input {
    String    samplename
    File      filtered_reads
    String?    artic_primer_version="V3"
    Int?      normalise=20000
    Int?      cpu=8
  }

  command{
    # version control
    echo "Medaka via $(artic -v)" | tee VERSION
    artic minion --medaka --normalise ${normalise} --threads ${cpu} --scheme-directory /artic-ncov2019/primer_schemes --read-file ${filtered_reads} nCoV-2019/${artic_primer_version} ${samplename}

    num_N=$( grep -v ">" ${samplename}.consensus.fasta | grep -o 'N' | wc -l )
    if [ -z "$num_N" ] ; then num_N="0" ; fi
    echo $num_N | tee NUM_N

    num_ACTG=$( grep -v ">" ${samplename}.consensus.fasta | grep -o -E "C|A|T|G" | wc -l )
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    echo $num_ACTG | tee NUM_ACTG

    num_degenerate=$( grep -v ">" ${samplename}.consensus.fasta | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
    echo $num_degenerate | tee NUM_DEGENERATE

    num_total=$( grep -v ">" ${samplename}.consensus.fasta | grep -o -E '[A-Z]' | wc -l )
    if [ -z "$num_total" ] ; then num_total="0" ; fi
    echo $num_total | tee NUM_TOTAL

    count_pool_1=$(samtools view ~{samplename}.primertrimmed.nCoV-2019_1.sorted.bam | wc -l)
    count_pool_2=$(samtools view ~{samplename}.primertrimmed.nCoV-2019_2.sorted.bam | wc -l)
    count_total=$(( $count_pool_1 + $count_pool_2 ))

    python -c "print ( round(($count_pool_1 / $count_total ) * 100, 2) )" | tee POOL1_PERCENT
    python -c "print ( round(($count_pool_2 / $count_total ) * 100, 2) )" | tee POOL2_PERCENT

    # clean up fasta header
    echo ">${samplename}" > ${samplename}.medaka.consensus.fasta
    grep -v ">" ${samplename}.consensus.fasta >> ${samplename}.medaka.consensus.fasta

  }

  output {
    File    consensus_seq = "${samplename}.medaka.consensus.fasta"
    File    sorted_bam = "${samplename}.trimmed.rg.sorted.bam"
    File    trim_sorted_bam = "${samplename}.primertrimmed.rg.sorted.bam"
    File    trim_sorted_bai = "${samplename}.primertrimmed.rg.sorted.bam.bai"
    Int     number_N = read_string("NUM_N")
    Int     number_ATCG = read_string("NUM_ACTG")
    Int     number_Degenerate = read_string("NUM_DEGENERATE")
    Int     number_Total = read_string("NUM_TOTAL")
    Float     pool1_percent = read_string("POOL1_PERCENT")
    Float     pool2_percent = read_string("POOL2_PERCENT")
    String  artic_pipeline_version = read_string("VERSION")
  }

  runtime {
    docker:       "theiagen/artic-ncov2019:1.1.3"
    memory:       "16 GB"
    cpu:          8
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
