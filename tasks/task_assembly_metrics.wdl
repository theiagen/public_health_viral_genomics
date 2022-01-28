version 1.0 

task stats_n_coverage {

  input {
    File        bamfile
    String      samplename
  }

  command{
    date | tee DATE
    samtools --version | head -n1 | tee VERSION

    samtools stats ${bamfile} > ${samplename}.stats.txt

    samtools coverage ${bamfile} -m -o ${samplename}.cov.hist
    samtools coverage ${bamfile} -o ${samplename}.cov.txt
    samtools flagstat ${bamfile} > ${samplename}.flagstat.txt
    
    coverage=$(cut -f 6 ${samplename}.cov.txt | tail -n 1)
    depth=$(cut -f 7 ${samplename}.cov.txt | tail -n 1)
    meanbaseq=$(cut -f 8 ${samplename}.cov.txt | tail -n 1)
    meanmapq=$(cut -f 9 ${samplename}.cov.txt | tail -n 1)
    
    samtools index ${bamfile} && samtools coverage -r "MN908947.3:21563-25384" ${bamfile} >> ${samplename}.cov.txt
    s_gene_depth=$(cut -f 7 ${samplename}.cov.txt | tail -n 1)

    if [ -z "$coverage" ] ; then coverage="0" ; fi
    if [ -z "s_gene_depth" ] ; then s_gene_depth="0"; fi
    if [ -z "$depth" ] ; then depth="0" ; fi
    if [ -z "$meanbaseq" ] ; then meanbaseq="0" ; fi
    if [ -z "$meanmapq" ] ; then meanmapq="0" ; fi

    echo $coverage | tee COVERAGE
    echo $depth | tee DEPTH 
    echo $s_gene_depth | tee S_GENE_DEPTH
    echo $meanbaseq | tee MEANBASEQ 
    echo $meanmapq | tee MEANMAPQ 
  }

  output {
    String     date = read_string("DATE")
    String     samtools_version = read_string("VERSION") 
    File       stats = "${samplename}.stats.txt"
    File       cov_hist = "${samplename}.cov.hist"
    File       cov_stats = "${samplename}.cov.txt"
    File       flagstat = "${samplename}.flagstat.txt"
    Float      coverage = read_string("COVERAGE")
    Float      depth = read_string("DEPTH")
    Float      s_gene_depth = read_string("S_GENE_DEPTH")
    Float      meanbaseq = read_string("MEANBASEQ")
    Float      meanmapq = read_string("MEANMAPQ")
  }

  runtime {
    docker:       "quay.io/staphb/samtools:1.10"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
