version 1.0 

task bedtools_cov {
  
  input {
    File       bamfile
    File       baifile
    String?    primer_bed = "/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019_amplicon.bed"
    String?    fail_threshold = 20
  }
  
  command <<<
    # date and version control
    date | tee DATE
    bedtools --version | tee VERSION
    cp ~{bamfile} ./
    cp ~{baifile} ./

    bedtools coverage -a ~{primer_bed} -b $(ls *bam) > amplicon_coverage.txt
    bedtools coverage -a ~{primer_bed} -b $(ls *bam) | cut -f 6 | awk '{if ( $1 < 20 ) print $0 }' | wc -l | tee AMP_FAIL
  >>>

  output {
    String     date = read_string("DATE")
    String     version = read_string("VERSION") 
    Int        amp_fail = read_string("AMP_FAIL")
    File       amp_coverage = "amplicon_coverage.txt"
  }

  runtime {
    docker:       "staphb/ivar:1.2.2_artic20200528"
    memory:       "2 GB"
    cpu:          1
    disks:        "local-disk 100 SSD"
    preemptible:  0      
  }
}

task bedtools_multicov {
  
  input {
    Array[File]  bamfiles
    Array[File]  baifiles
    Array[File]  primtrim_bamfiles
    Array[File]  primtrim_baifiles
    String?      primer_bed = "/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019_amplicon.bed"
  }
  
  command{
    # date and version control
    date | tee DATE
    bedtools --version | tee VERSION
    cp ${sep=" " bamfiles} ./
    cp ${sep=" " baifiles} ./
    cp ${sep=" " primtrim_bamfiles} ./
    cp ${sep=" " primtrim_baifiles} ./

    echo "primer" $(ls *bam) | tr ' ' '\t' > multicov.txt
    bedtools multicov -bams $(ls *bam) -bed ${primer_bed} | cut -f 4,6- >> multicov.txt
  }

  output {
    String     date = read_string("DATE")
    String     version = read_string("VERSION") 
    File       amp_coverage = "multicov.txt"
  }

  runtime {
    docker:       "staphb/ivar:1.2.2_artic20200528"
    memory:       "2 GB"
    cpu:          1
    disks:        "local-disk 100 SSD"
    preemptible:  0      
  }
}
