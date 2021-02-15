versoin 1.0 

task local_TSV_convert {

  input {
    String  samplename
    String  read1
    String  read2
  }

  command {
  	echo sample_info[0] | tee SAMPLENAME
    echo $PWD >> out.txt
    ls ./ >> out.txt
#  	ls ${read1} >> out.txt
#    ls ${read2} >> out.txt
  }

  output {
    String     sample_name = read_string("SAMPLENAME")
    File       files = "out.txt"
    # File       sample_R1 = read_string("READ1")
    # File       sample_R2 = read_string("READ2")
  }

  runtime {
    docker:       "staphb/ivar:1.2.2_artic20200528"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0      
  }
}

