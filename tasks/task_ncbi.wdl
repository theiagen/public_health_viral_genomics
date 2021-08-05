version 1.0

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. See https://github.com/ncbi/vadr/wiki/Coronavirus-annotation"
  }
  input {
    File   genome_fasta
    String vadr_opts="--glsearch -s -r --nomisc --mkey sarscov2 --alt_fail lowscore,fstukcnf,insertnn,deletinn --mdir /opt/vadr/vadr-models/"
    Int assembly_length_unambiguous
    Int skip_length=10000

    String  docker="staphb/vadr:1.2.1"
    Int minlen=50
    Int maxlen=30000
  }
  String out_base = basename(genome_fasta, '.fasta')
  command <<<
    set -e

  if [ ~{assembly_length_unambiguous} -gt ~{skip_length} ]; then 

      # remove terminal ambiguous nucleotides
      /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
        "~{genome_fasta}" \
        --minlen ~{minlen} \
        --maxlen ~{maxlen} \
        > "~{out_base}_trimmed.fasta"

      # run VADR
      v-annotate.pl \
        ~{vadr_opts} \
        "~{out_base}_trimmed.fasta" \
        "~{out_base}"


      # package everything for output
      tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .

      # prep alerts into a tsv file for parsing
      cat "~{out_base}/~{out_base}.vadr.alt.list" | cut -f 2 | tail -n +2 > "~{out_base}.vadr.alerts.tsv"
      cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS
      
    else
      echo "VADR skipped due to poor assembly; assembly length (unambiguous) = ~{assembly_length_unambiguous}" > NUM_ALERTS

    fi

  >>>
  output {
    File? feature_tbl  = "~{out_base}/~{out_base}.vadr.pass.tbl"
    String  num_alerts = read_string("NUM_ALERTS")
    File? alerts_list = "~{out_base}/~{out_base}.vadr.alt.list"
    File? outputs_tgz = "~{out_base}.vadr.tar.gz"
    String vadr_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "2 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}
