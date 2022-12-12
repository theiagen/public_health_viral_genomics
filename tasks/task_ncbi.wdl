version 1.0

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. See https://github.com/ncbi/vadr/wiki/Coronavirus-annotation"
  }
  input {
    File genome_fasta
    String vadr_opts = "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"
    Int assembly_length_unambiguous
    Int skip_length = 10000
    String docker = "staphb/vadr:1.5"
    Int minlen = 50
    Int maxlen = 30000
    Int cpu = 2
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
      # --split and --cpu must be used in conjuction
      v-annotate.pl \
        --split --cpu ~{cpu} \
        ~{vadr_opts} \
        "~{out_base}_trimmed.fasta" \
        "~{out_base}"

      # package everything for output
      tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .

      # package up FASTA files into zip file for output. Note: this will work whether the --out_allfasta flag is included or not (there are just more when the option is used)
      mkdir -v vadr_fasta_files
      cp -v ~{out_base}/*.fa vadr_fasta_files
      zip ~{out_base}_vadr-fasta-files.zip vadr_fasta_files/*.fa 

      # prep alerts into a tsv file for parsing
      cut -f 5 "~{out_base}/~{out_base}.vadr.alt.list" | tail -n +2 > "~{out_base}.vadr.alerts.tsv"
      cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS

    else
      echo "VADR skipped due to poor assembly; assembly length (unambiguous) = ~{assembly_length_unambiguous}" > NUM_ALERTS
    fi

  >>>
  output {
    File? feature_tbl = "~{out_base}/~{out_base}.vadr.pass.tbl"
    String num_alerts = read_string("NUM_ALERTS")
    File? alerts_list = "~{out_base}/~{out_base}.vadr.alt.list"
    File? outputs_tgz = "~{out_base}.vadr.tar.gz"
    File? vadr_fastas_zip_archive = "~{out_base}_vadr-fasta-files.zip"
    String vadr_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}
