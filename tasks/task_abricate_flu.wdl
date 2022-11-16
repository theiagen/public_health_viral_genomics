version 1.0

task abricate_flu {
  input {
    File assembly
    String samplename
    String database="insaflu"
    # Parameters 
    # --minid Minimum DNA %identity [80]
    # --mincov Minimum DNA %coverage [80]
    Int? minid
    Int? mincov
    Int cpu = 2
    Int memory = 4
    String docker = "staphb/abricate:1.0.1-insaflu-220727"
  }
  command <<<
    date | tee DATE    
    abricate -v | tee ABRICATE_VERSION
    # run abricate
    abricate \
      --db ~{database} \
      ~{'--minid ' + minid} \
      ~{'--mincov ' + mincov} \
      --threads ~{cpu} \
      --nopath \
      ~{assembly} > ~{samplename}_abricate_hits.tsv
    # capturing flu type (A or B) and subtype (e.g. H1 and N1)
    grep "M1" ~{samplename}_abricate_hits.tsv | awk -F '\t' '{ print $15 }' > FLU_TYPE
    HA_hit=$(grep "HA" ~{samplename}_abricate_hits.tsv | awk -F '\t' '{ print $15 }')
    NA_hit=$(grep 'NA' ~{samplename}_abricate_hits.tsv | awk -F '\t' '{ print $15 }')
    echo "${HA_hit}${NA_hit}" >  FLU_SUBTYPE
  >>>
  output {
      String abricate_flu_type = read_string("FLU_TYPE")
      String abricate_flu_subtype = read_string("FLU_SUBTYPE")
      File abricate_flu_results = "~{samplename}_abricate_hits.tsv"
      String abricate_flu_database = database
      String abricate_flu_version = read_string("ABRICATE_VERSION")
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk 100 SSD"
      preemptible:  0
  }
}