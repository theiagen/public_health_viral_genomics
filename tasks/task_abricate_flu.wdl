version 1.0

task abricate_flu {
  input {
    File assembly
    String samplename
    String database="insaflu"
    String nextclade_flu_h1n1_tag="2022-06-08T12:00:00Z"
    String nextclade_flu_h3n2_tag="2022-06-08T12:00:00Z"
    String nextclade_flu_vic_tag="2022-06-08T12:00:00Z"
    String nextclade_flu_yam_tag="2022-07-27T12:00:00Z"
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
    flu_subtype="${HA_hit}${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    # set nextstrain variables based on subptype
    run_nextstrain=true
    touch NEXTSTRAIN_REF NEXTSTRAIN_NAME NEXTSTRAIN_DS_TAG
    if [ "${flu_subtype}" == "H1N1" ]; then
      echo "flu_h1n1pdm_ha" > NEXTSTRAIN_NAME
      echo "CY121680" > NEXTSTRAIN_REF
      echo "~{nextclade_flu_h1n1_tag}" > NEXTSTRAIN_DS_TAG
      > nextstrain_dataset tag
    elif [ "${flu_subtype}" == "H3N2" ]; then
      echo "flu_h3n2_ha" > NEXTSTRAIN_NAME
      echo "CY163680" > NEXTSTRAIN_REF
      echo "~{nextclade_flu_h3n2_tag}" > NEXTSTRAIN_DS_TAG
    elif [ "${flu_subtype}" == "Victoria" ]; then
      echo "flu_vic_ha" > NEXTSTRAIN_NAME
      echo "KX058884" > NEXTSTRAIN_REF
      echo "~{nextclade_flu_vic_tag}" > NEXTSTRAIN_DS_TAG
    elif [ "${flu_subtype}" == "Yamagata" ]; then
      echo "flu_yam_ha" > NEXTSTRAIN_NAME
      echo "JN993010" > NEXTSTRAIN_REF
      echo "~{nextclade_flu_yam_tag}" > NEXTSTRAIN_DS_TAG 
    else 
      run_nextstrain=false 
    fi
    echo ${run_nextstrain} > RUN_NEXTSTRAIN
  >>>
  output {
      String abricate_flu_type = read_string("FLU_TYPE")
      String abricate_flu_subtype = read_string("FLU_SUBTYPE")
      File abricate_flu_results = "~{samplename}_abricate_hits.tsv"
      String abricate_flu_database = database
      String abricate_flu_version = read_string("ABRICATE_VERSION")
      Boolean run_nextstrain = read_boolean("RUN_NEXTSTRAIN")
      String nextstrain_ref = read_string("NEXTSTRAIN_REF")
      String nextstrain_name = read_string("NEXTSTRAIN_NAME")
      String nextstrain_tag = read_string("NEXTSTRAIN_DS_TAG")

  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk 100 SSD"
      preemptible:  0
  }
}