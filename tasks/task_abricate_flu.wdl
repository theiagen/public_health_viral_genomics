version 1.0

task abricate_flu {
  input {
    File assembly
    String samplename
    String database = "insaflu"
    String nextclade_flu_h1n1_tag
    String nextclade_flu_h3n2_tag
    String nextclade_flu_vic_tag
    String nextclade_flu_yam_tag
    Int minid = 70
    Int mincov =60
    Int cpu = 2
    Int memory = 4
    String docker = "staphb/abricate:1.0.1-insaflu-220727"
    Int disk_size = 100
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
    # capturing flu type (A or B based on M1 hit) and subtype (e.g. H1 and N1 based on HA/NA hits)
    ## awk for gene column ($6) to grab subtype ($15)
    cat ~{samplename}_abricate_hits.tsv | awk -F '\t' '{if ($6=="M1") print $15}' > FLU_TYPE
    HA_hit=$(cat ~{samplename}_abricate_hits.tsv | awk -F '\t' '{if ($6=="HA") print $15 }')
    NA_hit=$(cat ~{samplename}_abricate_hits.tsv | awk -F '\t' '{if ($6=="NA") print $15 }')
    flu_subtype="${HA_hit}${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    # set nextclade variables based on subptype
    run_nextclade=true
    touch NEXTCLADE_REF NEXTCLADE_NAME NEXTCLADE_DS_TAG
    if [ "${flu_subtype}" == "H1N1" ]; then
      echo "flu_h1n1pdm_ha" > NEXTCLADE_NAME
      echo "CY121680" > NEXTCLADE_REF
      echo "~{nextclade_flu_h1n1_tag}" > NEXTCLADE_DS_TAG
    elif [ "${flu_subtype}" == "H3N2" ]; then
      echo "flu_h3n2_ha" > NEXTCLADE_NAME
      echo "CY163680" > NEXTCLADE_REF
      echo "~{nextclade_flu_h3n2_tag}" > NEXTCLADE_DS_TAG
    elif [ "${flu_subtype}" == "Victoria" ]; then
      echo "flu_vic_ha" > NEXTCLADE_NAME
      echo "KX058884" > NEXTCLADE_REF
      echo "~{nextclade_flu_vic_tag}" > NEXTCLADE_DS_TAG
    elif [ "${flu_subtype}" == "Yamagata" ]; then
      echo "flu_yam_ha" > NEXTCLADE_NAME
      echo "JN993010" > NEXTCLADE_REF
      echo "~{nextclade_flu_yam_tag}" > NEXTCLADE_DS_TAG 
    else 
      run_nextclade=false 
    fi
    echo ${run_nextclade} > RUN_NEXTCLADE
  >>>
  output {
      String abricate_flu_type = read_string("FLU_TYPE")
      String abricate_flu_subtype = read_string("FLU_SUBTYPE")
      File abricate_flu_results = "~{samplename}_abricate_hits.tsv"
      String abricate_flu_database = database
      String abricate_flu_version = read_string("ABRICATE_VERSION")
      Boolean run_nextclade = read_boolean("RUN_NEXTCLADE")
      String nextclade_ref = read_string("NEXTCLADE_REF")
      String nextclade_name = read_string("NEXTCLADE_NAME")
      String nextclade_ds_tag = read_string("NEXTCLADE_DS_TAG")

  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB"
      preemptible:  0
  }
}