version 1.0

task abricate_flu {
  input {
    File assembly
    String samplename
    String database = "insaflu"
    String nextclade_flu_h1n1_ha_tag
    String nextclade_flu_h1n1_na_tag
    String nextclade_flu_h3n2_ha_tag
    String nextclade_flu_h3n2_na_tag
    String nextclade_flu_vic_ha_tag
    String nextclade_flu_vic_na_tag
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
    # capturing flu type (A or B) and subtype (e.g. H1 and N1)
    grep "M1" ~{samplename}_abricate_hits.tsv | awk -F '\t' '{ print $15 }' > FLU_TYPE
    HA_hit=$(grep "HA" ~{samplename}_abricate_hits.tsv | awk -F '\t' '{ print $15 }')
    NA_hit=$(grep 'NA' ~{samplename}_abricate_hits.tsv | awk -F '\t' '{ print $15 }')
    flu_subtype="${HA_hit}${NA_hit}" && echo "$flu_subtype" >  FLU_SUBTYPE
    # set nextclade variables based on subptype
    run_nextclade=true
    touch NEXTCLADE_REF_HA NEXTCLADE_REF_NA NEXTCLADE_NAME_HA NEXTCLADE_NAME_NA NEXTCLADE_DS_TAG_HA NEXTCLADE_DS_TAG_NA
    if [ "${flu_subtype}" == "H1N1" ]; then
      echo "flu_h1n1pdm_ha" > NEXTCLADE_NAME_HA
      echo "MW626062" > NEXTCLADE_REF_HA
      echo "~{nextclade_flu_h1n1_ha_tag}" > NEXTCLADE_DS_TAG_HA
      echo "flu_h1n1pdm_na" > NEXTCLADE_NAME_NA
      echo "MW626056" > NEXTCLADE_REF_NA
      echo "~{nextclade_flu_h1n1_na_tag}" > NEXTCLADE_DS_TAG_NA
    elif [ "${flu_subtype}" == "H3N2" ]; then
      echo "flu_h3n2_ha" > NEXTCLADE_NAME_HA
      echo "EPI1857216" > NEXTCLADE_REF_HA
      echo "~{nextclade_flu_h3n2_ha_tag}" > NEXTCLADE_DS_TAG_HA
      echo "flu_h3n2_na" > NEXTCLADE_NAME_NA
      echo "EPI1857215" > NEXTCLADE_REF_NA
      echo "~{nextclade_flu_h3n2_na_tag}" > NEXTCLADE_DS_TAG_NA
    elif [ "${flu_subtype}" == "Victoria" ]; then
      echo "flu_vic_ha" > NEXTCLADE_NAME_HA
      echo "KX058884" > NEXTCLADE_REF_HA
      echo "~{nextclade_flu_vic_ha_tag}" > NEXTCLADE_DS_TAG_HA
      echo "flu_vic_na" > NEXTCLADE_NAME_NA
      echo "CY073894" > NEXTCLADE_REF_NA
      echo "~{nextclade_flu_vic_na_tag}" > NEXTCLADE_DS_TAG_NA
    elif [ "${flu_subtype}" == "Yamagata" ]; then
      echo "flu_yam_ha" > NEXTCLADE_NAME_HA
      echo "JN993010" > NEXTCLADE_REF_HA
      echo "~{nextclade_flu_yam_tag}" > NEXTCLADE_DS_TAG_HA 
      # this makes no biological sense, but avoids errors with nextclade
      echo "flu_vic_na" > NEXTCLADE_NAME_NA
      echo "CY073894" > NEXTCLADE_REF_NA
      echo "~{nextclade_flu_vic_na_tag}" > NEXTCLADE_DS_TAG_NA
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
      String nextclade_ref_ha = read_string("NEXTCLADE_REF_HA")
      String nextclade_name_ha = read_string("NEXTCLADE_NAME_HA")
      String nextclade_ds_tag_ha = read_string("NEXTCLADE_DS_TAG_HA")
      String nextclade_ref_na = read_string("NEXTCLADE_REF_NA")
      String nextclade_name_na = read_string("NEXTCLADE_NAME_NA")
      String nextclade_ds_tag_na = read_string("NEXTCLADE_DS_TAG_NA")

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