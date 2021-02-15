version 1.0

task sample_metrics_v1 {

  input {
    String    samplename
    String    submission_id
    String    collection_date
    String    pangolin_lineage
    String    pangolin_aLRT
    String?    fastqc_raw_pairs
    String?    seqy_pairs
    String?    seqy_percent
    String?    kraken_human
    String?    kraken_sc2
    String    variant_num
    String    number_N
    String    number_ATCG
    String    number_Degenerate
    String    number_Total
    String    coverage
    String    depth
    String    meanbaseq_trim
    String    meanmapq_trim
    String    coverage_trim
    String    depth_trim
    String    amp_fail
    Float?    coverage_threshold = 95.00
    Float?    meanbaseq_threshold = 30.00
    Float?    meanmapq_threshold = 30.00
  }

  command {
    echo "sample_id,deidentified_id,collection_date,\
    pangolin_lineage,pangolin_probability,\
    fastqc_raw_pairs,pairs_after_cleaning,percent_kept_after_cleaning,\
    depth_before_trimming,depth_after_trimming,coverage_before_trimming,coverage_after_trimming,\
    %_human_reads,%_SARS-COV-2_reads,num_failed_amplicons,num_variants,\
    num_N,num_degenerate,num_ACTG,num_total,meanbaseq_trim,meanmapq_trim,assembly_status"

    # Determining assembly status based upon coverage & mapping quality provided
    cov=$(echo "${coverage} >= ${coverage_threshold}" | bc)
    if (($cov == 0)); then cov_status="Avg coverage < ${coverage_threshold}%;"; else cov_status="";fi

    baseq=$(echo "${meanbaseq_trim} >= ${meanbaseq_threshold}" | bc)
    if (($baseq == 0)); then baseq_status="Mean base quality < ${meanbaseq_threshold};"; else baseq_status=""; fi

    mapq=$(echo "${meanmapq_trim} >= ${meanmapq_threshold}" | bc)
    if (($mapq == 0)); then mapq_status="Mean map quality < ${meanmapq_threshold}"; else mapq_status=""; fi

    if (($cov == 1)) && (($baseq == 1)) && (($mapq == 1)); then
      assembly_status="PASS"
    else
      assembly_status=`echo "WARNING: $(cov_status)$(baseq_status)$(mapq_status)"`
    fi

    echo "${samplename},${submission_id},${collection_date},\
    ${pangolin_lineage},${pangolin_aLRT},\
    ${fastqc_raw_pairs},${seqy_pairs},${seqy_percent},\
    ${depth},${depth_trim},${coverage},${coverage_trim},\
    ${kraken_human},${kraken_sc2},${amp_fail},${variant_num},\
    ${number_N},${number_Degenerate},${number_ATCG},${number_Total},\
    ${meanbaseq_trim},${meanmapq_trim},$(assembly_status)" | tee SAMPLE_METRICS
  }

  output {
    String  single_metrics = read_string("SAMPLE_METRICS")
  }

  runtime {
      docker:       "staphb/multiqc:1.7"
      memory:       "1 GB"
      cpu:          1
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}

task sample_metrics {

  input {
    String    samplename
    String    submission_id
    String    collection_date
    String    pangolin_lineage
    Float     pangolin_aLRT
    String    nextclade_clade
    String    nextclade_aa_subs
    String    nextclade_aa_dels
    Int?       fastqc_raw_pairs
    Int?       seqy_pairs
    Float?     seqy_percent
    Float?     kraken_human
    Float?     kraken_sc2
    Int       variant_num
    Int       number_N
    Int       number_ATCG
    Int       number_Degenerate
    Int       number_Total
    Float     coverage
    Float     depth
    Float     meanbaseq_trim
    Float     meanmapq_trim
    Float     coverage_trim
    Float     depth_trim
    Int       amp_fail
    Float?    coverage_threshold = 95.00
    Float?    meanbaseq_threshold = 30.00
    Float?    meanmapq_threshold = 30.00
  }

  command <<<
  python3<<CODE

  if ~{coverage_trim} >= ~{coverage_threshold} and ~{meanbaseq_trim} >= ~{meanbaseq_threshold} and ~{meanmapq_trim} >= ~{meanmapq_threshold}:
    assembly_status = "PASS"
  else:
    assembly_status = "Warning: "

  if ~{coverage_trim} <= ~{coverage_threshold}:
      assembly_status += "Avg coverage < 95%; "
  if ~{meanbaseq_trim} <= ~{meanbaseq_threshold}:
      assembly_status += "Mean base quality < 30; "
  if ~{meanmapq_trim} <= ~{meanmapq_threshold}:
      assembly_status += "Mean map quality < 30"

  outstring="~{samplename},~{submission_id},~{collection_date},\
  ~{pangolin_lineage},~{pangolin_aLRT},\
  ~{nextclade_clade},~{nextclade_aa_subs},~{nextclade_aa_dels},\
  ~{fastqc_raw_pairs},~{seqy_pairs},~{seqy_percent},\
  ~{depth},~{depth_trim},~{coverage},~{coverage_trim},\
  ~{kraken_human},~{kraken_sc2},~{amp_fail},~{variant_num},\
  ~{number_N},~{number_Degenerate},~{number_ATCG},~{number_Total},\
  ~{meanbaseq_trim},~{meanmapq_trim}," + assembly_status

  print(outstring)

  CODE
  >>>

  output {
    String  single_metrics = read_string(stdout())
  }

  runtime {
      docker:       "staphb/multiqc:1.7"
      memory:       "1 GB"
      cpu:          1
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}


task merge_metrics {

  input {
    Array[String]   all_metrics
  }

  command {
    echo "sample_id,deidentified_id,collection_date,\
    pangolin_lineage,pangolin_probability,\
    nextclade_lineage,nextclade_aaSubstitutions,nextclade_aaDeletions,\
    raw_pairs,pairs_after_cleaning,percent_kept_after_cleaning,\
    depth_before_trimming,depth_after_trimming,coverage_before_trimming,coverage_after_trimming,\
    %_human_reads,%_SARS-COV-2_reads,num_failed_amplicons,num_variants,\
    num_N,num_degenerate,num_ACTG,num_total,\
    meanbaseq_trim,meanmapq_trim,assembly_status" >> run_results.csv

    echo "${sep="END" all_metrics}" >> run_results.csv
    sed -i "s/END/\n/g" run_results.csv
  }

  output {
    File    run_results = "run_results.csv"
  }

  runtime {
      docker:       "staphb/seqyclean:1.10.09"
      memory:       "1 GB"
      cpu:          1
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}
