version 1.0

task vadr {

  input {
    File genome
    String  samplename

  }
  command {
  predicted_genus=""
  genus_threshold=""
  predicted_species=""
  species_threshold=""

  midas query ~{genome} > midas_out.csv
,top_genus,,top_species,species_threshol
  python3 <<CODE
  with open ("midas_out.csv", 'r') as csv_file:
    csv_reader = list(csv.DictReader(csv_file, delimiter=","))
    for line in csv_reader:
      predicted_genus=line["predicted_genus"]""
      genus_threshold=line["genus_threshold"]
      predicted_species=line["predicted_species"]
      species_threshold=""



  CODE
  cp report.pdf ${cluster_name}_cluster_analysis.pdf
  cp SNP_heatmap.png ${cluster_name}_SNP_heatmap.png
  cp pairwise_snp_list.csv ${cluster_name}_pairwise_snp_list.csv
}
  output {
    Int  num_alerts = read_int("NUM_ALERTS")
    File? vadr_failed = "~{text}.fail"
    File? vadr_passed = "~{text}.pass"
  }
  runtime {
    docker: "staphb/vadr:1.1.2"
    memory: "64 GB"
    cpu: 8
  }
}
