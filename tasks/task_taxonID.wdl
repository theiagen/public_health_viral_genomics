version 1.0

task kraken2 {
  input {
    File read1
    File? read2
    String samplename
    String kraken2_db = "/kraken2-db"
    Int cpu = 4
    String? target_org
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    kraken2 --version | head -n1 | tee VERSION
    num_reads=$(ls *fastq.gz 2> /dev/nul | wc -l)
    if ! [ -z ~{read2} ]; then
      mode="--paired"
    fi
    echo $mode
    kraken2 $mode \
      --threads ~{cpu} \
      --db ~{kraken2_db} \
      ~{read1} ~{read2} \
      --report ~{samplename}_kraken2_report.txt >/dev/null

    percentage_human=$(grep "Homo sapiens" ~{samplename}_kraken2_report.txt | cut -f 1)
     # | tee PERCENT_HUMAN
    percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ~{samplename}_kraken2_report.txt | cut -f1 )
     # | tee PERCENT_COV
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN
    echo $percentage_sc2 | tee PERCENT_SC2
    # capture target org percentage 
    if [ ! -z "~{target_org}" ]; then 
      echo "Target org designated: ~{target_org}"
      percent_target_org=$(grep "~{target_org}" ~{samplename}_kraken2_report.txt | cut -f1 | head -n1 )
      if [-z "$percent_target_org" ] ; then percent_target_org="0" ; fi
    else 
      percent_target_org=""
    fi
    echo $percent_target_org | tee PERCENT_TARGET_ORG

  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File kraken_report = "~{samplename}_kraken2_report.txt"
    Float percent_human = read_float("PERCENT_HUMAN")
    Float percent_sc2 = read_float("PERCENT_SC2")
    String? percent_target_org = read_string("PERCENT_TARGET_ORG")
    String? kraken_target_org = target_org
  }
  runtime {
    docker: "quay.io/staphb/kraken2:2.0.8-beta_hv"
    memory: "8 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task pangolin3 {
  input {
    File fasta
    String samplename
    Int min_length = 10000
    Float max_ambig = 0.5
    String docker = "quay.io/staphb/pangolin:3.1.20-pangolearn-2022-02-02"
    String inference_engine = "usher"
    Int disk_size = 100
  }
  command <<<
    set -e
    # set inference inference_engine
    if [[ "~{inference_engine}" == "usher" ]]
    then
      pango_inference="--usher"
    elif [[ "~{inference_engine}" == "pangolearn" ]]
    then
      pango_inference=""
    else
      echo "unknown inference_engine designated: ~{inference_engine}; must be usher or pangolearn" >&2
      exit 1
    fi

    # date and version capture
    date | tee DATE

    { pangolin --all-versions && usher --version; } | tr '\n' ';'  | cut -f -6 -d ';' | tee VERSION_PANGOLIN_ALL

    pangolin "~{fasta}" $pango_inference \
       --outfile "~{samplename}.pangolin_report.csv" \
       --min-length ~{min_length} \
       --max-ambig ~{max_ambig} \
       --verbose

    python3 <<CODE
    import csv
    #grab output values by column header
    with open("~{samplename}.pangolin_report.csv",'r') as csv_file:
      csv_reader = list(csv.DictReader(csv_file, delimiter=","))
      for line in csv_reader:
        with open("PANGO_ASSIGNMENT_VERSION", 'wt') as lineage:
          pangolin_version=line["pangolin_version"]
          version=line["version"]
          lineage.write(f"pangolin {pangolin_version}; {version}")
        with open("PANGOLIN_LINEAGE", 'wt') as lineage:
          lineage.write(line["lineage"])
        with open("PANGOLIN_CONFLICTS", 'wt') as lineage:
          lineage.write(line["conflict"])
        with open("PANGOLIN_NOTES", 'wt') as lineage:
          lineage.write(line["note"])
    CODE
  >>>
  output {
    String date = read_string("DATE")
    String pangolin_lineage = read_string("PANGOLIN_LINEAGE")
    String pangolin_conflicts = read_string("PANGOLIN_CONFLICTS")
    String pangolin_notes = read_string("PANGOLIN_NOTES")
    String pangolin_assignment_version = read_string("PANGO_ASSIGNMENT_VERSION")
    String pangolin_versions = read_string("VERSION_PANGOLIN_ALL")
    String pangolin_docker = docker
    File pango_lineage_report = "~{samplename}.pangolin_report.csv"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task pangolin4 {
  input {
    File fasta
    String samplename
    Int min_length = 10000
    Float max_ambig = 0.5
    String docker = "staphb/pangolin:4.2-pdata-1.18.1.1"
    String? analysis_mode
    Boolean expanded_lineage=true
    Boolean skip_scorpio=false
    Boolean skip_designation_cache=false
    String? pangolin_arguments
    Int disk_size = 100
  }
  command <<<
    set -e

    # date and version capture
    date | tee DATE

    { pangolin --all-versions && usher --version; } | tr '\n' ';'  | cut -f -6 -d ';' | tee VERSION_PANGOLIN_ALL

    pangolin "~{fasta}" \
       ~{'--analysis-mode ' + analysis_mode} \
       ~{'--min-length ' + min_length} \
       ~{'--max-ambig ' + max_ambig} \
       ~{true='--expanded-lineage' false='' expanded_lineage} \
       ~{true='--skip-scorpio' false='' skip_scorpio} \
       ~{true='--skip-designation-cache' false='' skip_designation_cache} \
       --outfile "~{samplename}.pangolin_report.csv" \
       --verbose \
       ~{pangolin_arguments}


    python3 <<CODE
    import csv
    #grab output values by column header
    with open("~{samplename}.pangolin_report.csv",'r') as csv_file:
      csv_reader = list(csv.DictReader(csv_file, delimiter=","))
      for line in csv_reader:
        with open("PANGO_ASSIGNMENT_VERSION", 'wt') as assignment_version:
          pangolin_version=line["pangolin_version"]
          version=line["version"]
          assignment_version.write(f"pangolin {pangolin_version}; {version}")
        with open("PANGOLIN_LINEAGE", 'wt') as lineage:
          lineage.write(line["lineage"])
        with open("PANGOLIN_CONFLICTS", 'wt') as pangolin_conflicts:
          pangolin_conflicts.write(line["conflict"])
        with open("PANGOLIN_NOTES", 'wt') as pangolin_notes:
          pangolin_notes.write(line["note"])
        with open("EXPANDED_LINEAGE", "wt") as expanded_lineage:
          if "expanded_lineage" in line:
            expanded_lineage.write(line["expanded_lineage"])
          else:
            expanded_lineage.write("NA")
    CODE
  >>>
  output {
    String date = read_string("DATE")
    String pangolin_lineage = read_string("PANGOLIN_LINEAGE")
    String pangolin_lineage_expanded = read_string("EXPANDED_LINEAGE")
    String pangolin_conflicts = read_string("PANGOLIN_CONFLICTS")
    String pangolin_notes = read_string("PANGOLIN_NOTES")
    String pangolin_assignment_version = read_string("PANGO_ASSIGNMENT_VERSION")
    String pangolin_versions = read_string("VERSION_PANGOLIN_ALL")
    String pangolin_docker = docker
    File pango_lineage_report = "~{samplename}.pangolin_report.csv"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task pangolin_update_log {
  input {
    String samplename
    String current_lineage
    String current_pangolin_docker
    String current_pangolin_assignment_version
    String current_pangolin_versions
    String updated_lineage
    String updated_pangolin_docker
    String updated_pangolin_assignment_version
    String updated_pangolin_versions
    String? timezone
    File? lineage_log
    Int disk_size = 100
  }
  command <<<
    # set timezone for date outputs
    ~{default='' 'export TZ=' + timezone}
    DATE=$(date +"%Y-%m-%d")

    #check if lineage has been modified
    if [[ "~{current_lineage}" == "~{updated_lineage}" ]]
    then
      UPDATE_STATUS="pango lineage unchanged: ~{updated_lineage}"
    else
      UPDATE_STATUS="pango lineage modified: ~{current_lineage} -> ~{updated_lineage}"
    fi

    #if a lineage log not provided, create one with headers
    lineage_log_file="~{samplename}_pango_lineage_log.tsv"

    if [ -s "~{lineage_log}" ]
    then
      echo "Lineage log provided"

      if grep -q "previous_pangolin_assignment_version" ~{lineage_log}
      then
        mv "~{lineage_log}" ${lineage_log_file}
      else
        echo "pangolin log file provided not compatible with current PHVG version"
        exit 1
      fi
   else
     echo "Creating new lineage log file as none was provided"
     echo -e "analysis_date\tmodification_status\tprevious_lineage\tprevious_pangolin_docker\tprevious_pangolin_assignment_version\tprevious_pangolin_versions\tupdated_lineage\tupdated_pangolin_docker\tupdated_pangolin_assignment_version\tupdated_pangolin_versions" > ${lineage_log_file}
   fi

     #populate lineage log file
     echo -e "${DATE}\t${UPDATE_STATUS}\t~{current_lineage}\t~{current_pangolin_docker}\t~{current_pangolin_assignment_version}\t~{current_pangolin_versions}\t~{updated_lineage}\t~{updated_pangolin_docker}\t~{updated_pangolin_assignment_version}\t~{updated_pangolin_versions}" >> "${lineage_log_file}"

    echo "${UPDATE_STATUS} (${DATE})"  | tee PANGOLIN_UPDATE
  >>>
  output {
    String pangolin_updates = read_string("PANGOLIN_UPDATE")
    File pango_lineage_log = "~{samplename}_pango_lineage_log.tsv"
  }
  runtime {
    docker: "quay.io/theiagen/utility:1.1"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task nextclade_one_sample {
    meta {
        description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
      File genome_fasta
      File? root_sequence
      File? auspice_reference_tree_json
      File? qc_config_json
      File? gene_annotations_json
      File? pcr_primers_csv
      File? virus_properties
      String docker = "nextstrain/nextclade:2.11.0"
      String dataset_name
      String dataset_reference
      String dataset_tag
      Int disk_size = 50
    }
    String basename = basename(genome_fasta, ".fasta")
    command <<<
        NEXTCLADE_VERSION="$(nextclade --version)"
        echo $NEXTCLADE_VERSION > NEXTCLADE_VERSION

        nextclade dataset get --name="~{dataset_name}" --reference="~{dataset_reference}" --tag="~{dataset_tag}" -o nextclade_dataset_dir --verbose
        set -e
        nextclade run \
            --input-dataset=nextclade_dataset_dir/ \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            ~{"--input-virus-properties " + virus_properties}  \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json \
            --output-all=. \
            "~{genome_fasta}"
    >>>
    runtime {
      docker: "~{docker}"
      memory: "4 GB"
      cpu: 2
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries: 3 
    }
    output {
      String nextclade_version = read_string("NEXTCLADE_VERSION")
      File nextclade_json = "~{basename}.nextclade.json"
      File auspice_json = "~{basename}.nextclade.auspice.json"
      File nextclade_tsv = "~{basename}.nextclade.tsv"
      String nextclade_docker = docker
    }
}

task nextclade_output_parser_one_sample {
    meta {
      description: "Python and bash codeblocks for parsing the output files from Nextclade."
    }
    input {
      File nextclade_tsv
      String docker = "python:slim"
      Int disk_size = 50
      String? organism
      Boolean? NA
      String tamiflu_aa_substitutions = "NA:H275Y,NA:R292K"
    }
    command <<<
      # Set WDL input variable to input.tsv file
      cat "~{nextclade_tsv}" > input.tsv
      touch TAMIFLU_AASUBS

      # Parse outputs using python3
      python3 <<CODE
      import csv
      import codecs

      # list of aa substitutions linked with tamiflu resistance - with the possibility to
      # extend the list with other provived aa substitutions in the format "NA:V95A,NA:I97V"
      tamiflu_aa_subs = ["NA:V95A","NA:I97V","NA:E99A","NA:H101L","NA:G108E",
      "NA:Q116L","NA:V116A","NA:E119D","NA:E119G","NA:E119I","NA:E119V","NA:R136K",
      "NA:T146K","NA:T146P","NA:D151E","NA:N169S","NA:D179N","NA:D197N","NA:D198E",
      "NA:D198G","NA:D198N","NA:A200T","NA:I203M","NA:I203R","NA:I203V","NA:I221T",
      "NA:I222R","NA:I222V","NA:I223R","NA:I223V","NA:S227N","NA:S247N","NA:H255Y",
      "NA:E258Q","NA:H274N","NA:H274Y","NA:H275Y","NA:N275S","NA:H277Y","NA:R292K",
      "NA:N294S","NA:S334N","NA:R371K","NA:D432G","NA:H439P","NA:H439R"]
      
      tamiflu_aa_subs = tamiflu_aa_subs + "~{tamiflu_aa_substitutions}".split(',')

      print(tamiflu_aa_subs)

      def intersection(lst1, lst2):
        # returns intersection between nextclade identified aa substitutions and
        # tamiflu associated aa substitutions
        return list(set(lst1) & set(lst2))

      with codecs.open("./input.tsv",'r') as tsv_file:
        tsv_reader=csv.reader(tsv_file, delimiter="\t")
        tsv_data=list(tsv_reader)

        if len(tsv_data)==1:
          tsv_data.append(['NA']*len(tsv_data[0]))
        tsv_dict=dict(zip(tsv_data[0], tsv_data[1]))
        # parse 'clade_legacy' column if sars-cov-2, if false then parse 'clade' column
        if ("~{organism}" == "sars-cov-2"):
          with codecs.open ("NEXTCLADE_CLADE", 'wt') as Nextclade_Clade:
            nc_clade=tsv_dict['clade_legacy']
            if nc_clade=='':
              nc_clade='NA'
            else:
              nc_clade=nc_clade
            Nextclade_Clade.write(nc_clade)
        else:
          with codecs.open ("NEXTCLADE_CLADE", 'wt') as Nextclade_Clade:
            nc_clade=tsv_dict['clade']
            if nc_clade=='':
              nc_clade='NA'
            else:
              nc_clade=nc_clade
            Nextclade_Clade.write(nc_clade)
        with codecs.open ("NEXTCLADE_AASUBS", 'wt') as Nextclade_AA_Subs:
          nc_aa_subs=tsv_dict['aaSubstitutions']
          if nc_aa_subs=='':
            nc_aa_subs='NA'
          else:
            nc_aa_subs=nc_aa_subs
            # if organism is flu, return list of aa subs associated with tamiflu resistance
            if ("~{organism}" == "flu" and "~{NA}" == "true"):
              tamiflu_subs = intersection(tamiflu_aa_subs,nc_aa_subs.split(','))
              with codecs.open ("TAMIFLU_AASUBS", 'wt') as Tamiflu_AA_Subs:
                Tamiflu_AA_Subs.write(",".join(tamiflu_subs))
          Nextclade_AA_Subs.write(nc_aa_subs)
        with codecs.open ("NEXTCLADE_AADELS", 'wt') as Nextclade_AA_Dels:
          nc_aa_dels=tsv_dict['aaDeletions']
          if nc_aa_dels=='':
            nc_aa_dels='NA'
          else:
            nc_aa_dels=nc_aa_dels
          Nextclade_AA_Dels.write(nc_aa_dels)
        with codecs.open ("NEXTCLADE_LINEAGE", 'wt') as Nextclade_Lineage:
          if 'lineage' in tsv_dict:
            nc_lineage=tsv_dict['lineage']
            if nc_lineage is None:
              nc_lineage=""
            else:
              nc_lineage=nc_lineage
          else:
            nc_lineage=""
          Nextclade_Lineage.write(nc_lineage)
      CODE
    >>>
    runtime {
      docker: "~{docker}"
      memory: "4 GB"
      cpu: 2
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries: 3
    }
    output {
      String nextclade_clade = read_string("NEXTCLADE_CLADE")
      String nextclade_aa_subs = read_string("NEXTCLADE_AASUBS")
      String nextclade_tamiflu_aa_subs = read_string("TAMIFLU_AASUBS")
      String nextclade_aa_dels = read_string("NEXTCLADE_AADELS")
      String nextclade_lineage = read_string("NEXTCLADE_LINEAGE")
    }
}

task freyja_one_sample {
  input {
    File primer_trimmed_bam
    String samplename
    File reference_genome
    File? freyja_usher_barcodes
    File? freyja_lineage_metadata
    Float? eps
    Boolean update_db = false
    Boolean confirmed_only = false
    Boolean bootstrap = false
    Int? number_bootstraps
    Int memory = 4
    String docker = "staphb/freyja:1.3.11"
    Int disk_size = 100
  }
  command <<<
  # capture version
  freyja --version | tee FREYJA_VERSION
  # update freyja reference files if specified
  if ~{update_db}; then 
      freyja update 2>&1 | tee freyja_update.log
      # check log files to ensure update did not fail
      if grep "FileNotFoundError.*lineagePaths.*" freyja_update.log
      then 
        echo "Error in attempting to update Freyja files. Try increasing memory"
        >&2 echo "Killed"
        exit 1
      fi

      # can't update barcodes in freyja 1.3.2; will update known issue is closed (https://github.com/andersen-lab/Freyja/issues/33)
      freyja_usher_barcode_version="freyja update: $(date +"%Y-%m-%d")"
      freyja_metadata_version="freyja update: $(date +"%Y-%m-%d")"
  else
  # configure barcode    
    if [[ ! -z "~{freyja_usher_barcodes}" ]]; then
      echo "User freyja usher barcodes identified; ~{freyja_usher_barcodes} will be utilized fre freyja demixing"
      freyja_usher_barcode_version=$(basename -- "~{freyja_usher_barcodes}")
    else
      freyja_usher_barcode_version="unmodified from freyja container: ~{docker}"  
    fi
    # configure lineage metadata
    if [[ ! -z "~{freyja_lineage_metadata}" ]]; then
      echo "User lineage metadata; ~{freyja_lineage_metadata} will be utilized fre freyja demixing"
      freyja_metadata_version=$(basename -- "~{freyja_lineage_metadata}")
    else
      freyja_metadata_version="unmodified from freyja container: ~{docker}"
    fi
  fi
  # Capture reference file versions
  echo ${freyja_usher_barcode_version} | tee FREYJA_BARCODES
  echo ${freyja_metadata_version} | tee FREYJA_METADATA
  # Call variants and capture sequencing depth information
  echo "Running: freyja variants ~{primer_trimmed_bam} --variants ~{samplename}_freyja_variants.tsv --depths ~{samplename}_freyja_depths.tsv --ref ~{reference_genome}"
  freyja variants \
    ~{primer_trimmed_bam} \
    --variants ~{samplename}_freyja_variants.tsv \
    --depths ~{samplename}_freyja_depths.tsv \
    --ref ~{reference_genome}
  # Calculate Boostraps, if specified
  if ~{bootstrap}; then
    freyja boot \
    ~{"--eps " + eps} \
    ~{"--meta " + freyja_lineage_metadata} \
    ~{"--barcodes " + freyja_usher_barcodes} \
    ~{"--nb " + number_bootstraps } \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output_base ~{samplename} \
    --boxplot pdf
  fi
  # Demix variants 
  echo "Running: freyja demix --eps ~{eps} ${freyja_barcode} ${freyja_metadata} ~{samplename}_freyja_variants.tsv ~{samplename}_freyja_depths.tsv --output ~{samplename}_freyja_demixed.tmp"
  freyja demix \
    ~{'--eps ' + eps} \
    ~{'--meta ' + freyja_lineage_metadata} \
    ~{'--barcodes ' + freyja_usher_barcodes} \
    ~{true='--confirmedonly' false='' confirmed_only} \
    ~{samplename}_freyja_variants.tsv \
    ~{samplename}_freyja_depths.tsv \
    --output ~{samplename}_freyja_demixed.tmp
  # Adjust output header
  echo -e "\t/~{samplename}" > ~{samplename}_freyja_demixed.tsv
  tail -n+2 ~{samplename}_freyja_demixed.tmp >> ~{samplename}_freyja_demixed.tsv
  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: 2
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    maxRetries: 3
  }
  output {
    File freyja_variants = "~{samplename}_freyja_variants.tsv"
    File freyja_depths = "~{samplename}_freyja_depths.tsv"
    File freyja_demixed = "~{samplename}_freyja_demixed.tsv"
    File? freyja_update_log = "freyja_update.log"
    File? freyja_bootstrap_lineages = "~{samplename}_lineages.csv"
    File? freyja_bootstrap_lineages_pdf = "~{samplename}_lineages.pdf"
    File? freyja_bootstrap_summary = "~{samplename}_summarized.csv"
    File? freyja_bootstrap_summary_pdf = "~{samplename}_summarized.pdf"
    String freyja_barcode_version = read_string("FREYJA_BARCODES")
    String freyja_metadata_version = read_string("FREYJA_METADATA")
    String freyja_version = read_string("FREYJA_VERSION")
  }
}
