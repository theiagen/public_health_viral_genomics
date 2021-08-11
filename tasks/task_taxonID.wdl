version 1.0

task kraken2 {
  input {
  	File        read1
	  File? 		  read2
	  String      samplename
	  String?     kraken2_db = "/kraken2-db"
    Int?        cpus=4
  }

  command{
    # date and version control
    date | tee DATE
    kraken2 --version | head -n1 | tee VERSION
    num_reads=$(ls *fastq.gz 2> /dev/nul | wc -l)
    if ! [ -z ${read2} ]; then
      mode="--paired"
    fi
    echo $mode
    kraken2 $mode \
      --threads ${cpus} \
      --db ${kraken2_db} \
      ${read1} ${read2} \
      --report ${samplename}_kraken2_report.txt >/dev/null

    percentage_human=$(grep "Homo sapiens" ${samplename}_kraken2_report.txt | cut -f 1)
     # | tee PERCENT_HUMAN
    percentage_sc2=$(grep "Severe acute respiratory syndrome coronavirus 2" ${samplename}_kraken2_report.txt | cut -f1 )
     # | tee PERCENT_COV
    if [ -z "$percentage_human" ] ; then percentage_human="0" ; fi
    if [ -z "$percentage_sc2" ] ; then percentage_sc2="0" ; fi
    echo $percentage_human | tee PERCENT_HUMAN
    echo $percentage_sc2 | tee PERCENT_SC2
  }

  output {
    String     date          = read_string("DATE")
    String     version       = read_string("VERSION")
    File 	     kraken_report = "${samplename}_kraken2_report.txt"
    Float 	   percent_human = read_string("PERCENT_HUMAN")
    Float 	   percent_sc2   = read_string("PERCENT_SC2")
  }

  runtime {
    docker:       "staphb/kraken2:2.0.8-beta_hv"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task pangolin {
  input {
    File        fasta
    String      samplename

  }

  command{
    # date and version control
    date | tee DATE
    pangolin --version | head -n1 | tee VERSION

    pangolin --outdir ${samplename} ${fasta}
    pangolin_lineage=$(tail -n 1 ${samplename}/lineage_report.csv | cut -f 2 -d "," | grep -v "lineage")

    pangolin_aLRT=$(tail -n 1 ${samplename}/lineage_report.csv | cut -f 3 -d "," )
    pangolin_stats=$(tail -n 1 ${samplename}/lineage_report.csv | cut -f 4 -d "," )
    mv ${samplename}/lineage_report.csv ${samplename}_pango_lineage.csv

    echo $pangolin_lineage | tee PANGOLIN_LINEAGE
    echo $pangolin_aLRT | tee PANGOLIN_aLRT
    echo $pangolin_stats | tee PANGOLIN_STATS
  }

  output {
    String     date                 = read_string("DATE")
    String     version              = read_string("VERSION")
    String     pangolin_lineage     = read_string("PANGOLIN_LINEAGE")
    Float      pangolin_aLRT        = read_string("PANGOLIN_aLRT")
    Float      pangolin_stats       = read_string("PANGOLIN_STATS")
    File       pango_lineage_report = "${samplename}_pango_lineage.csv"
  }

  runtime {
    docker:       "staphb/pangolin:1.1.14"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task pangolin2 {
  input {
    File        fasta
    String      samplename
    Int         min_length=10000
    Float       max_ambig=0.5
    String      docker
  }

  command <<<
    # date and version control
    date | tee DATE
    set -e

    pangolin "~{fasta}" \
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
        with open("VERSION", 'wt') as lineage:
          pangolin_version=line["pangolin_version"]
          pangoLEARN_version=line["pangoLEARN_version"]
          lineage.write(f"pangolin {pangolin_version}; pangoLEARN {pangoLEARN_version}")
        with open("PANGOLIN_LINEAGE", 'wt') as lineage:
          lineage.write(line["lineage"])
        with open("PANGOLIN_CONFLICTS", 'wt') as lineage:
          lineage.write(line["conflict"])
        with open("PANGOLIN_NOTES", 'wt') as lineage:
          lineage.write(line["note"])
    CODE

  >>>

  output {
    String     date                 = read_string("DATE")
    String     version              = read_string("VERSION")
    String     pangolin_lineage     = read_string("PANGOLIN_LINEAGE")
    String     pangolin_conflicts    = read_string("PANGOLIN_CONFLICTS")
    String     pangolin_notes       = read_string("PANGOLIN_NOTES")
    String     pangolin_docker      = docker
    File       pango_lineage_report = "${samplename}.pangolin_report.csv"
  }

  runtime {
    docker:     "~{docker}"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task pangolin3 {
  input {
    File        fasta
    String      samplename
    Int         min_length=10000
    Float       max_ambig=0.5
    String      docker
    String      inference_engine="usher"
  }

  command <<<
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
    # date and version control
    date | tee DATE
    conda list -n pangolin | grep "usher" | awk -F ' +' '{print$1, $2}'| tee PANGO_USHER_VERSION 
    set -e

    echo "pangolin ~{fasta} ${pango_inference}  --outfile ~{samplename}.pangolin_report.csv  --min-length ~{min_length} --max-ambig ~{max_ambig} --verbose"

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
        with open("VERSION", 'wt') as lineage:
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
    String     date                 = read_string("DATE")
    String     version              = read_string("VERSION")
    String     pangolin_lineage     = read_string("PANGOLIN_LINEAGE")
    String     pangolin_conflicts    = read_string("PANGOLIN_CONFLICTS")
    String     pangolin_notes       = read_string("PANGOLIN_NOTES")
    String     pangolin_usher_version = read_string("PANGO_USHER_VERSION")
    String     pangolin_docker      = docker
    File       pango_lineage_report = "${samplename}.pangolin_report.csv"
  }

  runtime {
    docker:     "~{docker}"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}
task pangolin_update_log {
  input {
    String samplename        
    String current_lineage
    String current_pangolin_docker
    String current_pangolin_version
    String updated_lineage
    String updated_pangolin_docker
    String updated_pangolin_version
    String? timezone
    File?  lineage_log
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
      mv "~{lineage_log}" ${lineage_log_file}
     else 
       echo "Creating new lineage log file as none was provided"
       echo -e "analysis_date\tmodification_status\tprevious_lineage\tprevious_pangolin_docker\tprevious_pangolin_version\tupdated_lineage\tupdated_pangolin_docker\tupdated_pangolin_version" > ${lineage_log_file}
     fi
     
     #populate lineage log file
     echo -e "${DATE}\t${UPDATE_STATUS}\t~{current_lineage}\t~{current_pangolin_docker}\t~{current_pangolin_version}\t~{updated_lineage}\t~{updated_pangolin_docker}\t~{updated_pangolin_version}" >> "${lineage_log_file}"
     
    echo "${UPDATE_STATUS} (${DATE})"  | tee PANGOLIN_UPDATE       

  >>>

  output {
    String     pangolin_updates = read_string("PANGOLIN_UPDATE")
    File       pango_lineage_log = "~{samplename}_pango_lineage_log.tsv"
  }

  runtime {
    docker:     "theiagen/utility:1.1"
    memory:       "8 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}

task nextclade_one_sample {
    meta {
        description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        File   genome_fasta
        File?  root_sequence
        File?  auspice_reference_tree_json
        File?  qc_config_json
        File?  gene_annotations_json
        File?  pcr_primers_csv
        String docker = "nextstrain/nextclade:0.14.4"
    }
    String basename = basename(genome_fasta, ".fasta")
    command {
        set -e
        nextclade.js --version > VERSION
        nextclade.js \
            --input-fasta "~{genome_fasta}" \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json
        cp "~{basename}".nextclade.tsv input.tsv
        
       
        python3 <<CODE
        # transpose table
        with open('input.tsv', 'r', encoding='utf-8') as inf:
            with open('transposed.tsv', 'w', encoding='utf-8') as outf:
                for c in zip(*(l.rstrip().split('\t') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
        
        # set output files as NA to ensure task doesn't fail if no relevant outputs available in Nextclade report
        echo "NA" | tee NEXTCLADE_CLADE NEXTCLADE_AASUBS NEXTCLADE_AADELS
        
        # parse transposed report file if relevant outputs are available
        if [[ $(wc -l ~{basename}.nextclade.tsv) -ge 1 ]]
        then
          grep ^clade transposed.tsv | cut -f 2 | grep -v clade > NEXTCLADE_CLADE
          grep ^aaSubstitutions transposed.tsv | cut -f 2 | grep -v aaSubstitutions | sed 's/,/|/g' > NEXTCLADE_AASUBS
          grep ^aaDeletions transposed.tsv | cut -f 2 | grep -v aaDeletions | sed 's/,/|/g' > NEXTCLADE_AADELS
        fi
    }
    runtime {
        docker: "~{docker}"
        memory: "3 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
	maxRetries:   3
    }
    output {
        String nextclade_version  = read_string("VERSION")
        File   nextclade_json     = "~{basename}.nextclade.json"
        File   auspice_json       = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv      = "~{basename}.nextclade.tsv"
        String nextclade_clade    = read_string("NEXTCLADE_CLADE")
        String nextclade_aa_subs  = read_string("NEXTCLADE_AASUBS")
        String nextclade_aa_dels  = read_string("NEXTCLADE_AADELS")
    }
}
