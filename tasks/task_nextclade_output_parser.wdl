version 1.0

task nextclade_output_parser_one_sample {
    meta {
        description: "Python and bash codeblocks for parsing the output files from Nextclade."
    }
    input {
        File   nextclade_tsv
        String docker = "quay.io/theiagen/utility:1.1"
    }
    command {
      python3 <<CODE
      # transpose table
      with open(~{nextclade_tsv}, 'r', encoding='utf-8') as inf:
          with open('transposed.tsv', 'w', encoding='utf-8') as outf:
              for c in zip(*(l.rstrip().split('\t') for l in inf)):
                  outf.write('\t'.join(c)+'\n')
      CODE

      # set output files as NA to ensure task doesn't fail if no relevant outputs available in Nextclade report
      echo "NA" | tee NEXTCLADE_CLADE NEXTCLADE_AASUBS NEXTCLADE_AADELS

      # parse transposed report file if relevant outputs are available
      if [[ $(wc -l ~{nextclade_tsv}) -ge 1 ]]
      then
        grep ^clade transposed.tsv | cut -f 2 | grep -v clade > NEXTCLADE_CLADE
        grep ^aaSubstitutions transposed.tsv | cut -f 2 | grep -v aaSubstitutions | sed 's/,/|/g' > NEXTCLADE_AASUBS
        grep ^aaDeletions transposed.tsv | cut -f 2 | grep -v aaDeletions | sed 's/,/|/g' > NEXTCLADE_AADELS
      fi
    }
    runtime {
        docker: "~{docker}"
        memory: "4 GB"
        cpu:    2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        String nextclade_clade    = read_string("NEXTCLADE_CLADE")
        String nextclade_aa_subs  = read_string("NEXTCLADE_AASUBS")
        String nextclade_aa_dels  = read_string("NEXTCLADE_AADELS")
    }
}
