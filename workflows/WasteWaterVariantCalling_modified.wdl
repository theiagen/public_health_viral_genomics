version 1.0

workflow WasteWaterVariantCalling {
  meta {
      description: "Modified version of the CDPHE's WasteWaterVariantCalling WDL Worfklow to performs variant calling on SARS-CoV-2 in waster water samples and identifies mutations in the Spike gene associated with known VOCs and VUIs: https://github.com/CDPHE/WasteWaterVariantCalling)."
      author: "Kevin Libuit"
      email:  "kevin.libuit@theiagen.com"
  }
    input {
        Array[File] sorted_bam
        File covid_genome
        File spike_bed
        File spike_annotations
        Array[String] sample_id
    }

    scatter (id_bam in zip(sample_id, sorted_bam)) {
        call add_RG {
            input:
                sample_id = id_bam.left,
                bam = id_bam.right
        }
        call variant_calling {
            input:
                bam = add_RG.rgbam,
                ref = covid_genome,
                sample_id = id_bam.left
        }
        call sort_vcf {
            input:
                vcf = variant_calling.vcf,
                sample_id = id_bam.left
        }
        call sample_spike {
            input:
                vcf = sort_vcf.sorted_vcf,
                bed = spike_bed,
                sample_id = id_bam.left
        }
        call vcf2tsv {
            input:
                vcf = sample_spike.sample_spike_vcf,
                sample_id = id_bam.left,
                bed = spike_bed
        }
        call fill_NA {
            input:
                tsv = vcf2tsv.sample_spike_tsv,
                sample_id = id_bam.left,
                spike_bed = spike_bed
        }
        call allele_freq {
            input:
                tsv = fill_NA.fill_NA_tsv,
                sample_id = id_bam.left
        }
        call reformat_tsv {
            input:
                tsv = allele_freq.allele_freq_tsv,
                sample_id = id_bam.left
        }
        call summary_prep {
            input:
                tsv = reformat_tsv.reformat_tsv_tsv,
                sample_id = id_bam.left,
                spike_annotations = spike_annotations
        }
    }
    call dashboard_tsv {
        input:
            tsv = summary_prep.sample_spike_tsv_summary,
            tsv_dash = summary_prep.sample_spike_tsv_dash,
            tsv_counts = summary_prep.sample_spike_tsv_counts,
            spike_annotations = spike_annotations
    }
    call summary_tsv {
        input:
            tsv = dashboard_tsv.spike_summary_temp
    }
    
    output {
        Array[File] addrg_bam = add_RG.rgbam
        Array[File] variants = variant_calling.vcf
        Array[File] sorted_vcf = sort_vcf.sorted_vcf
        Array[File] sample_spike_vcf = sample_spike.sample_spike_vcf
        Array[File] sample_spike_tsv = vcf2tsv.sample_spike_tsv
        Array[File] sample_spike_tsv_summary = summary_prep.sample_spike_tsv_summary
        Array[File] sample_spike_tsv_dash = summary_prep.sample_spike_tsv_dash
        Array[File] fill_NA_tsv = fill_NA.fill_NA_tsv
        Array[File] allele_freq_tsv = allele_freq.allele_freq_tsv
        Array[File] reformat_tsv_tsv = reformat_tsv.reformat_tsv_tsv
        Array[File] sample_spike_tsv_counts = summary_prep.sample_spike_tsv_counts
        File spike_summary_temp = dashboard_tsv.spike_summary_temp
        File spike_summary = summary_tsv.spike_summary
        File spike_dashboard = dashboard_tsv.spike_dashboard
        File spike_counts = dashboard_tsv.spike_counts
        
    }
}

task add_RG {
    input {
        String sample_id
        File bam
    }

    command <<<
                        
        samtools addreplacerg -r ID:~{sample_id} -r LB:L1 -r SM:~{sample_id} -o ~{sample_id}_addRG.bam ~{bam}
                        
    >>>

    output {
        File rgbam = "${sample_id}_addRG.bam"
    }

    runtime {
        docker: "staphb/samtools:1.10"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}

task variant_calling {
    input {
        String sample_id
        File bam
        File ref

    }

    command <<<
        
        freebayes -f ~{ref} --haplotype-length 0 --min-alternate-count 3 --min-alternate-fraction 0.05 --min-mapping-quality 20 --min-base-quality 20 --min-coverage 10 --use-duplicate-reads --report-monomorphic --pooled-continuous ~{bam} > ~{sample_id}_variants.vcf
        
    >>>

    output {
        File vcf = "${sample_id}_variants.vcf"
    }

    runtime {
        docker: "wgspipeline/freebayes:v0.0.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 200 SSD"
    }
}

task sort_vcf {
    input {
        String sample_id
        File vcf

    }

    command <<<
        
        bgzip -c ~{vcf} > ~{sample_id}_variants.vcf.gz
        
        tabix -p vcf ~{sample_id}_variants.vcf.gz
        
        bcftools sort -O z ~{sample_id}_variants.vcf.gz > ~{sample_id}_sorted.vcf.gz
        
    >>>

    output {
        File sorted_vcf = "${sample_id}_sorted.vcf.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }
}

task sample_spike {
    input {
        File vcf
        File bed
        String sample_id
    }

    command <<<
        
        tabix -p vcf ~{vcf}
        
        bcftools view --regions-file ~{bed} --output-type v --output-file ~{sample_id}_spike_mutations.vcf ~{vcf}
        
    >>>

    output {
         File sample_spike_vcf = "${sample_id}_spike_mutations.vcf"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task vcf2tsv {
    input {
        File vcf
        File bed
        String sample_id
    }

    command <<<
    
        bgzip -c ~{vcf} > ~{sample_id}_spike_mutations.vcf.gz
        
        tabix -p vcf ~{sample_id}_spike_mutations.vcf.gz
        
        bcftools query --regions-file ~{bed} --format '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%RO\t%AO]\n' ~{sample_id}_spike_mutations.vcf.gz > ~{sample_id}_spike_mutations.tsv
        
    >>>

    output {
         File sample_spike_tsv = "${sample_id}_spike_mutations.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

task fill_NA {
    input {
        File tsv
        String sample_id
        File spike_bed
    }

    command <<<    
        
        # create key of unique locations
        cat ~{spike_bed} | cut -f 1,2 | tr "\t" "_" | sort | uniq > keys.txt
        
        # add headers to tsv and use key to fill in missing values
        echo -e "CHROM\tPOS\tREF\t~{sample_id}_ALT\t~{sample_id}_DP\t~{sample_id}_RO\t~{sample_id}_AO" | cat - ~{tsv} | sed 's/\t/_/' | sort -t $'\t' -k1,1 > ~{sample_id}_spike_mutations_temp1.tsv
        
        # get the filled columns we want
        join -t $'\t' -e NA -a 1 -1 1 -2 1 -o "1.1,2.3,2.4,2.6" keys.txt "~{sample_id}_spike_mutations_temp1.tsv" > ~{sample_id}_spike_fill_NA.tsv
            
    >>>

    output {
         File fill_NA_tsv = "${sample_id}_spike_fill_NA.tsv"
    }

    runtime {
        docker: "quay.io/theiagen/utility:1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task allele_freq {
    input {
        File tsv
        String sample_id
    }

    command <<<    
        
        # separate the comma separated alleles into separate rows (might need to fix delimiters)
        awk '{split($2,a,","); split($4,b,","); for(i in a){print $1,a[i],$3,b[i]}}' ~{tsv} > ~{sample_id}_spike_mutations_temp2.tsv
        
        # use AO and DP fields to calculate ALT allele frequency, fix delimiters, change -nan allele frequencies to NA
        awk '$3~"^NA"||$4~"^NA"{$5="NA";print;next}{$5=$4/$3}1' ~{sample_id}_spike_mutations_temp2.tsv | sed 's/ /\t/g' | awk '$5 == "-nan" {$5="NA"} 1' OFS="\t" > ~{sample_id}_spike_allele_freq.tsv
    
    >>>

    output {
         File allele_freq_tsv = "${sample_id}_spike_allele_freq.tsv"
    }

    runtime {
        docker: "quay.io/theiagen/utility:1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task reformat_tsv {
    input {
        File tsv
        String sample_id
    }

    command <<<    

        # combine the rows based on matching nucl location
        
        awk '{f2[$1]=f2[$1] sep[$1] $2; 
              f3[$1]=f3[$1] sep[$1] $3;
              f4[$1]=f4[$1] sep[$1] $4;
              f5[$1]=f5[$1] sep[$1] $5; 
              sep[$1]=";"}
         END {for(k in f2) print k,f2[k],f3[k],f4[k],f5[k]}' ~{tsv} > ~{sample_id}_spike_mutations_temp3.tsv
         
        # fix delimiters, add a column containing the sample ids
        sed 's/ /\t/g' ~{sample_id}_spike_mutations_temp3.tsv | awk 'NF=NF+1{$NF="~{sample_id}"}1' > ~{sample_id}_spike_mutations_temp4.tsv
        
        # fix the column headers, convert from space to tab delimited and then sort by col1
        echo -e "CHROMPOS ~{sample_id}_ALT ~{sample_id}_DP ~{sample_id}_AO ~{sample_id}_ALTfreq sample_id" | cat - ~{sample_id}_spike_mutations_temp4.tsv | sed 's/ /\t/g' | sort -t $'\t' -k 1,1 -V > ~{sample_id}_spike_reformat.tsv
        
    >>>

    output {
         File reformat_tsv_tsv = "${sample_id}_spike_reformat.tsv"
    }

    runtime {
        docker: "quay.io/theiagen/utility:1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task summary_prep {
    input {
        File tsv
        String sample_id
        File spike_annotations
    }

    command <<<    
        
        # cut the columns we want for the results summary and make output file
        cut -f2,5 ~{tsv} > ~{sample_id}_spike_mutations_forsummary.tsv
        
        # cut the columns we want for the dashboard summary
        awk '{print $6 "\t" $2 "\t" $5}' ~{tsv} > ~{sample_id}_spike_mutations_temp5.tsv
        
        # add annotations to the dashboard summary, reorder the dashboard summary columns, fix the dashboard summary column headers and make output file
        paste ~{spike_annotations} ~{sample_id}_spike_mutations_temp5.tsv | awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6}' | awk 'BEGIN{FS=OFS="\t"; print "sample_id", "AA_change", "Nucl_change", "Lineages", "ALT", "ALTfreq"} NR>1{print $1, $2, $3, $4, $5, $6}' > ~{sample_id}_spike_mutations_fordash.tsv
    
        # cut the columns we want for the counts summary
        awk '{print $6 "\t" $2 "\t" $3 "\t" $4}' ~{tsv} > ~{sample_id}_spike_mutations_temp6.tsv
        
        # add annotations to the counts summary, reorder the dashboard summary columns, fix the dashboard summary column headers and make output file
        paste ~{spike_annotations} ~{sample_id}_spike_mutations_temp6.tsv | awk '{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7}' | awk 'BEGIN{FS=OFS="\t"; print "sample_id", "AA_change", "Nucl_change", "Lineages", "ALT", "Total_count", "ALT_count"} NR>1{print $1, $2, $3, $4, $5, $6, $7}' > ~{sample_id}_spike_mutations_counts.tsv
   
   >>>

    output {
         File sample_spike_tsv_summary = "${sample_id}_spike_mutations_forsummary.tsv"
         File sample_spike_tsv_dash = "${sample_id}_spike_mutations_fordash.tsv"
         File sample_spike_tsv_counts = "${sample_id}_spike_mutations_counts.tsv"
    }

    runtime {
        docker: "quay.io/theiagen/utility:1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 2500 HDD"
    }
}

task dashboard_tsv {
    input {
        Array[File] tsv
        Array[File] tsv_dash
        Array[File] tsv_counts
        File spike_annotations
    }

    command <<<
        
        # concatenate the tsvs and make the dashboard summary output
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' tsv_dash} >> spike_mutations_dashboard.tsv
        
        # concatenate the tsvs and make the dashboard summary output
        awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=' ' tsv_counts} >> spike_mutations_counts.tsv
        
        # fix delimiters in annotations file
        sed 's/ /\t/g' ~{spike_annotations} > spike_annotations.tsv
        
        # concatentate tsvs for sequencing and bioinformatics team summary file and make output
        paste spike_annotations.tsv ~{sep=' ' tsv} > spike_mutations_summary_temp.tsv

    >>>

    output {
        File spike_summary_temp = "spike_mutations_summary_temp.tsv"
        File spike_dashboard = "spike_mutations_dashboard.tsv"
        File spike_counts = "spike_mutations_counts.tsv"
    }

    runtime {
        docker: "quay.io/theiagen/utility:1.1"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}

task summary_tsv {
    input {
        File tsv
    }

    command <<<
        
        # datamash to tranpose results summary
        datamash -H transpose < ~{tsv} > spike_mutations_summary.tsv

    >>>

    output {
        File spike_summary = "spike_mutations_summary.tsv"
    }

    runtime {
        docker: "rapatsky/debian"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}
