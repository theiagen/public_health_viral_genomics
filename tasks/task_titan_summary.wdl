version 1.0

task titan_summary {
    input {
        String samplename
        String seq_platform
        Float  percent_reference_coverage
        Int    number_N
        String pango_lineage
        String pangolin_conflicts
        String pangolin_notes
        String pangolin_version
        String pangolin_docker
        String nextclade_clade
        String nextclade_aa_subs
        String nextclade_aa_dels
        Int    vadr_num_alerts
        Int    assembly_length_unambiguous
        Float  assembly_mean_coverage
        String assembly_method
        Int    number_Degenerate
        Int    number_Total
        Float  meanbaseq_trim
        Float  meanmapq_trim
        Int    fastqc_clean1
        Int    fastqc_clean2 = ""
        String fastqc_clean_pairs = ""
        Int    fastqc_raw1
        Int    fastqc_raw2 = ""
        String fastqc_raw_pairs = ""
        Float  kraken_human
        Float  kraken_human_dehosted = ""
        String kraken_report
        String kraken_report_dehosted = ""
        Float  kraken_sc2
        Float  kraken_sc2_dehosted = ""
        Float  pool1_percent = ""
        Float  pool2_percent = ""
        Float  primer_trimmed_read_percent
        String artic_version
        String bbduk_docker = ""
        String bedtools_version = ""
        String bwa_version = ""
        String fastqc_version = ""
        String ivar_variant_version = ""
        String ivar_version_consensus = ""
        String ivar_version_primtrim = ""
        String kraken_version = ""
        String nextclade_version = ""
        String samtools_version
        String samtools_version_consensus = ""
        String samtools_version_primtrim = ""
        String samtools_version_stats = ""
        String trimmomatic_version = ""
        String vadr_docker = ""
    }

    command <<<
        python3<<CODE
        data = [
            '~{samplename}',
            '~{seq_platform}',
            '~{percent_reference_coverage}',
            '~{number_N}',
            '~{pango_lineage}',
            '~{pangolin_conflicts}',
            '~{pangolin_notes}',
            '~{pangolin_version}',
            '~{pangolin_docker}',
            '~{nextclade_clade}',
            '~{nextclade_aa_subs}',
            '~{nextclade_aa_dels}',
            '~{vadr_num_alerts}',
            '~{assembly_length_unambiguous}',
            '~{assembly_mean_coverage}',
            '~{assembly_method}',
            '~{number_Degenerate}',
            '~{number_Total}',
            '~{meanbaseq_trim}',
            '~{meanmapq_trim}',
            '~{fastqc_clean1}',
            '~{fastqc_clean2}',
            '~{fastqc_clean_pairs}',
            '~{fastqc_raw1}',
            '~{fastqc_raw2}',
            '~{fastqc_raw_pairs}',
            '~{kraken_human}',
            '~{kraken_human_dehosted}',
            '~{kraken_report}',
            '~{kraken_report_dehosted}',
            '~{kraken_sc2}',
            '~{kraken_sc2_dehosted}',
            '~{pool1_percent}',
            '~{pool2_percent}',
            '~{primer_trimmed_read_percent}',
            '~{artic_version}',
            '~{bbduk_docker}',
            '~{bwa_version}',
            '~{fastqc_version}',
            '~{ivar_variant_version}',
            '~{ivar_version_consensus}',
            '~{ivar_version_primtrim}',
            '~{kraken_version}',
            '~{nextclade_version}',
            '~{samtools_version}',
            '~{samtools_version_consensus}',
            '~{samtools_version_primtrim}',
            '~{samtools_version_stats}',
            '~{trimmomatic_version}',
            '~{vadr_docker}'
        ]
        print("\t".join(data))
        CODE
    >>>

    output {
        String summary = read_string(stdout())
    }

    runtime {
        docker: "python:slim"
        memory: "1 GB"
        cpu: 1
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task merge_titan_summary {
    input {
        Array[String] summaries
    }

    command <<<
        python3<<CODE
        with open("titan-results.tsv", 'wt') as tsv_fh:
            header = [
                'samplename',
                'seq_platform',
                'percent_reference_coverage',
                'number_N',
                'pango_lineage',
                'pangolin_conflicts',
                'pangolin_notes',
                'pangolin_version',
                'pangolin_docker',
                'nextclade_clade',
                'nextclade_aa_subs',
                'nextclade_aa_dels',
                'vadr_num_alerts',
                'assembly_length_unambiguous',
                'assembly_mean_coverage',
                'assembly_method',
                'number_Degenerate',
                'number_Total',
                'meanbaseq_trim',
                'meanmapq_trim',
                'fastqc_clean1',
                'fastqc_clean2',
                'fastqc_clean_pairs',
                'fastqc_raw1',
                'fastqc_raw2',
                'fastqc_raw_pairs',
                'kraken_human',
                'kraken_human_dehosted',
                'kraken_report',
                'kraken_report_dehosted',
                'kraken_sc2',
                'kraken_sc2_dehosted',
                'pool1_percent',
                'pool2_percent',
                'primer_trimmed_read_percent',
                'artic_version',
                'bbduk_docker',
                'bwa_version',
                'fastqc_version',
                'ivar_variant_version',
                'ivar_version_consensus',
                'ivar_version_primtrim',
                'kraken_version',
                'nextclade_version',
                'samtools_version',
                'samtools_version_consensus',
                'samtools_version_primtrim',
                'samtools_version_stats',
                'trimmomatic_version',
                'vadr_docker'
            ]
            header = "\t".join(header)
            tsv_fh.write(f'{header}\n')
            summaries = '~{sep="SUMMARY" summaries}'.split('SUMMARY')
            for summary in summaries:
                tsv_fh.write(f'{summary}\n')
        CODE
    >>>

    output {
        File merged_summaries = "titan-results.tsv"
    }

    runtime {
        docker: "python:slim"
        memory: "1 GB"
        cpu: 1
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}
