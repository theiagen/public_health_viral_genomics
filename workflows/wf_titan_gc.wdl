version 1.0

# Workflows from Theiagen's public_health_viral_genomics
# Source: https://github.com/theiagen/public_health_viral_genomics
import "wf_titan_clearlabs.wdl" as clearlabs
import "wf_titan_illumina_pe.wdl" as illumina_pe
import "wf_titan_illumina_se.wdl" as illumina_se
import "wf_titan_ont.wdl" as ont
import "../tasks/task_titan_summary.wdl" as summary

struct parseJSON {
    String samplename
    String run_id
    String platform
    File   r1
    File   r2
    File   primer_bed
}

workflow titan_gc {
    meta {
        description: "Incorporates each of the Titan workflows (clearlabs, illumina_pe, illumina_se, ont) into a single run."
        author: "Robert A. Petit III"
        email:  "robert.petit@theiagen.com"
    }

    input {
        Array[parseJSON] samples
        String  pangolin_docker_image = "staphb/pangolin:3.1.3-pangolearn-2021-06-15"
    }

    scatter (sample in samples) {
        if (sample.platform == "clearlabs") {
            call clearlabs.titan_clearlabs as titan_clearlabs { 
                input:
                    samplename = sample.samplename,
                    clear_lab_fastq = sample.r1,
                    primer_bed = sample.primer_bed,
                    pangolin_docker_image = pangolin_docker_image
            }

           call summary.titan_summary as clearlabs_summary {
                input:
                    samplename = sample.samplename,
                    titan_workflow = 'titan_clearlabs',
                    titan_version = titan_clearlabs.titan_clearlabs_version,
                    titan_analysis_date = titan_clearlabs.titan_clearlabs_analysis_date,
                    seq_platform = titan_clearlabs.seq_platform,
                    primer_bed_name = titan_clearlabs.primer_bed_name,
                    percent_reference_coverage = titan_clearlabs.percent_reference_coverage,
                    number_N = titan_clearlabs.number_N,
                    pango_lineage = titan_clearlabs.pango_lineage,
                    pangolin_conflicts = titan_clearlabs.pangolin_conflicts,
                    pangolin_notes = titan_clearlabs.pangolin_notes,
                    pangolin_version = titan_clearlabs.pangolin_version,
                    pangolin_docker = titan_clearlabs.pangolin_docker,
                    pangolin_usher_version = titan_clearlabs.pangolin_usher_version,
                    nextclade_clade = titan_clearlabs.nextclade_clade,
                    nextclade_aa_subs = titan_clearlabs.nextclade_aa_subs,
                    nextclade_aa_dels = titan_clearlabs.nextclade_aa_dels,
                    vadr_num_alerts = titan_clearlabs.vadr_num_alerts,
                    assembly_length_unambiguous = titan_clearlabs.assembly_length_unambiguous,
                    assembly_mean_coverage = titan_clearlabs.assembly_mean_coverage,
                    assembly_method = titan_clearlabs.assembly_method,
                    number_Degenerate = titan_clearlabs.number_Degenerate,
                    number_Total = titan_clearlabs.number_Total,
                    meanbaseq_trim = titan_clearlabs.meanbaseq_trim,
                    meanmapq_trim = titan_clearlabs.meanmapq_trim,
                    fastqc_clean1 = titan_clearlabs.fastqc_clean,
                    fastqc_raw1 = titan_clearlabs.fastqc_raw,
                    kraken_human = titan_clearlabs.kraken_human,
                    kraken_human_dehosted = titan_clearlabs.kraken_human_dehosted,
                    kraken_sc2 = titan_clearlabs.kraken_sc2,
                    kraken_sc2_dehosted = titan_clearlabs.kraken_sc2_dehosted,
                    artic_version = titan_clearlabs.artic_version,
                    kraken_version = titan_clearlabs.kraken_version,
                    nextclade_version = titan_clearlabs.nextclade_version,
                    samtools_version = titan_clearlabs.samtools_version,
                    vadr_docker = titan_clearlabs.vadr_docker
            }
        }

        if (sample.platform == "illumina_pe") {
            call illumina_pe.titan_illumina_pe as titan_illumina_pe { 
                input:
                    samplename = sample.samplename,
                    read1_raw = sample.r1,
                    read2_raw = sample.r2,
                    primer_bed = sample.primer_bed,
                    pangolin_docker_image = pangolin_docker_image
            }

           call summary.titan_summary as illumina_pe_summary {
                input:
                    samplename = sample.samplename,
                    titan_workflow = 'titan_illumina_pe',
                    titan_version = titan_illumina_pe.titan_illumina_pe_version,
                    titan_analysis_date = titan_illumina_pe.titan_illumina_pe_analysis_date,
                    seq_platform = titan_illumina_pe.seq_platform,
                    primer_bed_name = titan_illumina_pe.primer_bed_name,
                    percent_reference_coverage = titan_illumina_pe.percent_reference_coverage,
                    number_N = titan_illumina_pe.number_N,
                    pango_lineage = titan_illumina_pe.pango_lineage,
                    pangolin_conflicts = titan_illumina_pe.pangolin_conflicts,
                    pangolin_notes = titan_illumina_pe.pangolin_notes,
                    pangolin_version = titan_illumina_pe.pangolin_version,
                    pangolin_docker = titan_illumina_pe.pangolin_docker,
                    pangolin_usher_version = titan_illumina_pe.pangolin_usher_version,
                    nextclade_clade = titan_illumina_pe.nextclade_clade,
                    nextclade_aa_subs = titan_illumina_pe.nextclade_aa_subs,
                    nextclade_aa_dels = titan_illumina_pe.nextclade_aa_dels,
                    vadr_num_alerts = titan_illumina_pe.vadr_num_alerts,
                    assembly_length_unambiguous = titan_illumina_pe.assembly_length_unambiguous,
                    assembly_mean_coverage = titan_illumina_pe.assembly_mean_coverage,
                    assembly_method = titan_illumina_pe.assembly_method,
                    number_Degenerate = titan_illumina_pe.number_Degenerate,
                    number_Total = titan_illumina_pe.number_Total,
                    meanbaseq_trim = titan_illumina_pe.meanbaseq_trim,
                    meanmapq_trim = titan_illumina_pe.meanmapq_trim,
                    fastqc_clean1 = titan_illumina_pe.fastqc_clean1,
                    fastqc_clean2 = titan_illumina_pe.fastqc_clean2,
                    fastqc_clean_pairs = titan_illumina_pe.fastqc_clean_pairs,
                    fastqc_raw1 = titan_illumina_pe.fastqc_raw1,
                    fastqc_raw2 = titan_illumina_pe.fastqc_raw2,
                    fastqc_raw_pairs = titan_illumina_pe.fastqc_raw_pairs,
                    kraken_human = titan_illumina_pe.kraken_human,
                    kraken_human_dehosted = titan_illumina_pe.kraken_human_dehosted,
                    kraken_sc2 = titan_illumina_pe.kraken_sc2,
                    kraken_sc2_dehosted = titan_illumina_pe.kraken_sc2_dehosted,
                    primer_trimmed_read_percent = titan_illumina_pe.primer_trimmed_read_percent,
                    bbduk_docker = titan_illumina_pe.bbduk_docker,
                    bwa_version = titan_illumina_pe.bwa_version,
                    fastqc_version = titan_illumina_pe.fastqc_version,
                    ivar_variant_version = titan_illumina_pe.ivar_variant_version,
                    ivar_version_consensus = titan_illumina_pe.ivar_version_consensus,
                    ivar_version_primtrim = titan_illumina_pe.ivar_version_primtrim,
                    kraken_version = titan_illumina_pe.kraken_version,
                    nextclade_version = titan_illumina_pe.nextclade_version,
                    samtools_version = titan_illumina_pe.samtools_version,
                    samtools_version_consensus = titan_illumina_pe.samtools_version_consensus,
                    samtools_version_primtrim = titan_illumina_pe.samtools_version_primtrim,
                    samtools_version_stats = titan_illumina_pe.samtools_version_stats,
                    trimmomatic_version = titan_illumina_pe.trimmomatic_version,
                    vadr_docker = titan_illumina_pe.vadr_docker
            }
        }
        
        if (sample.platform == "illumina_se") {
            call illumina_se.titan_illumina_se as titan_illumina_se { 
                input:
                    samplename = sample.samplename,
                    read1_raw  = sample.r1,
                    primer_bed = sample.primer_bed,
                    pangolin_docker_image = pangolin_docker_image
            }

            call summary.titan_summary as illumina_se_summary {
                input:
                    samplename = sample.samplename,
                    titan_workflow = 'titan_illumina_se',
                    titan_version = titan_illumina_se.titan_illumina_se_version,
                    titan_analysis_date = titan_illumina_se.titan_illumina_se_analysis_date,
                    seq_platform = titan_illumina_se.seq_platform,
                    primer_bed_name = titan_illumina_se.primer_bed_name,
                    percent_reference_coverage = titan_illumina_se.percent_reference_coverage,
                    number_N = titan_illumina_se.number_N,
                    pango_lineage = titan_illumina_se.pango_lineage,
                    pangolin_conflicts = titan_illumina_se.pangolin_conflicts,
                    pangolin_notes = titan_illumina_se.pangolin_notes,
                    pangolin_version = titan_illumina_se.pangolin_version,
                    pangolin_docker = titan_illumina_se.pangolin_docker,
                    pangolin_usher_version = titan_illumina_se.pangolin_usher_version,
                    nextclade_clade = titan_illumina_se.nextclade_clade,
                    nextclade_aa_subs = titan_illumina_se.nextclade_aa_subs,
                    nextclade_aa_dels = titan_illumina_se.nextclade_aa_dels,
                    vadr_num_alerts = titan_illumina_se.vadr_num_alerts,
                    assembly_length_unambiguous = titan_illumina_se.assembly_length_unambiguous,
                    assembly_mean_coverage = titan_illumina_se.assembly_mean_coverage,
                    assembly_method = titan_illumina_se.assembly_method,
                    number_Degenerate = titan_illumina_se.number_Degenerate,
                    number_Total = titan_illumina_se.number_Total,
                    meanbaseq_trim = titan_illumina_se.meanbaseq_trim,
                    meanmapq_trim = titan_illumina_se.meanmapq_trim,
                    fastqc_clean1 = titan_illumina_se.fastqc_clean,
                    fastqc_raw1 = titan_illumina_se.fastqc_raw,
                    kraken_human = titan_illumina_se.kraken_human,
                    kraken_sc2 = titan_illumina_se.kraken_sc2,
                    primer_trimmed_read_percent = titan_illumina_se.primer_trimmed_read_percent,
                    bbduk_docker = titan_illumina_se.bbduk_docker,
                    bwa_version = titan_illumina_se.bwa_version,
                    fastqc_version = titan_illumina_se.fastqc_version,
                    ivar_variant_version = titan_illumina_se.ivar_variant_version,
                    ivar_version_consensus = titan_illumina_se.ivar_version_consensus,
                    ivar_version_primtrim = titan_illumina_se.ivar_version_primtrim,
                    kraken_version = titan_illumina_se.kraken_version,
                    nextclade_version = titan_illumina_se.nextclade_version,
                    samtools_version = titan_illumina_se.samtools_version,
                    samtools_version_consensus = titan_illumina_se.samtools_version_consensus,
                    samtools_version_primtrim = titan_illumina_se.samtools_version_primtrim,
                    samtools_version_stats = titan_illumina_se.samtools_version_stats,
                    trimmomatic_version = titan_illumina_se.trimmomatic_version,
                    vadr_docker = titan_illumina_se.vadr_docker
            }
        }
        
        if (sample.platform == "ont") {
            call ont.titan_ont as titan_ont { 
                input:
                    samplename = sample.samplename,
                    demultiplexed_reads = sample.r1,
                    primer_bed = sample.primer_bed,
                    pangolin_docker_image = pangolin_docker_image
            }

            call summary.titan_summary as ont_summary {
                input:
                    samplename = sample.samplename,
                    titan_workflow = 'titan_ont',
                    titan_version = titan_ont.titan_ont_version,
                    titan_analysis_date = titan_ont.titan_ont_analysis_date,
                    seq_platform = titan_ont.seq_platform,
                    primer_bed_name = titan_ont.primer_bed_name,
                    percent_reference_coverage = titan_ont.percent_reference_coverage,
                    number_N = titan_ont.number_N,
                    pango_lineage = titan_ont.pango_lineage,
                    pangolin_conflicts = titan_ont.pangolin_conflicts,
                    pangolin_notes = titan_ont.pangolin_notes,
                    pangolin_version = titan_ont.pangolin_version,
                    pangolin_docker = titan_ont.pangolin_docker,
                    pangolin_usher_version = titan_ont.pangolin_usher_version,
                    nextclade_clade = titan_ont.nextclade_clade,
                    nextclade_aa_subs = titan_ont.nextclade_aa_subs,
                    nextclade_aa_dels = titan_ont.nextclade_aa_dels,
                    vadr_num_alerts = titan_ont.vadr_num_alerts,
                    assembly_length_unambiguous = titan_ont.assembly_length_unambiguous,
                    assembly_mean_coverage = titan_ont.assembly_mean_coverage,
                    assembly_method = titan_ont.assembly_method,
                    number_Degenerate = titan_ont.number_Degenerate,
                    number_Total = titan_ont.number_Total,
                    meanbaseq_trim = titan_ont.meanbaseq_trim,
                    meanmapq_trim = titan_ont.meanmapq_trim,
                    fastqc_clean1 = titan_ont.fastqc_clean,
                    fastqc_raw1 = titan_ont.fastqc_raw,
                    kraken_human = titan_ont.kraken_human,
                    kraken_human_dehosted = titan_ont.kraken_human_dehosted,
                    kraken_sc2 = titan_ont.kraken_sc2,
                    kraken_sc2_dehosted = titan_ont.kraken_sc2_dehosted,
                    artic_version = titan_ont.artic_version,
                    kraken_version = titan_ont.kraken_version,
                    nextclade_version = titan_ont.nextclade_version,
                    samtools_version = titan_ont.samtools_version,
                    vadr_docker = titan_ont.vadr_docker
            }
        }
    }

    call summary.merge_titan_summary {
        input:
            clearlabs_summaries = clearlabs_summary.summary,
            illumina_pe_summaries = illumina_pe_summary.summary,
            illumina_se_summaries = illumina_se_summary.summary,
            ont_summaries = ont_summary.summary
    }

    output {
        # Titan outputs
        File        summaries_tsv         = merge_titan_summary.summaries_tsv
        File        summaries_json        = merge_titan_summary.summaries_json
        Array[File] reads_dehosted        = flatten([
                                                select_all(titan_clearlabs.reads_dehosted), select_all(titan_illumina_pe.read1_dehosted), 
                                                select_all(titan_illumina_pe.read2_dehosted), select_all(titan_ont.reads_dehosted)
                                            ])
        Array[File] aligned_bam           = flatten([
                                                select_all(titan_clearlabs.aligned_bam), select_all(titan_illumina_pe.aligned_bam),
                                                select_all(titan_illumina_se.aligned_bam), select_all(titan_ont.aligned_bam)
                                            ])
        Array[File] aligned_bai           = flatten([
                                                select_all(titan_clearlabs.aligned_bai), select_all(titan_illumina_pe.aligned_bai),
                                                select_all(titan_illumina_se.aligned_bai), select_all(titan_ont.aligned_bai)
                                            ])
        Array[File] assembly_fasta        = flatten([
                                                select_all(titan_clearlabs.assembly_fasta), select_all(titan_illumina_pe.assembly_fasta),
                                                select_all(titan_illumina_se.assembly_fasta), select_all(titan_ont.assembly_fasta)
                                            ])
        Array[File] consensus_stats       = flatten([
                                                select_all(titan_clearlabs.consensus_stats), select_all(titan_illumina_pe.consensus_stats),
                                                select_all(titan_illumina_se.consensus_stats), select_all(titan_ont.consensus_stats)
                                            ])
        Array[File] consensus_flagstat    = flatten([
                                                select_all(titan_clearlabs.consensus_flagstat), select_all(titan_illumina_pe.consensus_flagstat),
                                                select_all(titan_illumina_se.consensus_flagstat), select_all(titan_ont.consensus_flagstat)
                                            ])
        Array[File] pango_lineage_report  = flatten([
                                                select_all(titan_clearlabs.pango_lineage_report), select_all(titan_illumina_pe.pango_lineage_report),
                                                select_all(titan_illumina_se.pango_lineage_report), select_all(titan_ont.pango_lineage_report)
                                            ])
        Array[File] nextclade_json        = flatten([
                                                select_all(titan_clearlabs.nextclade_json), select_all(titan_illumina_pe.nextclade_json),
                                                select_all(titan_illumina_se.nextclade_json), select_all(titan_ont.nextclade_json)
                                            ])
        Array[File] auspice_json          = flatten([
                                                select_all(titan_clearlabs.auspice_json), select_all(titan_illumina_pe.auspice_json),
                                                select_all(titan_illumina_se.auspice_json), select_all(titan_ont.auspice_json)
                                            ])
        Array[File] nextclade_tsv         = flatten([
                                                select_all(titan_clearlabs.nextclade_tsv), select_all(titan_illumina_pe.nextclade_tsv),
                                                select_all(titan_illumina_se.nextclade_tsv), select_all(titan_ont.nextclade_tsv)
                                            ])
        Array[File] vadr_alerts_list      = flatten([
                                                select_all(titan_clearlabs.vadr_alerts_list), select_all(titan_illumina_pe.vadr_alerts_list),
                                                select_all(titan_illumina_se.vadr_alerts_list), select_all(titan_ont.vadr_alerts_list)
                                            ])
        Array[File] kraken_report         = flatten([
                                                select_all(titan_clearlabs.kraken_report), select_all(titan_illumina_pe.kraken_report),
                                                select_all(titan_illumina_se.kraken_report), select_all(titan_ont.kraken_report)
                                            ])
        Array[File] kraken_report_dehosted = flatten([
                                                select_all(titan_clearlabs.kraken_report_dehosted), select_all(titan_illumina_pe.kraken_report_dehosted),
                                                select_all(titan_ont.kraken_report_dehosted)
                                            ])
        Array[File] json_summary           = flatten([
                                                select_all(clearlabs_summary.summary), select_all(illumina_pe_summary.summary),
                                                select_all(illumina_se_summary.summary), select_all(ont_summary.summary)
                                            ])             
    }
}
