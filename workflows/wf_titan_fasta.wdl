version 1.0

import "../tasks/task_ont_medaka.wdl" as medaka
import "../tasks/task_assembly_metrics.wdl" as assembly_metrics
import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_qc_utils.wdl" as qc_utils
import "../tasks/task_versioning.wdl" as versioning

workflow titan_clearlabs {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }

  input {
    String  samplename
    File    assembly_fasta
    String  seq_method  
    String  input_assembly_method
    String  nextclade_dataset_name = "sars-cov-2"
    String  nextclade_dataset_reference = "MN908947"
    String  nextclade_dataset_tag = "2021-06-25T00:00:00Z"
  }
  call qc_utils.consensus_qc {
    input:
      assembly_fasta = assembly_fasta
  }
  call taxon_ID.pangolin3 {
    input:
      samplename = samplename,
      fasta = assembly_fasta
  }
  call taxon_ID.nextclade_one_sample {
    input:
      genome_fasta = assembly_fasta,
      dataset_name = nextclade_dataset_name,
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
  }
  call taxon_ID.nextclade_output_parser_one_sample {
    input:
      nextclade_tsv = nextclade_one_sample.nextclade_tsv
  }
  call ncbi.vadr {
    input:
      genome_fasta = assembly_fasta,
      assembly_length_unambiguous = consensus_qc.number_ATCG
  }
  call versioning.version_capture{
    input:
  }
  output {
    String titan_fasta_version          = version_capture.phvg_version
    String titan_fasta_analysis_date    = version_capture.date
    String seq_platform                 = seq_method
    String assembly_method              = input_assembly_method
    
    Int    number_N                         = consensus_qc.number_N
    Int    assembly_length_unambiguous      = consensus_qc.number_ATCG
    Int    number_Degenerate                = consensus_qc.number_Degenerate
    Int    number_Total                     = consensus_qc.number_Total
    Float  percent_reference_coverage       = consensus_qc.percent_reference_coverage

    String pango_lineage                    = pangolin3.pangolin_lineage
    String pangolin_conflicts               = pangolin3.pangolin_conflicts
    String pangolin_notes                   = pangolin3.pangolin_notes
    String pangolin_assignment_version                 = pangolin3.pangolin_assignment_version
    File   pango_lineage_report             = pangolin3.pango_lineage_report
    String pangolin_docker                  = pangolin3.pangolin_docker
    String pangolin_versions           = pangolin3.pangolin_versions

    File   nextclade_json                   = nextclade_one_sample.nextclade_json
    File   auspice_json                     = nextclade_one_sample.auspice_json
    File   nextclade_tsv                    = nextclade_one_sample.nextclade_tsv
    String nextclade_version                = nextclade_one_sample.nextclade_version

    String nextclade_clade                  = nextclade_output_parser_one_sample.nextclade_clade
    String nextclade_aa_subs                = nextclade_output_parser_one_sample.nextclade_aa_subs
    String nextclade_aa_dels                = nextclade_output_parser_one_sample.nextclade_aa_dels

    File?  vadr_alerts_list                 = vadr.alerts_list
    String vadr_num_alerts                  = vadr.num_alerts
    String vadr_docker                      = vadr.vadr_docker
  }
}
