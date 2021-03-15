version 1.0

import "wf_sarscov2_nextstrain_noSub.wdl" as augur_noSub
import "../tasks/task_phylo.wdl" as phylo
import "../tasks/task_data_vis.wdl" as vis

workflow titan_augur_run {
    meta {
        description: "Align assemblies, build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+    assembly_fastas
        Array[File]+    sample_metadata_tsvs

        String          build_name
#        File            builds_yaml

        Array[String]?  ancestral_traits_to_infer

        File?           auspice_config
        File?           ref_fasta
        File?           clades_tsv
        File?           lat_longs_tsv
        File?           render_template

        Int             min_unambig_genome = 27000
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two. They may be compressed (gz, bz2, zst, lz4), uncompressed, or a mixture.",
          patterns: ["*.fasta", "*.fa", "*.fasta.gz", "*.fasta.zst"]
        }
        sample_metadata_tsvs: {
            description: "Tab-separated metadata file that contain binning variables and values. Must contain all samples: output will be filtered to the IDs present in this file.",
            patterns: ["*.txt", "*.tsv"]
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
        min_unambig_genome: {
          description: "Minimum number of called bases in genome to pass prefilter."
        }
        ancestral_traits_to_infer: {
          description: "A list of metadata traits to use for ancestral node inference (see https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html). Multiple traits may be specified; must correspond exactly to column headers in metadata file. Omitting these values will skip ancestral trait inference, and ancestral nodes will not have estimated values for metadata."
        }
        clades_tsv: {
          description: "A TSV file containing clade mutation positions in four columns: [clade  gene    site    alt]; see: https://nextstrain.org/docs/tutorials/defining-clades",
          patterns: ["*.tsv", "*.txt"]
        }
    }

    call augur_noSub.sarscov2_nextstrain_noSub {
    input:
        assembly_fastas=assembly_fastas,
        sample_metadata_tsvs=sample_metadata_tsvs,
        build_name=build_name
    }

    call phylo.snp_dists {
      input:
        cluster_name = build_name,
        alignment = sarscov2_nextstrain_noSub.mafft_alignment
    }
    call vis.cluster_render {
      input:
        cluster_name = build_name,
        snp_matrix = snp_dists.snp_matrix,
        ml_tree = sarscov2_nextstrain_noSub.ml_tree,
        render_template = render_template
    }

    output {
      File  combined_assemblies   = sarscov2_nextstrain_noSub.combined_assemblies
      File  MAFFT_alignment    = sarscov2_nextstrain_noSub.mafft_alignment
      File  unmasked_snps         = sarscov2_nextstrain_noSub.unmasked_snps

      File  metadata_merged       = sarscov2_nextstrain_noSub.metadata_merged
      File  keep_list             = sarscov2_nextstrain_noSub.keep_list
      File  time_tree             = sarscov2_nextstrain_noSub.time_tree

      File  auspice_input_json    = sarscov2_nextstrain_noSub.auspice_input_json
      File      analysis_doc = cluster_render.analysis_doc
      File      snp_list     = cluster_render.snp_list
      File      snp_matrix   = snp_dists.snp_matrix
    }
}
