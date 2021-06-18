version 1.0

import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_versioning.wdl" as versioning

workflow pangolin_update {
    input {
        String  samplename
        File    assembly
        String  updated_pangolin_docker
    }

    call taxon_ID.pangolin3 {
    input:
        samplename = samplename,
        fasta = assembly,
        docker = updated_pangolin_docker
    }
    call versioning.version_capture{
      input:
    }
    output {
        String pangolin_update_version            = version_capture.phvg_version
        String pangolin_update_analysis_date      = version_capture.date
        String pango_lineage        = pangolin3.pangolin_lineage
        String pangolin_conflicts   = pangolin3.pangolin_conflicts
        String pangolin_notes       = pangolin3.pangolin_notes
        String pangolin_version     = pangolin3.version
        File   pango_lineage_report = pangolin3.pango_lineage_report
        String pangolin_docker      = pangolin3.pangolin_docker
    }
}
