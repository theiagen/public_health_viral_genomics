version 1.0

import "../tasks/task_taxonID.wdl" as taxon_ID
import "../tasks/task_versioning.wdl" as versioning

workflow pangolin_update {
    input {
        String  samplename
        File    assembly
        String current_lineage
        String current_pangolin_docker
        String current_pangolin_assignment_version
        String current_pangolin_usher_version
        String current_pangolin_versions
        String  updated_pangolin_docker
        String? timezone
        File? lineage_log
    }

    call taxon_ID.pangolin3 {
      input:
        samplename = samplename,
        fasta = assembly,
        docker = updated_pangolin_docker
    }
    call taxon_ID.pangolin_update_log {
      input:
        samplename = samplename,
        current_lineage = current_lineage,
        current_pangolin_docker = current_pangolin_docker,
        current_pangolin_assignment_version = current_pangolin_assignment_version,
        current_pangolin_usher_version = current_pangolin_usher_version,
        current_pangolin_versions = current_pangolin_versions,
        updated_lineage = pangolin3.pangolin_lineage,
        updated_pangolin_docker = pangolin3.pangolin_docker,
        updated_pangolin_assignment_version = pangolin3.pangolin_assignment_version,
        updated_pangolin_usher_version = pangolin3.pangolin_usher_version,
        updated_pangolin_versions = pangolin3.pangolin_versions,
        timezone = timezone,
        lineage_log = lineage_log
    }
    call versioning.version_capture{
      input:
        timezone = timezone
    }
    output {
        String pangolin_update_version            = version_capture.phvg_version
        String pangolin_update_analysis_date      = version_capture.date
        
        String pango_lineage        = pangolin3.pangolin_lineage
        String pangolin_conflicts   = pangolin3.pangolin_conflicts
        String pangolin_notes       = pangolin3.pangolin_notes
        String pangolin_assignment_version = pangolin3.pangolin_assignment_version
        String pangolin_versions = pangolin3.pangolin_versions
        File   pango_lineage_report = pangolin3.pango_lineage_report
        String pangolin_docker      = pangolin3.pangolin_docker
        
        String pangolin_updates = pangolin_update_log.pangolin_updates
        File pango_lineage_log = pangolin_update_log.pango_lineage_log
    }
}
