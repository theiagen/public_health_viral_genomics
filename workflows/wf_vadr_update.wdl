version 1.0

import "../tasks/task_ncbi.wdl" as ncbi
import "../tasks/task_versioning.wdl" as versioning

workflow vadr_update {
    input {
        File genome_fasta
        String docker
    }
    call ncbi.vadr {
    input:
      genome_fasta = genome_fasta,
            docker = docker
    }
    call versioning.version_capture{
      input:
    }
    output {
        String vaddr_update_version            = version_capture.phvg_version
        String vadr_update_analysis_date      = version_capture.date
        File? vadr_alerts_list = vadr.alerts_list
        String vadr_num_alerts = vadr.num_alerts
        String vadr_docker = vadr.vadr_docker
    }
}
