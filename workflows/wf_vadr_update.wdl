version 1.0

import "../tasks/task_ncbi.wdl" as ncbi

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

    output {
        File vadr_alerts_list = vadr.alerts_list
        Int vadr_num_alerts = vadr.num_alerts
        String vadr_docker = vadr.vadr_docker
    }
}
