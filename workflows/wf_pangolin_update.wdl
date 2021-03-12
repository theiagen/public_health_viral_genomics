version 1.0

import "../tasks/task_taxonID.wdl" as taxon_ID

workflow pangolin_update {
	input {
		String    samplename
		String 		assembly
		String    current_pangolin_docker
		String    updated_pangolin_docker
	}

	call check_version {
		input:
			current_pangolin_docker = current_pangolin_docker,
			updated_pangolin_docker = updated_pangolin_docker
	}

	if (check_version.update_pangolin) {
	call taxon_ID.pangolin2 {
    input:
      samplename = samplename,
      fasta = assembly,
			docker = updated_pangolin_docker
	}
	}
	output {
			String?  pangolin_lineage       = pangolin2.pangolin_lineage
			Float?   pangolin_aLRT          = pangolin2.pangolin_aLRT
			File?    pango_lineage_report   = pangolin2.pango_lineage_report
			String?  pangolin_version       = pangolin2.version
			String? 	pangolin_docker 			 = pangolin2.pangolin_docker
	}
}

task check_version {
	input {
		String current_pangolin_docker
		String updated_pangolin_docker
	}
	command <<<
		if ~{current_pangolin_docker} != ~{updated_pangolin_docker}; then
			echo true > update_pangolin
		else
			echo false > update_pangolin
		fi
	>>>
	output {
		Boolean update_pangolin = read_boolean("update_pangolin")
	}
	runtime {
		docker:       "theiagen/utility:1.0"
		memory:       "1 GB"
    cpu:          1
    disks:        "local-disk 10 SSD"
    preemptible:  0
}
}
