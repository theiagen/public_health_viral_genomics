version 1.0

import "../tasks/task_taxonID.wdl" as taxon_ID

workflow pangolin_update {
	input {
		String    samplename
		String 		assembly
		String    updated_pangolin_docker
	}

	call taxon_ID.pangolin2 {
    input:
      samplename = samplename,
      fasta = assembly,
			docker = updated_pangolin_docker
	}

	output {
			String  pangolin_lineage       = pangolin2.pangolin_lineage
			Float   pangolin_aLRT          = pangolin2.pangolin_aLRT
			File    pango_lineage_report   = pangolin2.pango_lineage_report
			String  pangolin_version       = pangolin2.version
			String	pangolin_docker 			 = pangolin2.pangolin_docker
	}
}
