version 1.0

import "../tasks/task_taxonID.wdl" as taxon_ID

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

	output {
		String	pango_lineage	=	pangolin3.pangolin_lineage
    String	pangolin_conflicts	=	pangolin3.pangolin_conflicts
    String pangolin_notes = pangolin3.pangolin_notes
  	String	pangolin_version	=	pangolin3.version
  	File	pango_lineage_report	=	pangolin3.pango_lineage_report
  	String	pangolin_docker	=	pangolin3.pangolin_docker
		String	pangolin_usher_version = pangolin3.pangolin_usher_version
	}
}
