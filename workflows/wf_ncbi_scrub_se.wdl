version 1.0

import "../tasks/task_read_clean.wdl" as read_clean

workflow dehost_se {
	input {
		String    samplename
		File 		reads
	}

	call read_clean.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = reads
  }

	output {
			File	reads_dehosted	=	ncbi_scrub_pe.read1_dehosted
			String	ncbi_scrub_docker	=	ncbi_scrub_pe.ncbi_scrub_docker
	}
}
