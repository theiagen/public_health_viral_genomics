version 1.0

import "../tasks/task_read_clean.wdl" as read_clean

workflow dehost_pe {
	input {
		String    samplename
		File 		read1
		File    read2
	}

	call read_clean.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2
  }

	output {
			File	read1_dehosted	=	ncbi_scrub_pe.read1_dehosted
			File 	read2_dehosted	=	ncbi_scrub_pe.read2_dehosted
			String	ncbi_scrub_docker	=	ncbi_scrub_pe.ncbi_scrub_docker
	}
}
