version 1.0

import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_taxonID.wdl" as taxonID

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
	  call taxonID.kraken2 {
	    input:
	      samplename = samplename,
	      read1 = ncbi_scrub_pe.read1_dehosted,
				read2 = ncbi_scrub_pe.read2_dehosted
	  }

	output {
			File	read1_dehosted	=	ncbi_scrub_pe.read1_dehosted
			File 	read2_dehosted	=	ncbi_scrub_pe.read2_dehosted
			Int		read1_human_sports_removed = ncbi_scrub_pe.read1_human_spots_removed
			Int 	read2_human_sports_removed = ncbi_scrub_pe.read2_human_spots_removed
			String	ncbi_scrub_docker	=	ncbi_scrub_pe.ncbi_scrub_docker

			Float	kraken_human_dehosted	=	kraken2.percent_human
			Float	kraken_sc2_dehosted	=	kraken2.percent_sc2
			String	kraken_report_dehosted	=	kraken2.kraken_report
			String	kraken_version_dehosted	=	kraken2.version
	}
}
