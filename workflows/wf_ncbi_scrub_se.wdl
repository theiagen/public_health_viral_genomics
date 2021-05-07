version 1.0

import "../tasks/task_read_clean.wdl" as read_clean
import "../tasks/task_taxonID.wdl" as taxonID

workflow dehost_se {
	input {
		String    samplename
		File 		reads
	}

	call read_clean.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = reads
  }
	call taxonID.kraken2 {
		input:
			samplename = samplename,
			read1 = ncbi_scrub_se.read1_dehosted
	}
	output {
			File	reads_dehosted	=	ncbi_scrub_se.read1_dehosted
			String	ncbi_scrub_docker	=	ncbi_scrub_se.ncbi_scrub_docker
			Int		human_sports_removed = ncbi_scrub_se.read1_human_spots_removed

			Float	kraken_human_dehosted	=	kraken2.percent_human
			Float	kraken_sc2_dehosted	=	kraken2.percent_sc2
			String	kraken_version_dehosted	=	kraken2.version
			String	kraken_report_dehosted	=	kraken2.kraken_report
	}
}
