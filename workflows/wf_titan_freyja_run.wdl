version 1.0

import "../tasks/task_taxonID.wdl" as taxon_id
import "../tasks/task_versioning.wdl" as versioning

workflow titan_freyja_run {
	input {
		File primer_trimmed_bam
    String samplename
	}
	call taxon_id.freyja_variants_one_sample as variants {
		input:
			primer_trimmed_bam = primer_trimmed_bam,
      samplename = samplename
	}
  call taxon_id.freyja_demix_one_sample as demix {
    input:
      freyja_variants = variants.freyja_variants,
      freyja_depths = variants.freyja_depths,
      samplename = samplename
  }
	call versioning.version_capture{
    input:
  }
  output {
    String titan_freyja_run_wf_version = version_capture.phvg_version
    String titan_freyja_run_wf_analysis_date = version_capture.date
		
		File freyja_variants = variants.freyja_variants
		File freyja_depths = variants.freyja_depths
    File freyja_demixed = demix.freyja_demixed
		
		
    }
}