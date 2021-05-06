version 1.0

import "wf_titan_illumina_se.wdl" as assembly
import "../tasks/task_amplicon_metrics.wdl" as assembly_metrics
import "../tasks/task_sample_metrics.wdl" as summary


workflow nCoV19_pipeline {
input {
  Array[Pair[Array[String], File]] inputSamples
  File primer_bed
}
  scatter (sample in inputSamples) {
    call assembly.titan_illumina_se {
      input:
        samplename = sample.left[0],
        primer_bed = primer_bed,
        read1_raw  = sample.right
    }
  }
}
