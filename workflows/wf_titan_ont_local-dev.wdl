version 1.0

import "wf_titan_clearlabs.wdl" as assembly

workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], File]] inputSamples
  }
  scatter (sample in inputSamples) {
    call assembly.titan_clearlabs {
      input:
        samplename = sample.left[0],
        clear_lab_fastq  = sample.right
    }
  }
}
