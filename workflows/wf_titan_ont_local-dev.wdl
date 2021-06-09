version 1.0

import "wf_titan_ont.wdl" as assembly

workflow nCoV19_pipeline {
  input {
    Array[Pair[Array[String], File]] inputSamples
  }
  scatter (sample in inputSamples) {
    call assembly.titan_ont {
      input:
        samplename = sample.left[0],
        demultiplexed_reads  = sample.right
    }
  }
}
