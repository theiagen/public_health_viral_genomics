# Public Health Viral Genomics
Bioinformatics workflows for genomic characterization and epidemiology of viral pathogens of concern, especially the SARS-CoV-2 virus, and for preparation for public database submission. These workflows are written in [WDL](https://github.com/openwdl/wdl), a language for specifying data processing workflows with a human-readable and writeable syntax. They have been developed by [Theiagen genomics](https://theiagen.com/) to primarily run on the [Terra.bio](https://terra.bio/) platform but can be used with a command-line interface running WDL. More information about the steps undertaken in these workflows is available via the [Public Health Viral Genomics ReadTheDocs page](https://public-health-viral-genomics-theiagen.readthedocs.io/en/latest/overview.html). Support for running these workflows can be sought by raising a [GitHub issue](https://github.com/theiagen/public_health_viral_genomics/issues/new) or by contacting Theiagen [directly](https://theiagen.com/contact/).


### Contributors & Influence
* Based on collaborative work with Andrew Lang, PhD & his [Genomic Analysis WDL workflows](https://github.com/AndrewLangvt/genomic_analyses)
* Workflows and task development influenced by The Broad's [Viral Pipes](https://github.com/broadinstitute/viral-pipelines)
* TheiaCoV workflows for genomic characterization influenced by UPHL's [Cecret](https://github.com/UPHL-BioNGS/Cecret) & StaPH-B's [Monroe](https://staph-b.github.io/staphb_toolkit/workflow_docs/monroe/)
* The TheiaCoV workflow for waste water variant calling (TheiaCoV_WWVC) incorporates a modified version of the [CDPHE's WasteWaterVariantCalling WDL Worfklow](https://github.com/CDPHE/WasteWaterVariantCalling).
* The TheiaCoV workflow user community. To provide feedback, please raise a GitHub issue [GitHub issue](https://github.com/theiagen/public_health_viral_genomics/issues/new).

###Contributing to the PHVG workflows
Contributions to the workflows contained in this repository are warmly welcomed. Our style guide may be found below for convenience of formatting.

 
2-space indents (no tabs), braces on same line, single space when defining input/output variables & runtime attributes, single-line breaks between non-intended constructs, and task commands enclosed with triple braces (`<<< ... >>>`). 

<em>E.g.</em>:
```
workflow w {
  input {
    String input
  }
  call task_01 {
    input:
      input = input
  }
  call task_02 {
    input: 
      input = input
  }
  output {
    File task_01_out = task_01.output
    File task_02_out = task_02.output 
  }      
}

task task1 {
  input {
    String input
    String docker = "theiagen/utility:1.1"
  }
  command <<<
    echo '~{input}' > output.txt
  >>>
  output {
    File output = "output.txt"
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 0
  }
}

task task_02 {
  input{
    String input
    String docker = "theiagen/utility:1.1"
  }
  command <<<
    echo '~{input}' > output.txt
  >>>
  output {
    File output = "output.txt"
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 0
  }
}
```
Style guide inspired by [scottfrazer](https://gist.github.com/scottfrazer)'s [WDL Best Pratcices Style Guide](https://gist.github.com/scottfrazer/aa4ab1945a6a4c331211)
