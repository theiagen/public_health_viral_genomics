import "../tasks/task_pub_repo_submission.wdl" as submission

workflow submission_files {
    input {
        Array[String] inputdata
        Array[File]   inputfiles
    }

    call submission.deidentify {
        input:
            samplename = inputdata[0],
            submission_id = inputdata[1],
            collection_date = inputdata[2],
            sequence = inputfiles[0],
            read1 = inputfiles[1],
            read2 = inputfiles[2]
    }
    output {
        File   read1_submission = deidentify.read1_submission
        File   read2_submission = deidentify.read2_submission
        File   deID_assembly    = deidentify.deID_assembly
        File?  genbank_assembly = deidentify.genbank_assembly
        File?  gisaid_assembly  = deidentify.gisaid_assembly
    }
}


