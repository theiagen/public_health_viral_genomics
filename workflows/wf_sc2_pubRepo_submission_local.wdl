version 1.0 

import "wf_sc2_pubRepo_submission.wdl" as submission

workflow SC2_submission_files_local {
    input {
        Array[String]  inputdata
        Array[File]    inputfiles
    }

    call submission.SC2_submission_files {
        input:
            samplename      = inputdata[0],
            submission_id   = inputdata[1],
            collection_date = inputdata[2],
            coverage        = inputdata[3],
            number_N        = inputdata[4],
            number_ATCG     = inputdata[5],
            number_Total    = inputdata[6],
            sequence        = inputfiles[0],
            read1           = inputfiles[1],
            read2           = inputfiles[2]
    }
    output {
        File      deID_assembly    = SC2_submission_files.deID_assembly
        File      read1_submission = SC2_submission_files.read1_submission
        File      read2_submission = SC2_submission_files.read2_submission
        File?     genbank_assembly = SC2_submission_files.genbank_assembly
        File?     gisaid_assembly  = SC2_submission_files.gisaid_assembly
    }
}


