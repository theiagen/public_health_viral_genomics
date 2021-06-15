version 1.0
import "../tasks/task_pub_repo_submission.wdl" as submission

workflow submission_files {
    input {
        String  samplename
        String  submission_id
        String  collection_date
        File    sequence
        File    read
    }

    call submission.deidentify {
        input:
            samplename = samplename,
            submission_id = submission_id,
            collection_date = collection_date,
            sequence = sequence,
            read = read
    }

    output {
        File   read_submission  = deidentify.read_submission
        File   deID_assembly    = deidentify.deID_assembly
        File?  genbank_assembly = deidentify.genbank_assembly
        File?  gisaid_assembly  = deidentify.gisaid_assembly
    }
}
