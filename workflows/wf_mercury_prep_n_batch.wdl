version 1.0

import "../tasks/task_mercury_file_wrangling.wdl" as submission


workflow mercury_prep_n_batch {
  input {

  }
  call submission.sm_metadata_wrangling {


  }
}