version 1.0

import "../tasks/task_test.wdl" as test

workflow test {

  input {
    File          text
  }

  call test.vadr {
    input:
      text = text
  }

  output {
    File?      test_fail = vadr.vadr_failed
    File?      test_pass     = vadr.vadr_passed
    Int num_alerts = vadr.num_alerts
  }
}
