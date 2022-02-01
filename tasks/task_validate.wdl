version 1.0

task export_two_tsvs {

  input {
    String        terra_project
    String        terra_workspace
    String        datatable1
    String        datatable2
  }
  command {
    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project} --workspace ~{terra_workspace} --entity_type ~{datatable1} --tsv_filename ~{datatable1}

    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project} --workspace ~{terra_workspace} --entity_type ~{datatable2} --tsv_filename ~{datatable2}
  }
  runtime {
      docker: "broadinstitute/terra-tools:tqdm"
      memory: "1 GB"
      cpu: 1
      disks: "local-disk 10 HDD"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries: 3
  }
  output {
      File   datatable1_tsv     = "~{datatable1}.tsv"
      File   datatable2_tsv     = "~{datatable2}.tsv"
  }
}

task compare_two_tsvs {
  input {
    File  datatable1_tsv
    File  datatable2_tsv
    String out_dir
    String out_prefix
  }
  command{
    python3 compare_data_tables.py ~{datatable1_tsv} ~{datatable2_tsv} --outdir ~{out_dir} --prefix ~{out_prefix}
  }
  runtime {
      docker: "quay.io/theiagen/utility:1.2"
      memory: "4 GB"
      cpu:    2
      disks: "local-disk 50 HDD"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries:   3
  }
  output {
      File  pdf_report = "~{out_dir}/~{out_prefix}.pdf"
      File xl_report = "~{out_dir}/~{out_prefix}.xlsx"
  }
}
