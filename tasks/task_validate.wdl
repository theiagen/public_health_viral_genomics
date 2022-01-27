version 1.0

task export_two_tsvs {

  input {
    String        terra_project
    String        terra_workspace
    String        datatable1
    String        datatable2
    String        docker="broadinstitute/terra-tools:tqdm"
  }
  command {
    python3 scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project} --workspace ~{terra_workspace} --entity_type ~{datatable1} --tsv_filename ~{datatable1}

    python3 scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project} --workspace ~{terra_workspace} --entity_type ~{datatable2} --tsv_filename ~{datatable2}
  }
  runtime {
      docker: "~{docker}"
      memory: "4 GB"
      cpu:    2
      disks: "local-disk 50 HDD"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries:   3
  }
  output {
      File   datatable1_tsv     = "datatable1.tsv"
      File   datatable2_tsv     = "datatable2.tsv"
  }
}

task compare_two_tsv {
  input {
    File  tsv1_clean
    File  tsv2_clean
    String out_prefix
  }
  command{
    python3 compare_data_tables.py ~{tsv1_clean} ~{tsv2_clean} --out_file ~{out_prefix}
  }
  runtime {
      docker: Special_pdfkit_docker
      memory: "4 GB"
      cpu:    2
      disks: "local-disk 50 HDD"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries:   3
  }
  output {
      File  pdf_report     = ~{out_prefix}
  }
}

task drop_unnecessary_columns {
  input {
    File  tsv1
    String  tsv1_out
    File  tsv2
    String  tsv2_out
  }
  command {
    python3 <<CODE
    import csv
    import pandas as pd
    import codecs
    keep_list=['assembly_length_unambiguous', 'assembly_mean_coverage', 'kraken_human', 'kraken_sc2', 'nextclade_aa_dels', 'nextclade_clade','number_Degenerate', 'number_N', 'number_Total', 'pango_lineage', 'vadr_num_alerts']
    with codecs.open("~{tsv1}",'r') as tsv_file:
      tsv_reader1=csv.reader(tsv_file, delimiter="\t")
      df1=pd.read_csv(tsv_reader1, sep='\t')
      drop_list1 = []
      for i in df1.columns.values:
      	if i not in keep_list:
      		drop_list1.append(i)
      df1.drop(drop_list1, axis='columns', inplace=True)
      df1.to_csv(~{tsv1_out}, sep="\t", index=False)
    with codecs.open("~{tsv2}",'r') as tsv_file:
      tsv_reader1=csv.reader(tsv_file, delimiter="\t")
      df1=pd.read_csv(tsv_reader2, sep='\t')
      drop_list2 = []
      for i in df2.columns.values:
      	if i not in keep_list:
      		drop_list2.append(i)
      df2.drop(drop_list2, axis='columns', inplace=True)
      df2.to_csv(~{tsv2_out}, sep="\t", index=False)
    CODE
  }
  runtime {
      docker: "~{docker}"
      memory: "4 GB"
      cpu:    2
      disks: "local-disk 50 HDD"
      dx_instance_type: "mem1_ssd1_v2_x2"
      maxRetries:   3
  }
  output {
      File   datatable1_tsv_clean     = ~{tsv1_out}
      File   datatable2_tsv_clean     = ~{tsv2_out}
  }
}
