version 1.0

task snp_dists {
  input {
    File alignment
    String cluster_name
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    snp-dists -v | tee VERSION

    snp-dists ~{alignment} > ~{cluster_name}_snp_distance_matrix.tsv
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File snp_matrix = "~{cluster_name}_snp_distance_matrix.tsv"
  }
  runtime {
    docker: "quay.io/staphb/snp-dists:0.6.2"
    memory: "2 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task reorder_matrix {
  input {
    File tree
    File matrix
    String cluster_name
    Int disk_size = 100
  }
  command <<<
    python3 <<CODE
    from Bio import Phylo
    import pandas as pd

    # read in newick tree
    tree = Phylo.read("~{tree}", "newick")
    
    # extract ordered terminal ends
    term_names = [term.name for term in tree.get_terminals()]

    # read in matrix into pandas data frame
    snps = pd.read_csv("~{matrix}", header=0, index_col=0, delimiter="\t")

    # ensure all header and index values are strings for proper reindexing
    # this is because if sample_name is entirely composed of integers, pandas 
    # auto-casts them as integers; get_terminals() interprets those as strings. 
    # this incompatibility leads to failure and an empty ordered SNP matrix
    snps.columns = snps.columns.astype(str)
    snps.index = snps.index.astype(str)

    # reorder matrix according to terminal ends
    snps = snps.reindex(index=term_names, columns=term_names)

    # add phandango suffix to ensure continuous coloring
    snps_out1 = snps.add_suffix(":c1")

    # write out reordered matrix to a file
    snps_out1.to_csv("~{cluster_name}_ordered_snp_distance_matrix.csv", sep=",")

    # reroot tree with midpoint
    tree.root_at_midpoint()

    # re-extract ordered terminal ends of rerooted tree
    term_names = [term.name for term in tree.get_terminals()]

    # reorder matrix with re-ordered terminal ends
    snps = snps.reindex(index=term_names, columns=term_names)

    # add phandango suffix to ensure continuous coloring
    snps_out2 = snps.add_suffix(":c1")

    # write out reordered matrix of rerooted tree to a file
    snps_out2.to_csv("~{cluster_name}_midpoint_snp_distance_matrix.csv", sep=",")

    # write rerooted tree to a file
    Phylo.write(tree, "~{cluster_name}_midpoint_tree.nwk", "newick")

    CODE
  >>>
  output{
    File ordered_matrix = "~{cluster_name}_ordered_snp_distance_matrix.csv"
    File ordered_midpoint_matrix = "~{cluster_name}_midpoint_snp_distance_matrix.csv"
    File midpoint_rooted_tree = "~{cluster_name}_midpoint_tree.nwk"
  }
  runtime {
    docker: "staphb/mykrobe:0.12.1" # used because it contains both biopython and pandas
    memory: "2 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
   # maxRetries: 3
    preemptible: 0
  }
}

task visualize_matrix {
  input {
    File matrix
    String cluster_name
    Int disk_size = 100
  }
  command <<<
    python3 <<CODE

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
    import numpy as np
    import pandas as pd

    # read in matrix into pandas data frame
    snps = pd.read_csv("~{matrix}", header=0, index_col=0)

    np_snps = snps.to_numpy()

    fig, ax = plt.subplots()
    im = ax.imshow(np_snps)

    # remove the ':c1' string from column labels
    snps.columns = [s.replace(':c1', '') for s in snps.columns]

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(snps.columns)), labels=snps.columns)
    ax.set_yticks(np.arange(len(snps.index)), labels=snps.index)

    # create an Axes on the right side of ax. padding between cax and ax will be fixed at 0.05 inch.
    # dynamically adjust cax
    size = len(snps.columns)*0.005
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=size, pad=0.05)

    # draw colorbar
    plt.colorbar(im, cax=cax)

    # Rotate the x-axis tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
      rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(snps.columns)):
      for j in range(len(snps.index)):
          text = ax.text(j, i, np_snps[i, j],
            ha="center", va="center", color="w")

    ax.set_title("~{cluster_name} SNP Matrix")
    # dynamically scale image size to number of samples
    fig.set_size_inches(len(snps.columns)/2,len(snps.index)/2) 
    # ensure all tick labels fit in chart area
    fig.tight_layout()
    plt.savefig("~{cluster_name}_matrix.png")

    CODE
  >>>
  output{
    File snp_matrix_plot = "~{cluster_name}_matrix.png"
  }
  runtime {
    docker: "staphb/freyja:1.3.9" 
    memory: "2 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
   # maxRetries: 3
    preemptible: 0
  }
}

task assign_clusters {
  input {
    File r_script = "gs://theiagen-public-files-rp/terra/sars-cov-2-files/assign_clusters.R"
    File matrix
    String cluster_name
    File 
    Int cluster_snp_threshold = "2"
    Int disk_size = 100
  }
  command <<<
  
  Rscript ~{r_script} ~(cluster_snp_threshold) ~{matrix} ~{merged_metadata} 
  
  >>>
  output{
    File cluster_report = "cluster_report.csv"
  }
  runtime {
    docker: "rocker/tidyverse:4.2" 
    memory: "2 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
   # maxRetries: 3
    preemptible: 0
  }
}

task iqtree {
  input {
    File alignment
    String cluster_name
    String? iqtree_model = "GTR+G4"
    String? iqtree_bootstraps = 1000
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    iqtree --version | grep version | sed 's/.*version/version/;s/ for Linux.*//' | tee VERSION

    numGenomes=`grep -o '>' ~{alignment} | wc -l`
    if [ $numGenomes -gt 3 ]
    then
      cp ~{alignment} msa.fasta
      iqtree -nt AUTO -s msa.fasta -m ~{iqtree_model} -bb ~{iqtree_bootstraps}
      cp msa.fasta.contree ~{cluster_name}_msa.tree
    fi
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File ml_tree = "~{cluster_name}_msa.tree"
  }
  runtime {
    docker: "quay.io/staphb/iqtree:1.6.7"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
