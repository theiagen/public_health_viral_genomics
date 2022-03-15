version 1.0

task stats_n_coverage {
  input {
    File bamfile
    String samplename
    Int s_gene_start = 21563
    Int s_gene_stop = 25384
    Int? min_depth
  }
  command <<<
    date | tee DATE
    samtools --version | head -n1 | tee VERSION

    samtools stats ~{bamfile} > ~{samplename}.stats.txt

    samtools coverage ~{bamfile} -m -o ~{samplename}.cov.hist
    samtools coverage ~{bamfile} -o ~{samplename}.cov.txt
    samtools flagstat ~{bamfile} > ~{samplename}.flagstat.txt

    coverage=$(cut -f 6 ~{samplename}.cov.txt | tail -n 1)
    depth=$(cut -f 7 ~{samplename}.cov.txt | tail -n 1)
    meanbaseq=$(cut -f 8 ~{samplename}.cov.txt | tail -n 1)
    meanmapq=$(cut -f 9 ~{samplename}.cov.txt | tail -n 1)

    samtools index ~{bamfile}
    chr=$(samtools idxstats ~{bamfile} | cut -f 1 | head -1)
    samtools coverage -r "${chr}:~{s_gene_start}-~{s_gene_stop}" ~{bamfile} >> ~{samplename}.cov.txt
    s_gene_depth=$(cut -f 7 ~{samplename}.cov.txt | tail -n 1)

    # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
    #  above min_depth, then wc -l counts them all
    orf1ab=$(samtools depth -J -r "${chr}:21563-25384" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth} print;}' | wc -l )
    sgene=$(samtools depth -J -r "${chr}:~{s_gene_start}-~{s_gene_stop}" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    orf3a=$(samtools depth -J -r "${chr}:25393-26220" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    egene=$(samtools depth -J -r "${chr}:26245-26472" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    mgene=$(samtools depth -J -r "${chr}:26523-27191" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    orf6=$(samtools depth -J -r "${chr}:27202-27387" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    orf7a=$(samtools depth -J -r "${chr}:27394-27759" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    orf7b=$(samtools depth -J -r "${chr}:27756-27887" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
		orf8=$(samtools depth -J -r "${chr}:27894-28259" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ngene=$(samtools depth -J -r "${chr}:28274-29533" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
		orf10=$(samtools depth -J -r "${chr}:29558-29674" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )

    orf1ab_pc=$(python3 -c "print ( round( ($orf1ab / 21290 ) * 100, 2 ) )")
    sgene_pc=$(python3 -c "print ( round( ($sgene / 3822 ) * 100, 2 ) )")
    orf3a_pc=$(python3 -c "print ( round( ($orf3a / 828 ) * 100, 2 ) )")
    egene_pc=$(python3 -c "print ( round( ($egene / 228 ) * 100, 2 ) )")
    mgene_pc=$(python3 -c "print ( round( ($mgene / 669 ) * 100, 2 ) )")
    orf6_pc=$(python3 -c "print ( round( ($orf6 / 186 ) * 100, 2 ) )")
    orf7a_pc=$(python3 -c "print ( round( ($orf7a / 266 ) * 100, 2 ) )")
    orf7b_pc=$(python3 -c "print ( round( ($orf7b / 132 ) * 100, 2 ) )")
    orf8_pc=$(python3 -c "print ( round( ($orf8 / 366 ) * 100, 2 ) )")
    ngene_pc=$(python3 -c "print ( round( ($ngene / 1260 ) * 100, 2 ) )")
    orf10_pc=$(python3 -c "print ( round( ($orf6 / 117 ) * 100, 2 ) )")

    echo -e "Gene\tPercent_Coverage" > ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF1ab\t" $orf1ab_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "S_gene\t" $sgene_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF3a\t" $orf3a_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "E_gene\t" $egene_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "M_gene\t" $mgene_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF6\t" $orf6_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF7a\t" $orf7a_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF7b\t" $orf7b_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF8\t" $orf8_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "N_gene\t" $ngene_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ORF10\t" $orf10_pc >> ~{samplename}.percent_gene_coverage.tsv


    if [ -z "$coverage" ] ; then coverage="0" ; fi
    if [ -z "s_gene_depth" ] ; then s_gene_depth="0"; fi
    if [ -z "$depth" ] ; then depth="0" ; fi
    if [ -z "$meanbaseq" ] ; then meanbaseq="0" ; fi
    if [ -z "$meanmapq" ] ; then meanmapq="0" ; fi

    echo $coverage | tee COVERAGE
    echo $depth | tee DEPTH
    echo $s_gene_depth | tee S_GENE_DEPTH
    echo $sgene_pc | tee S_GENE_PC
    echo $meanbaseq | tee MEANBASEQ
    echo $meanmapq | tee MEANMAPQ
  >>>
  output {
    String date = read_string("DATE")
    String samtools_version = read_string("VERSION")
    File stats = "~{samplename}.stats.txt"
    File cov_hist = "~{samplename}.cov.hist"
    File cov_stats = "~{samplename}.cov.txt"
    File flagstat = "~{samplename}.flagstat.txt"
    Float coverage = read_string("COVERAGE")
    Float depth = read_string("DEPTH")
    Float s_gene_depth = read_string("S_GENE_DEPTH")
    Float s_gene_percent_coverage = read_string("S_GENE_PC")
    File percent_gene_coverage = "~{samplename}.percent_gene_coverage.tsv"
    Float meanbaseq = read_string("MEANBASEQ")
    Float meanmapq = read_string("MEANMAPQ")
  }
  runtime {
    docker: "quay.io/staphb/samtools:1.10"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
