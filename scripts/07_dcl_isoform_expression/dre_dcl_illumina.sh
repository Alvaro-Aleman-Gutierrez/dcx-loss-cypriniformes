#!/bin/bash
# Uniones de empalme de dclk1a en cerebro de pez cebra (SRR37646263, R1) contra una ventana genomica de 800 kb
# (NC_141030.1:37.9-38.7 Mb). minimap2 -x splice, solo primarias, PAF con CIGAR.
MM=~/bin_extra/minimap2-2.28_x64-linux/minimap2; O=${DATA}/rnaseq/dre_expr
F=${DATA}/rnaseq/pezcebra_illumina_cerebro/SRR37646263_1.fasta
cd $O; grep -c '^>' $F > n_R1.txt
$MM -x splice -c -t 2 --secondary=no dre_win.fa $F 2> mm.err > dre_R1.paf && echo "[$(date)] FIN $(wc -l < dre_R1.paf)" > dre.log
