#!/bin/bash
# Aclarar el HSP ambiguo de I. furcatus: NC_071271.1:24530778-24531116 con mejor hit DCX de pollo.
# Hipotesis: es un exon del dominio DC de dclk1b, no un resto de DCX.
set -uo pipefail
P=${PROJECT}; O=${DATA}/furcatus; B=${DATA}/p1b/genomes
mkdir -p $O; cd $O
echo "NC_071271.1:24515000-24545000" > r.list
seqkit faidx -l r.list $B/Ictalurus_furcatus.fna > region.fna
# que genes de la familia mapean en esta ventana de 30 kb
tblastn -query $P/census/queries.faa -subject region.fna -evalue 1e-3 -seg no -max_hsps 20 \
  -outfmt "6 qseqid pident length qstart qend sstart send evalue bitscore" > fam_en_region.tsv
{ echo "== Genes de la familia con HSPs en NC_071271.1:24515000-24545000 (I. furcatus)";
  sort -k9,9nr fam_en_region.tsv | awk '!seen[$1]++ {printf "  %-12s mejor HSP: q%s-%s len=%s pid=%s%% bit=%s\n",$1,$4,$5,$3,$2,$9}'
  echo; echo "== Orden de los HSPs a lo largo de la region (gen : posicion):";
  sort -k6,6n fam_en_region.tsv | awk '$3>=30 {printf "  %9d  %-12s q%s-%s len=%s pid=%s%%\n",$6,$1,$4,$5,$3,$2}'
} > furcatus_results.txt
cat furcatus_results.txt
