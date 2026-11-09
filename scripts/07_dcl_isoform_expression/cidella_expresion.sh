#!/bin/bash
# Expresion de DCX y de las isoformas de dclk1a (larga vs DCL) en 6 tejidos de C. idella (PRJNA1306682).
# Mapeo de R1 con minimap2 -x sr contra mRNA RefSeq; salida PAF (solo primarias). Un proceso a la vez, 2 hilos.
MM=~/bin_extra/minimap2-2.28_x64-linux/minimap2
D=${DATA}/rnaseq; O=$D/cid_expr
cd $O
for f in $D/cidella_tejidos/*/*_1.fasta; do
  t=$(basename $(dirname $f)); r=$(basename $f _1.fasta)
  [ -s $O/$t.$r.paf.ok ] && continue
  n=$(grep -c '^>' $f)
  $MM -x sr -t 2 --secondary=no $O/ref_cuant.fa $f 2>/dev/null > $O/$t.$r.paf
  echo $n > $O/$t.$r.paf.ok
  echo "[$(date)] OK $t $r lecturas=$n mapeadas=$(wc -l < $O/$t.$r.paf)" >> $O/cid_expr.log
done
echo "[$(date)] FIN" >> $O/cid_expr.log
