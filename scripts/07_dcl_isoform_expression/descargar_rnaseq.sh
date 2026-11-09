#!/bin/bash
# Descargas de RNA-seq para Q4 (expresion de la isoforma DCL de dclk1a y de DCX). Secuencial, con reintentos y
# verificacion ("Written" en el registro de fastq-dump). Reanudable: se salta lo que ya tiene <run>.ok.
#  1) Pez cebra (sin DCX), Iso-Seq de cerebro: SRR32588721 (PRJNA1232602, Sequel II), completo.
#  2) Pez cebra, RNA-seq Illumina de cerebro completo, silvestre: SRR37646263 (PRJNA1438185, WT-1), 20 M pares.
#  3) Ctenopharyngodon idella (con DCX), 6 tejidos x 2 individuos (PRJNA1306682, grupo 'tolerance'), 10 M pares c/u.
# Ya descargado: Paedocypris micromegethes Iso-Seq (rnaseq_dcl/Pmic_SRR11971539.fq, 19,7 GB).
set -uo pipefail
O=${DATA}/rnaseq; mkdir -p $O; cd $O
bajar() {   # bajar <dir> <run> <spots|full> <single|paired>
  d=$1; run=$2; n=$3; tipo=$4; mkdir -p $O/$d
  [ -e $O/$d/$run.ok ] && { echo "[$(date)] ya estaba $d/$run"; return 0; }
  X=""; [ "$n" != "full" ] && X="-X $n"
  S=""; [ "$tipo" = "paired" ] && S="--split-files"
  for i in 1 2 3 4; do
    nice -n 15 fastq-dump --fasta 0 --skip-technical $S $X -O $O/$d $run > $O/$d/$run.log 2>&1
    if grep -q "^Written" $O/$d/$run.log; then touch $O/$d/$run.ok; echo "[$(date)] OK $d/$run $(du -sh $O/$d/${run}* | head -1 | cut -f1)"; return 0; fi
    echo "[$(date)] reintento $i $d/$run"; rm -f $O/$d/${run}*.fasta; sleep $((60 * i))
  done
  echo "[$(date)] FALLO $d/$run"
}
bajar pezcebra_isoseq_cerebro SRR32588721 full single
bajar pezcebra_illumina_cerebro SRR37646263 20000000 paired
while read run tejido; do bajar cidella_tejidos/$tejido $run 10000000 paired; done <<'L'
SRR35031399 cerebro
SRR35031400 cerebro
SRR35029734 branquias
SRR35029735 branquias
SRR35030603 intestino
SRR35030604 intestino
SRR35032171 rinon
SRR35032172 rinon
SRR35030834 higado
SRR35030835 higado
SRR35031011 bazo
SRR35031012 bazo
L
echo "[$(date)] FIN descargas RNA-seq"
