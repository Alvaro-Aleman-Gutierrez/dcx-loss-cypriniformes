#!/bin/bash
# Confirmar en lecturas crudas la ausencia de DCX en una segunda especie por linaje.
# Motivo (leccion de I. furcatus, 2026-11-05): un ensamblado cromosomico puede tener huecos.
# Danio: 8/8 solo por ensamblado (salvo D. rerio T2T). Paedocypris: solo micromegethes en lecturas.
set -uo pipefail
P=${PROJECT}; DD=${DATA}
LOG=$DD/reads_extra.log; OUT=$DD/reads_extra_results.txt
XSPOTS=${XSPOTS:-40000000}
: > $OUT; echo "[$(date)] inicio" > $LOG
for item in "reads_tinwini:ERR3284969:D_tinwini" "reads_carbunculus:SRR11971490:P_carbunculus"; do
  IFS=: read dir run tag <<< "$item"; d=$DD/$dir; mkdir -p $d; ln -sfn $d $P/$dir
  if [ ! -s $d/$run.fasta ]; then
    echo "[$(date)] descargando $run ($tag)" >> $LOG
    nice -n 15 fastq-dump --fasta 0 --split-spot --skip-technical -X $XSPOTS -O $d $run >> $LOG 2>&1 \
      || { echo "[$(date)] FALLO descarga $run" >> $LOG; continue; }
  fi
  echo "[$(date)] $tag descargado: $(du -h $d/$run.fasta | cut -f1)" >> $LOG
  if ! ls $d/reads_db*.nsq >/dev/null 2>&1; then
    nice -n 15 makeblastdb -in $d/$run.fasta -dbtype nucl -out $d/reads_db >> $LOG 2>&1 \
      || { echo "[$(date)] FALLO makeblastdb $run" >> $LOG; continue; }
  fi
  echo "[$(date)] prueba $tag" >> $LOG
  THREADS=3 nice -n 15 bash $P/scripts/reads_test.sh $d $tag >> $OUT 2>> $LOG && echo "" >> $OUT
done
echo "FIN reads_extra" >> $OUT; echo "[$(date)] COMPLETO" >> $LOG
