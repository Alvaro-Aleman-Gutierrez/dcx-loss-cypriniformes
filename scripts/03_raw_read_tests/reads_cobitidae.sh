#!/bin/bash
# Prueba en lecturas crudas de la ausencia de DCX en Cobitidae (Misgurnus, Paramisgurnus).
# Descarga acotada (-X) a ~10 Gb por especie (~10x para genoma ~1 Gb), directo a D.
# Controles internos en la misma corrida: dclk1a, dclk2a, dclk3 (presentes en el ensamblado, p1b).
set -uo pipefail
P=${PROJECT}; DD=${DATA}
LOG=$DD/cobitidae_reads.log; OUT=$DD/cobitidae_reads_results.txt
XSPOTS=${XSPOTS:-35000000}
echo "[$(date)] inicio" >> $LOG
for item in "reads_misgurnus:SRR36982626:M_anguillicaudatus" "reads_paramisgurnus:SRR23104325:P_dabryanus"; do
  IFS=: read dir run tag <<< "$item"
  d=$DD/$dir; mkdir -p $d; ln -sfn $d $P/$dir
  if [ ! -s $d/$run.fasta ]; then
    echo "[$(date)] descargando $run ($tag)" >> $LOG
    nice -n 15 fastq-dump --fasta 0 --split-spot --skip-technical -X $XSPOTS -O $d $run >> $LOG 2>&1 \
      || { echo "[$(date)] FALLO descarga $run" >> $LOG; continue; }
  fi
  if ! ls $d/reads_db*.nsq >/dev/null 2>&1; then
    echo "[$(date)] makeblastdb $tag" >> $LOG
    nice -n 15 makeblastdb -in $d/$run.fasta -dbtype nucl -out $d/reads_db >> $LOG 2>&1 \
      || { echo "[$(date)] FALLO makeblastdb $run" >> $LOG; continue; }
  fi
  echo "[$(date)] prueba $tag" >> $LOG
  THREADS=3 nice -n 15 bash $P/scripts/reads_test.sh $d $tag >> $OUT 2>> $LOG && echo "" >> $OUT
  echo "[$(date)] fin $tag" >> $LOG
done
echo "FIN cobitidae" >> $OUT; echo "[$(date)] COMPLETO" >> $LOG
