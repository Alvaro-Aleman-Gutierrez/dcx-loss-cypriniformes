#!/bin/bash
# Verificar en lecturas crudas la ausencia de DCX en Ictalurus furcatus (Siluriformes).
# Importante: I. punctatus (especie hermana) SI tiene DCX -> hay que descartar hueco de ensamblado.
# Controles internos en la misma corrida: dclk1a, dclk2a, dclk3.
set -uo pipefail
P=${PROJECT}; DD=${DATA}
LOG=$DD/furcatus_reads.log; OUT=$DD/furcatus_reads_results.txt
XSPOTS=${XSPOTS:-40000000}
run=SRR18966194; tag=I_furcatus; d=$DD/reads_furcatus
mkdir -p $d; ln -sfn $d $P/reads_furcatus
echo "[$(date)] inicio $run" >> $LOG
if [ ! -s $d/$run.fasta ]; then
  nice -n 15 fastq-dump --fasta 0 --split-spot --skip-technical -X $XSPOTS -O $d $run >> $LOG 2>&1 \
    || { echo "[$(date)] FALLO descarga" >> $LOG; exit 1; }
fi
echo "[$(date)] descarga lista: $(du -h $d/$run.fasta | cut -f1)" >> $LOG
if ! ls $d/reads_db*.nsq >/dev/null 2>&1; then
  nice -n 15 makeblastdb -in $d/$run.fasta -dbtype nucl -out $d/reads_db >> $LOG 2>&1 \
    || { echo "[$(date)] FALLO makeblastdb" >> $LOG; exit 1; }
fi
echo "[$(date)] prueba $tag" >> $LOG
THREADS=4 nice -n 15 bash $P/scripts/reads_test.sh $d $tag >> $OUT 2>> $LOG
echo "FIN furcatus" >> $OUT; echo "[$(date)] COMPLETO" >> $LOG
