#!/bin/bash
# Espera a que cada base de lecturas esté lista y ejecuta reads_test.sh; resultados -> resultados_lecturas.txt
P=${PROJECT}
OUT=$P/resultados_lecturas.txt
LOG=$P/auto_pipeline.log
echo "[$(date)] inicio" >> $LOG
ready() { # base lista = existe y no hay makeblastdb/fastq-dump activos sobre ese directorio
  local d=$1
  ls $d/reads_db*.nsq >/dev/null 2>&1 || return 1
  pgrep -f "fastq-dump.*$2" >/dev/null && return 1
  pgrep -af "makeblastdb -in $2" >/dev/null && return 1
  return 0
}
for item in "reads_cerebrum_ctrl:ERR3284980:D_cerebrum_CONTROL_POSITIVO" "reads_dracula:ERR2639754:D_dracula" "reads_paedocypris:SRR11971489:P_micromegethes"; do
  IFS=: read dir run tag <<< "$item"; d=$P/$dir
  n=0
  until ready $d $run; do sleep 120; n=$((n+1)); [ $n -gt 540 ] && { echo "[$(date)] TIMEOUT $tag" >> $LOG; continue 2; }; done
  echo "[$(date)] prueba $tag" >> $LOG
  bash $P/scripts/reads_test.sh $d $tag >> $OUT 2>> $LOG && echo "" >> $OUT
  echo "[$(date)] fin $tag" >> $LOG
done
{ echo; echo "## Prueba en lecturas crudas ($(date +%F))"; cat $OUT; } >> $P/piloto_resumen.md
echo "[$(date)] PIPELINE COMPLETO" >> $LOG
