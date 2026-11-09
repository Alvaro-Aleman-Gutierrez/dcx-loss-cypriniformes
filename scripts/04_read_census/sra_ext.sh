#!/bin/bash
# Ampliacion del muestreo de DCX con lecturas WGS de SRA (sobre todo PRJNA1483025, NMBE; y PRJEB19273).
# Dos procesos: descarga (15 M pares por especie, FASTA) y analisis DIAMOND (reads_diamond.sh); se borran las lecturas.
set -uo pipefail
P=${PROJECT}; O=${DATA}/sra_ext; cd $O
LISTA=$O/lista.tsv; RES=$O/sra_ext_results.tsv; N=${NSPOTS:-15000000}
descargar() {
  while IFS=$'\t' read run sp; do
    { [ -e $O/$run.done ] || [ -e $O/$run.fasta.ok ] || [ -e $O/$run.fasta ]; } && continue
    for i in 1 2 3; do
      fastq-dump --fasta 0 --split-spot --skip-technical -X $N --stdout $run > $O/$run.fasta 2> $O/$run.dump.err
      cat $O/$run.dump.err >> $O/descarga.log
      grep -q "^Written $N spots" $O/$run.dump.err && break      # solo se acepta la descarga completa
      rm -f $O/$run.fasta
      echo "[$(date)] reintento $i $run" >> $O/descarga.log; sleep 30
    done
    grep -q "^Written $N spots" $O/$run.dump.err 2>/dev/null && touch $O/$run.fasta.ok || echo "[$(date)] FALLO $run" >> $O/descarga.log
    # no adelantarse mas de 2 especies al analisis (disco)
    while [ $(ls $O/*.fasta.ok 2>/dev/null | wc -l) -ge ${GATE:-3} ]; do sleep 30; done
  done < $LISTA
  touch $O/descarga.fin
}
analizar() {
  while true; do
    ok=$(ls $O/*.fasta.ok 2>/dev/null | head -1)
    if [ -z "$ok" ]; then [ -e $O/descarga.fin ] && break; sleep 30; continue; fi
    run=$(basename $ok .fasta.ok); sp=$(awk -F'\t' -v r=$run '$1==r{print $2}' $LISTA | tr ' ' '_')
    THREADS=${THREADS:-4} bash $P/scripts/reads_diamond.sh $O/$run.fasta "$sp|$run" $O/$run.dmnd.tsv >> $RES 2>> $O/analisis.log \
      && { rm -f $O/$run.fasta; touch $O/$run.done; } || echo "[$(date)] FALLO analisis $run" >> $O/analisis.log
    rm -f $ok
  done
  echo "FIN" >> $RES
}
descargar & analizar; wait
