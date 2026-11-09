#!/bin/bash
# Reanudar el trabajo tras suspender o apagar el equipo.
# PROBLEMA QUE RESUELVE: si un proceso muere a medias deja ficheros PARCIALES que los scripts
# tomarian por buenos (una base BLAST incompleta o un .tbn truncado dan resultados silenciosamente
# erroneos: menos lecturas, menos loci). Aqui se verifica la integridad y se borra lo dudoso.
set -uo pipefail
P=${PROJECT}; DD=${DATA}
echo "=== 1) bases de lecturas BLAST: verificar integridad"
for d in $DD/reads_*; do
  [ -d "$d" ] || continue
  ls $d/reads_db*.nsq >/dev/null 2>&1 || { echo "  $(basename $d): sin base (se creara)"; continue; }
  if blastdbcmd -db $d/reads_db -info >/dev/null 2>&1; then
    n=$(blastdbcmd -db $d/reads_db -info 2>/dev/null | grep -oE "^\s+[0-9,]+ sequences" | tr -d ' ,' | head -1)
    echo "  $(basename $d): OK ($n secuencias)"
  else
    echo "  $(basename $d): BASE CORRUPTA -> se borra para rehacerla"
    rm -f $d/reads_db.* 
  fi
done
echo "=== 2) caches tblastn del bloque: borrar los truncados"
B=$DD/syn2/bloque
if [ -d "$B" ]; then
  for f in $B/*.tbn.tsv; do
    [ -e "$f" ] || continue
    # un fichero valido termina en linea completa con 8 campos
    if [ -s "$f" ] && [ "$(tail -1 "$f" | awk -F'\t' '{print NF}')" = "8" ]; then
      echo "  $(basename $f): OK ($(wc -l < $f) lineas)"
    else
      echo "  $(basename $f): TRUNCADO -> se borra"; rm -f "$f"
    fi
  done
fi
echo "=== 3) relanzar lo que falte"
if [ ! -s $B/bloque_results.txt ] || ! grep -q "FIN bloque" $B/bloque_results.txt 2>/dev/null; then
  echo "  bloque sintenico: relanzando"
  THREADS=6 nice -n 10 python3 -I $P/scripts/bloque_dcx.py > $DD/bloque.log 2>&1 &
fi
if ! grep -q "FIN reads_extra" $DD/reads_extra_results.txt 2>/dev/null; then
  echo "  lecturas P. carbunculus: relanzando"
  nice -n 15 bash $P/scripts/reads_extra.sh >> $DD/reads_extra.log 2>&1 &
fi
wait
echo "=== LISTO. Revisa:"
echo "  $B/bloque_results.txt"
echo "  $DD/reads_extra_results.txt"
