#!/bin/bash
# Prueba de presencia/ausencia de DCX en lecturas WGS crudas, version rapida con DIAMOND blastx.
# Equivalente a reads_test.sh (misma base de la familia completa, misma asignacion por mejor hit, mismos controles
# monocopia dclk1a/dclk2a/dclk3), pero sin tBLASTn previo ni base de datos de lecturas.
# uso: reads_diamond.sh <reads.fasta> <etiqueta> <salida.tsv>
set -euo pipefail
F=$1; TAG=$2; OUT=$3; T=${THREADS:-4}; DB=${DATA}/sra_ext/fam.dmnd
[ -s $OUT ] || ~/bin_extra/diamond blastx -q $F -d $DB -e 1e-5 -k 5 --max-hsps 1 -p $T -b 0.5 -c 1 --quiet \
  -f 6 qseqid sseqid pident length bitscore -o $OUT
python3 ${PROJECT}/scripts/parse_dmnd.py "$OUT" "$TAG"
