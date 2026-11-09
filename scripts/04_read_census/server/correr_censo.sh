#!/bin/bash
# Censo de DCX en lecturas WGS crudas de Cypriniformes (SRA). Una tarea por especie, N en paralelo.
# Por especie: descarga los primeros 15 M pares (FASTA) -> DIAMOND blastx contra la familia DCX/DCLK -> parse_dmnd.py
# -> borra las lecturas. Reanudable: las especies con resultados/<run>.done se saltan.
# uso: ./correr_censo.sh [N_paralelo=20] [hilos_por_especie=4]
set -uo pipefail
cd "$(dirname "$0")"; H=$PWD; P=${1:-20}; T=${2:-4}; NSPOTS=${NSPOTS:-15000000}
SRA=$(ls -d $H/bin/sratoolkit*/bin | head -1); mkdir -p resultados tmp
una() {
  run=$1; sp=$2; R=$H/resultados; TMP=$H/tmp/$run
  [ -e $R/$run.done ] && return 0
  mkdir -p $TMP
  for i in 1 2 3 4; do
    $SRA/fastq-dump --fasta 0 --split-spot --skip-technical -X $NSPOTS --stdout $run > $TMP/r.fasta 2> $R/$run.dump.err
    grep -q "^Written" $R/$run.dump.err && break
    sleep $((60 * i))
  done
  if ! grep -q "^Written" $R/$run.dump.err; then echo "FALLO_DESCARGA" > $R/$run.fallo; rm -rf $TMP; return 1; fi
  grep "^Written" $R/$run.dump.err > $R/$run.spots
  $H/bin/diamond blastx -q $TMP/r.fasta -d $H/fam.dmnd -e 1e-5 -k 5 --max-hsps 1 -p $T -b ${BLOCK:-1} -c 1 --quiet \
     -f 6 qseqid sseqid pident length bitscore -o $R/$run.dmnd.tsv --tmpdir $TMP \
     || { echo "FALLO_DIAMOND" > $R/$run.fallo; rm -rf $TMP; return 1; }
  python3 $H/parse_dmnd.py $R/$run.dmnd.tsv "$sp|$run" > $R/$run.res && touch $R/$run.done
  gzip -f $R/$run.dmnd.tsv; rm -rf $TMP $R/$run.fallo
}
export -f una; export H SRA T NSPOTS
date > resultados/inicio.txt
LISTA=${LISTA:-especies.tsv}
awk -F'\t' '{gsub(/[ \x27]/,"_",$2); print $1, $2}' $LISTA | xargs -P $P -n 2 bash -c 'una "$0" "$1"'
cat resultados/*.res > censo_DCX_resultados.tsv
echo "Especies con resultado: $(ls resultados/*.done | wc -l) de $(wc -l < especies.tsv); fallidas: $(ls resultados/*.fallo 2>/dev/null | wc -l)"
date > resultados/fin.txt
