#!/bin/bash
# Punto 1b: verificar en el genoma las especies RefSeq sin DCX anotado (lochas Cobitidae + Ictalurus furcatus).
# Trabajo temporal en C (~/tmp_dcx), resultados y genomas finales en D.
set -uo pipefail
P=${PROJECT}; T=${HOME}/tmp_dcx/p1b; D=${DATA}/p1b
mkdir -p $T $D; cd $T
cat > list.txt <<EOF
GCF_027580225.2 Misgurnus_anguillicaudatus
GCF_030506205.2 Paramisgurnus_dabryanus
GCF_023375685.1 Ictalurus_furcatus
EOF
cp $P/genomes/dl.sh $T/dl.sh 2>/dev/null || cp ${DATA}/genomes/dl.sh $T/dl.sh
bash $T/dl.sh $T/list.txt
# consultas: DCX (Carassius, Astyanax) + vecinos sinténicos
cat $P/tbn/q_dcx.faa $P/syn/block_q.faa > q.faa
: > $D/p1b_results.txt
while read acc sp; do
  [ -s $T/$sp.fna ] || { echo "FALTA $sp" >> $D/p1b_results.txt; continue; }
  tblastn -query $P/tbn/q_dcx.faa -db $T/db_$sp -evalue 1e-5 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq" \
    -num_threads 2 -max_target_seqs 50 > $sp.tbn
  awk -v sp=$sp '$1 ~ /Caraur/ && $4>=40 {gsub("-","",$11); gsub("\\*","X",$11); print ">"sp"|"$2":"$7"-"$8"|qpos"$5"-"$6"|pid"$3"\n"$11}' $sp.tbn > $sp.hsp.faa
  blastp -query $sp.hsp.faa -db $P/tbn/famdb -outfmt "6 qseqid sseqid pident bitscore" -max_target_seqs 1 -max_hsps 1 2>/dev/null | sort -u -k1,1 > $sp.rbh
  { echo "== $sp ($acc)"; echo "-- HSPs de DCX (Carassius) y su mejor hit recíproco:";
    awk '{split($1,a,"|"); print "  ",a[2],a[3],a[4],"->",$2,$3"%"}' $sp.rbh
    echo "-- vecinos sinténicos (mejor locus):"
    tblastn -query $P/syn/block_q.faa -db $T/db_$sp -evalue 1e-10 -outfmt "6 qseqid sseqid pident length sstart bitscore" -num_threads 2 -max_target_seqs 3 -max_hsps 1 2>/dev/null \
      | sort -k1,1 -k6,6nr | awk '{c[$1]++} c[$1]<=1 {print "  ",$1,$2,$3"%",$5}'; } >> $D/p1b_results.txt
done < list.txt
# mover genomas y bases a D (verificando) y limpiar temporal en C
mkdir -p $D/genomes
for f in $T/*.fna $T/db_*; do cp "$f" $D/genomes/ && cmp -s "$f" "$D/genomes/$(basename $f)" && rm -f "$f"; done
cp $T/*.tbn $T/*.rbh $D/ 2>/dev/null
echo "FIN p1b $(date)" >> $D/p1b_results.txt
rmdir $T 2>/dev/null; ls $T 2>/dev/null | head
