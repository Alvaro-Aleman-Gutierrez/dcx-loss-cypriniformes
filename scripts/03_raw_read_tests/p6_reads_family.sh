#!/bin/bash
# Punto 6: verificar en lecturas crudas las ausencias de dcdc2 (Paedocypris) y dcdc2c (Danionella).
# Controles cruzados: dcdc2c presente en Paedocypris; dcdc2 presente en Danionella; dcdc2b y dclk3 presentes en todos.
# Clasificación por BLASTx contra las 16 proteínas de zebrafish/consulta (etiquetas limpias).
set -euo pipefail
P=${PROJECT}; O=${DATA}/p6; mkdir -p $O; cd $O
T=${THREADS:-2}
python3 -I - <<EOF
from Bio import SeqIO
keep={"Q_dcdc2","Q_dcdc2c","Q_dcdc2b","Q_dclk3"}
SeqIO.write([r for r in SeqIO.parse("$P/census/queries.faa","fasta") if r.id in keep],"q.faa","fasta")
EOF
makeblastdb -in $P/census/queries.faa -dbtype prot -out qdb >/dev/null
: > p6_results.txt
for d in reads_paedocypris reads_cerebrum_ctrl reads_dracula; do
  RD=$P/$d; mkdir -p $d
  [ -s $d/tbn.tsv ] || tblastn -query q.faa -db $RD/reads_db -evalue 1e-5 -outfmt "6 qseqid sseqid pident length" \
      -max_target_seqs 100000 -num_threads $T -seg no > $d/tbn.tsv
  cut -f2 $d/tbn.tsv | sort -u > $d/ids.txt
  seqkit grep -j $T -f $d/ids.txt $(ls $RD/*.fasta) | seqkit rename > $d/reads.fna
  blastx -query $d/reads.fna -db qdb -evalue 1e-5 -outfmt "6 qseqid sseqid pident length bitscore" -max_target_seqs 3 -max_hsps 1 -num_threads $T > $d/bx.tsv
  python3 -I - $d $d/bx.tsv >> p6_results.txt <<'PY'
import sys,collections
tag,fn=sys.argv[1],sys.argv[2]; best={}
for l in open(fn):
    q,s,p,ln,b=l.split(); b=float(b)
    if q not in best or b>best[q][0]: best[q]=(b,s,float(p))
L={"Q_dcdc2":420,"Q_dcdc2c":385,"Q_dcdc2b":454,"Q_dclk3":721}
c=collections.Counter(v[1] for v in best.values()); pid=collections.defaultdict(list)
for v in best.values(): pid[v[1]].append(v[2])
print(f"## {tag}")
for g in L: n=c.get(g,0); print(f"{g[2:]:8s} reads={n:5d} reads/100aa={100*n/L[g]:6.2f} pid={sum(pid[g])/max(1,len(pid[g])):5.1f}")
PY
done
# resumen de modelos con posible pseudogenización (ψ) en el censo
python3 -I - > psi_models.tsv <<'EOF'
import csv
rows=list(csv.reader(open("${PROJECT}/census/census_table.tsv"),delimiter="\t"))
h=rows[0]; print("especie\tgen\testado(cobertura;identidad)")
for r in rows[1:]:
    for g,v in zip(h[1:],r[1:]):
        if v.startswith(("pseudo?","parcial")): print(r[0],g,v,sep="\t")
EOF
echo "FIN punto 6" >> p6_results.txt
