#!/bin/bash
# Candidato del censo v2: dcdc2c ausente en Ictalurus furcatus. Sospechoso, porque ese mismo ensamblado
# daba DCX como ausente y las lecturas demostraron que si lo tiene. Reutiliza la base de lecturas ya creada.
# Controles internos: dcdc2, dcdc2b, dclk3 (presentes en el censo v2 para I. furcatus).
set -uo pipefail
P=${PROJECT}; DD=${DATA}; O=$DD/furcatus_dcdc2c
RD=$DD/reads_furcatus; T=${THREADS:-4}; mkdir -p $O; cd $O
ls $RD/reads_db*.nsq >/dev/null 2>&1 || { echo "falta la base de lecturas de I. furcatus"; exit 1; }
python3 -I - <<PY
from Bio import SeqIO
keep={"Q_dcdc2c","Q_dcdc2","Q_dcdc2b","Q_dclk3"}
SeqIO.write([r for r in SeqIO.parse("$P/census/queries.faa","fasta") if r.id in keep],"q.faa","fasta")
PY
# base discriminante: toda la familia (16) -> un fragmento de dcdc2 no puede contarse como dcdc2c
makeblastdb -in $P/census/queries.faa -dbtype prot -out fam >/dev/null
[ -s tbn.tsv ] || nice -n 15 tblastn -query q.faa -db $RD/reads_db -evalue 1e-5 -seg no \
   -outfmt "6 qseqid sseqid pident length" -max_target_seqs 200000 -num_threads $T > tbn.tsv
cut -f2 tbn.tsv | sort -u > ids.txt
[ -s reads.fna ] || seqkit grep -j $T -f ids.txt $(ls $RD/*.fasta) | seqkit rename > reads.fna
[ -s bx.tsv ] || nice -n 15 blastx -query reads.fna -db fam -evalue 1e-5 \
   -outfmt "6 qseqid sseqid pident length bitscore" -max_target_seqs 5 -max_hsps 1 -num_threads $T > bx.tsv
python3 -I - > dcdc2c_results.txt <<'PY'
import collections
L={"Q_dcdc2c":385,"Q_dcdc2":420,"Q_dcdc2b":454,"Q_dclk3":721}
best={}
for l in open("bx.tsv"):
    q,s,pid,ln,bs=l.split(); bs=float(bs)
    if q not in best or bs>best[q][1]: best[q]=(s,bs,float(pid))
c=collections.Counter(); pid=collections.defaultdict(list)
for q,(s,bs,p) in best.items(): c[s]+=1; pid[s].append(p)
print("## I. furcatus: dcdc2c y controles (lecturas/100 aa)")
for g,ln in L.items():
    n=c.get(g,0); m=sum(pid[g])/len(pid[g]) if pid[g] else 0
    print(f"{g[2:]:8s} reads={n:5d}  reads/100aa={100*n/ln:6.2f}  pid_medio={m:5.1f}")
otros={k:v for k,v in c.items() if k not in L}
if otros: print("otros genes de la familia:",dict(sorted(otros.items(),key=lambda kv:-kv[1])[:6]))
PY
cat dcdc2c_results.txt
