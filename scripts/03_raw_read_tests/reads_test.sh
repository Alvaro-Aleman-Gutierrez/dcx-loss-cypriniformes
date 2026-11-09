#!/bin/bash
# Prueba de presencia/ausencia de DCX en lecturas WGS crudas.
# uso: reads_test.sh <dir_con_reads_db> <etiqueta>
# 1) tBLASTn de DCX (Carassius) + controles monocopia (dclk1a, dclk2a, dclk3 de Danio) contra lecturas
# 2) BLASTx de cada lectura con hit contra la familia completa -> asignación por mejor hit
# 3) lecturas por gen normalizadas por longitud de proteína (lecturas/100 aa)
set -euo pipefail
P=${PROJECT}
T=${THREADS:-4}
RD=$1; TAG=$2; OUT=$RD/test; mkdir -p $OUT; cd $OUT
python3 -I - <<EOF
from Bio import SeqIO
want={"Caraur_LOC113084161":"DCX","Danrer_dclk1a":"dclk1a","Danrer_dclk2a":"dclk2a","Danrer_dclk3":"dclk3"}
R=[]
for r in SeqIO.parse("$P/phylo/family_clean.faa","fasta"):
    for k,v in want.items():
        if r.id.startswith(k): r.id=v; r.description=""; R.append(r)
SeqIO.write(R,"queries.faa","fasta")
EOF
[ -s tbn.tsv ] || tblastn -query queries.faa -db $RD/reads_db -evalue 1e-5 -outfmt "6 qseqid sseqid pident length" \
  -max_target_seqs 100000 -num_threads ${T} -seg no > tbn.tsv
cut -f2 tbn.tsv | sort -u > hit_ids.txt
seqkit grep -j ${T} -f hit_ids.txt $(ls $RD/*.fasta) | seqkit rename > hit_reads.fna  # IDs de pares repetidos -> renombrados
blastx -query hit_reads.fna -db $P/tbn/famdb -evalue 1e-5 -outfmt "6 qseqid sseqid pident length bitscore" \
  -max_target_seqs 5 -max_hsps 1 -num_threads ${T} > bx.tsv
python3 -I - "$TAG" <<'EOF'
import sys,collections
tag=sys.argv[1]
DCX={"DCX","dcx","LOC103178947","LOC113084161","LOC103035707","LOC102692282","LOC102366980"}
L={"DCX":357,"dclk1a":770,"dclk2a":826,"dclk3":721}
best={}
for l in open("bx.tsv"):
    q,s,pid,ln,bs=l.split(); bs=float(bs)
    if q not in best or bs>best[q][1]: best[q]=(s,bs,float(pid))
cnt=collections.Counter(); pids=collections.defaultdict(list)
for q,(s,bs,pid) in best.items():
    parts=s.split("_"); gene=parts[1]
    if gene in DCX: g="DCX"
    elif s.startswith(("Danrer_dclk1a","Astmex_dclk1a","Lepocu_dclk1a","Calmil_dclk1a","Latcha_DCLK1A")): g="dclk1a"
    elif "dclk2a" in s or "DCLK2A" in s: g="dclk2a"
    elif "dclk3" in s.lower(): g="dclk3"
    else: g="other:"+gene
    cnt[g]+=1; pids[g].append(pid)
print(f"## {tag}: reads asignadas por mejor hit BLASTx")
for g in ["DCX","dclk1a","dclk2a","dclk3"]:
    n=cnt.get(g,0); m=sum(pids[g])/len(pids[g]) if pids[g] else 0
    print(f"{g:8s} reads={n:5d}  reads/100aa={100*n/L[g]:6.2f}  pid_medio={m:5.1f}")
oth=sum(v for k,v in cnt.items() if k.startswith("other")); print(f"otros (DCDC, dclk1b/2b, etc.) reads={oth}")
EOF
