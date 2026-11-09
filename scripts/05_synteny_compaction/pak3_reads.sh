#!/bin/bash
# Verificar pak3 (= pak3loc, el unico vecino verificado de DCX) en lecturas crudas.
# Base: el RBH de sintenia validó Pimephales XP_039517591.2 como el vecino de DCX.
# Dato clave de las anotaciones: Pimephales y Gobio tienen UNA copia de pak3; Danio rerio NO tiene pak3.
# Por eso la presencia/ausencia de pak3 en lecturas es informativa sin ambiguedad de parálogos.
# PRECAUCION: la base discriminante lleva las 69 proteinas PAK de los 3 proteomas + las 16 de la familia DCX,
# para que un fragmento de pak1/pak2 no se cuente como pak3 (mismo error que nos engañó 3 veces con el dominio DC).
set -uo pipefail
P=${PROJECT}; DD=${DATA}; O=$DD/pak3
T=${THREADS:-4}; LOG=$O/pak3.log; OUT=$O/pak3_results.txt
: > $OUT; echo "[$(date)] inicio" > $LOG
# consultas: pak3 + controles de cobertura de la familia DCX
python3 -I - <<PY
from Bio import SeqIO
keep={"Q_dclk1a","Q_dclk2a","Q_dclk3"}
R=[r for r in SeqIO.parse("$P/census/queries.faa","fasta") if r.id in keep]
R+=list(SeqIO.parse("$O/q_pak3.faa","fasta"))
SeqIO.write(R,"$O/q.faa","fasta")
PY
for item in "reads_cerebrum_ctrl:D_cerebrum_CONTROL(tiene DCX)" "reads_dracula:D_dracula" \
            "reads_paedocypris:P_micromegethes" "reads_misgurnus:M_anguillicaudatus" \
            "reads_paramisgurnus:P_dabryanus"; do
  IFS=: read dir tag <<< "$item"; RD=$DD/$dir; W=$O/$dir; mkdir -p $W
  ls $RD/reads_db*.nsq >/dev/null 2>&1 || { echo "## $tag: SIN base de lecturas" >> $OUT; continue; }
  echo "[$(date)] $tag tblastn" >> $LOG
  [ -s $W/tbn.tsv ] || nice -n 15 tblastn -query $O/q.faa -db $RD/reads_db -evalue 1e-5 -seg no \
      -outfmt "6 qseqid sseqid pident length" -max_target_seqs 200000 -num_threads $T > $W/tbn.tsv 2>>$LOG
  cut -f2 $W/tbn.tsv | sort -u > $W/ids.txt
  [ -s $W/reads.fna ] || seqkit grep -j $T -f $W/ids.txt $(ls $RD/*.fasta) | seqkit rename > $W/reads.fna 2>>$LOG
  echo "[$(date)] $tag blastx" >> $LOG
  [ -s $W/bx.tsv ] || nice -n 15 blastx -query $W/reads.fna -db $O/discdb -evalue 1e-5 \
      -outfmt "6 qseqid sseqid pident length bitscore" -max_target_seqs 5 -max_hsps 1 -num_threads $T > $W/bx.tsv 2>>$LOG
  python3 -I - "$tag" "$W/bx.tsv" "$O/discdb.faa" >> $OUT <<'PY'
import sys,collections,re
tag,fn,dbf=sys.argv[1],sys.argv[2],sys.argv[3]
desc={}
for l in open(dbf):
    if l.startswith(">"):
        i=l[1:].split()[0]; desc[i]=l[1:].rstrip()
def grupo(s):
    d=desc.get(s,s)
    m=re.search(r"kinase PAK (\d[ab]?)",d)
    if m: return "pak"+m.group(1)
    if "p21-activated" in d: return "pak_otros"
    if s.startswith("Q_"): return s[2:]
    return "otro"
best={}
for l in open(fn):
    q,s,pid,ln,bs=l.split(); bs=float(bs)
    if q not in best or bs>best[q][1]: best[q]=(s,bs,float(pid))
c=collections.Counter(); pids=collections.defaultdict(list)
for q,(s,bs,pid) in best.items():
    g=grupo(s); c[g]+=1; pids[g].append(pid)
L={"pak3":552,"pak1":545,"pak2a":524,"pak2b":524,"dclk1a":770,"dclk2a":826,"dclk3":721}
print(f"## {tag}")
for g in ["pak3","pak1","pak2a","pak2b","dclk1a","dclk2a","dclk3"]:
    n=c.get(g,0); ln=L.get(g,500); m=sum(pids[g])/len(pids[g]) if pids[g] else 0
    print(f"{g:8s} reads={n:5d}  reads/100aa={100*n/ln:6.2f}  pid_medio={m:5.1f}")
otros={k:v for k,v in c.items() if k not in L}
if otros: print("otros:",dict(sorted(otros.items(),key=lambda kv:-kv[1])[:5]))
print()
PY
  echo "[$(date)] $tag listo" >> $LOG
done
echo "FIN pak3" >> $OUT; echo "[$(date)] COMPLETO" >> $LOG
