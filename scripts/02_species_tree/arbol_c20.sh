#!/bin/bash
# Control de LBA con modelo heterogeneo por sitio: LG+C20+F+G (PMSF) en aminoacidos, matriz final de 24 taxones
# (arbol_ext/super.fas traducido por gen). ML libre + ML restringido a topologia D + prueba AU (10.000 RELL).
set -u
D=${DATA}; O=$D/c20; T=${THREADS:-6}; cd $O
python3 - <<'PY'
import re
from Bio import SeqIO
from Bio.Seq import Seq
D="${DATA}"
parts=[tuple(map(int,m.groups())) for m in re.finditer(r"charset \S+ = (\d+)-(\d+);", open(f"{D}/arbol_ext/partitions.nex").read())]
recs=list(SeqIO.parse(f"{D}/arbol_ext/super.fas","fasta"))
with open("aa.fas","w") as fo:
    for r in recs:
        s=str(r.seq).upper(); out=[]
        for a,b in parts:
            g=s[a-1:b]
            for i in range(0,len(g)-2,3):
                c=g[i:i+3]
                if c=="---": out.append("-")
                elif "-" in c or "N" in c: out.append("X")
                else:
                    p=str(Seq(c).translate()); out.append("X" if p=="*" else p)
        fo.write(f">{r.id}\n{''.join(out)}\n")
taxa=[r.id for r in recs]
PAE=[t for t in taxa if t.startswith("Paedocypris")]; DAN=[t for t in taxa if t.startswith("Danio_")]
DNL=[t for t in taxa if t.startswith("Danionella")]; resto=[t for t in taxa if t not in PAE+DAN+DNL]
cl=lambda xs:"("+",".join(xs)+")" if len(xs)>1 else xs[0]
open("topoD.nwk","w").write(f"(({cl(PAE)},{cl(DAN)}),{cl(DNL)},{cl(resto)});\n")
print("taxa",len(taxa),"sitios",len(out))
PY
IQ=/usr/bin/iqtree2
# 1) arbol guia LG+F+G4
$IQ -s aa.fas -m LG+F+G4 -T $T --prefix guia -redo -quiet
# 2) PMSF libre con UFBoot
$IQ -s aa.fas -m LG+C20+F+G -ft guia.treefile -B 1000 -T $T --prefix c20_ml -redo -quiet
# 3) PMSF restringido a topologia D (arbol guia tambien restringido)
$IQ -s aa.fas -m LG+F+G4 -g topoD.nwk -T $T --prefix guiaD -redo -quiet
$IQ -s aa.fas -m LG+C20+F+G -ft guiaD.treefile -g topoD.nwk -T $T --prefix c20_mlD -redo -quiet
# 4) AU con el modelo PMSF del arbol libre
cat c20_ml.treefile c20_mlD.treefile > cands.trees
$IQ -s aa.fas -m LG+C20+F+G -ft guia.treefile -z cands.trees -n 0 -zb 10000 -zw -au -T $T --prefix c20_test -redo -quiet
grep -A6 "Tree      logL" c20_test.iqtree
echo FIN_C20
