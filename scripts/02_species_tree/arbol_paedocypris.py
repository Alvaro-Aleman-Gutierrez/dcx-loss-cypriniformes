#!/usr/bin/env python3
"""Arbol de especies con los 49 genes aleatorios del punto 5, para CONTRASTAR la topologia D
(Paedocypris hermano exclusivo del clado Danio), la unica que reduce las perdidas de DCX de 4 a 3.
Salida: supermatriz concatenada + particiones por gen, ML con IQ-TREE2 (bootstrap ultrarrapido)
y prueba topologica (AU/KH/SH y Bayes factor) de la topologia D frente a la mejor topologia ML.
"""
import glob, os, subprocess, collections
from Bio import SeqIO
D = "${DATA}/p5/relax"
O = "${DATA}/arbol"; os.makedirs(O, exist_ok=True)
T = os.environ.get("THREADS", "8")
IQ = "/usr/bin/iqtree2"
MIN_TAXA = 12          # genes con al menos este numero de taxones
ficheros = sorted(glob.glob(f"{D}/*.codon.fas"))
seqs, longs, usados = collections.defaultdict(dict), {}, []
for f in ficheros:
    g = os.path.basename(f).replace(".codon.fas", "")
    d = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    if len(d) < MIN_TAXA: continue
    L = len(next(iter(d.values())))
    if any(len(v) != L for v in d.values()): continue
    usados.append(g); longs[g] = L
    for sp, s in d.items(): seqs[sp][g] = s
taxones = sorted(seqs)
print(f"genes usados: {len(usados)} de {len(ficheros)}   taxones: {len(taxones)}", flush=True)
# supermatriz: huecos como gaps
with open(f"{O}/super.fas", "w") as fo:
    for sp in taxones:
        fo.write(f">{sp}\n" + "".join(seqs[sp].get(g, "-" * longs[g]) for g in usados) + "\n")
ini = 1
with open(f"{O}/partitions.nex", "w") as fo:
    fo.write("#nexus\nbegin sets;\n")
    for g in usados:
        fo.write(f"  charset {g} = {ini}-{ini+longs[g]-1};\n"); ini += longs[g]
    fo.write("  charpartition mine = " + ", ".join(f"GTR+G:{g}" for g in usados) + ";\nend;\n")
print(f"supermatriz: {ini-1} pb", flush=True)
# topologias a contrastar (solo los taxones presentes)
def pod(nombres): return [t for t in taxones if t in nombres]
DANIO = pod({f"Danio_{x}" for x in ["rerio","aesculapii","kyathit","albolineatus","tinwini","choprae","jaintianensis","nigrofasciatus"]})
PAEDO = pod({"Paedocypris_micromegethes","Paedocypris_carbunculus","Paedocypris_sp_Pulau_Singkep"})
DNLA  = pod({"Danionella_cerebrum","Danionella_dracula"})
def cl(xs): return "(" + ",".join(xs) + ")" if len(xs) > 1 else xs[0]
restantes = [t for t in taxones if t not in DANIO + PAEDO + DNLA]
# D: Paedocypris hermano SOLO del clado Danio
topoD = f"(({cl(PAEDO)},{cl(DANIO)}),{cl(DNLA)},{cl(restantes)});"
# A: Paedocypris fuera de Danioninae (hermano del resto, referencia)
topoA = f"({cl(PAEDO)},({cl(DNLA)},{cl(DANIO)}),{cl(restantes)});"
open(f"{O}/topo_D.nwk", "w").write(topoD + "\n")
open(f"{O}/topo_A.nwk", "w").write(topoA + "\n")
print("topoD:", topoD[:160], flush=True)
# 1) ML sin restricciones
cmd = [IQ, "-s", f"{O}/super.fas", "-p", f"{O}/partitions.nex", "-m", "GTR+F+G4",
       "-B", "1000", "-alrt", "1000", "-T", T, "--prefix", f"{O}/ml", "-redo",
       "-o", "Astyanax_mexicanus"]
print("ML sin restricciones...", flush=True)
subprocess.run(cmd, check=False)
# 2) ML con la topologia D impuesta (-g restringe la busqueda)
print("ML con topologia D impuesta...", flush=True)
subprocess.run([IQ, "-s", f"{O}/super.fas", "-p", f"{O}/partitions.nex", "-m", "GTR+F+G4",
                "-g", f"{O}/topo_D.nwk", "-T", T, "--prefix", f"{O}/mlD", "-redo",
                "-o", "Astyanax_mexicanus"], check=False)
# 3) prueba topologica: ML libre vs ML restringido a D
if os.path.exists(f"{O}/ml.treefile") and os.path.exists(f"{O}/mlD.treefile"):
    with open(f"{O}/cands.trees", "w") as fo:
        fo.write(open(f"{O}/ml.treefile").read()); fo.write(open(f"{O}/mlD.treefile").read())
    print("prueba topologica AU/KH/SH...", flush=True)
    subprocess.run([IQ, "-s", f"{O}/super.fas", "-p", f"{O}/partitions.nex", "-m", "GTR+F+G4",
                    "-z", f"{O}/cands.trees", "-n", "0", "-zb", "10000", "-zw", "-au",
                    "-T", T, "--prefix", f"{O}/test", "-redo"], check=False)
print("FIN arbol", flush=True)
