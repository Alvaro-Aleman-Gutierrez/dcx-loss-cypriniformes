#!/usr/bin/env python3
"""Arbol de especies para la Figura 1 con la matriz EXTENDIDA (matriz_ext/aln: 50 genes, 24 taxones,
incluye Cobitidae -Misgurnus, Paramisgurnus- e Ictalurus furcatus). El arbol anterior (arbol/ml.treefile)
no tenia Cobitidae, asi que la figura solo podia mostrar 3 de las 4 perdidas.
Metodo identico a arbol_paedocypris.py: supermatriz de codones, particion por gen GTR+F+G4,
IQ-TREE2, 1000 UFBoot + 1000 SH-aLRT. Grupo externo: Astyanax + Ictalurus (Characiphysi).
"""
import glob, os, subprocess, collections
from Bio import SeqIO
D = "${DATA}/matriz_ext/aln"
O = "${DATA}/arbol_ext"; os.makedirs(O, exist_ok=True)
T = os.environ.get("THREADS", "6"); IQ = "/usr/bin/iqtree2"; MIN_TAXA = 12
ficheros = sorted(glob.glob(f"{D}/*.codon.fas"))
seqs, longs, usados = collections.defaultdict(dict), {}, []
for f in ficheros:
    g = os.path.basename(f).replace(".codon.fas", "")
    d = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    if len(d) < MIN_TAXA: continue
    L = len(next(iter(d.values())))
    if any(len(v) != L for v in d.values()): print("longitudes desiguales, se descarta", g); continue
    usados.append(g); longs[g] = L
    for sp, s in d.items(): seqs[sp][g] = s
taxones = sorted(seqs)
print(f"genes usados: {len(usados)} de {len(ficheros)}   taxones: {len(taxones)}", flush=True)
for sp in taxones: print(f"   {sp:32s} {len(seqs[sp])} genes", flush=True)
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
subprocess.run([IQ, "-s", f"{O}/super.fas", "-p", f"{O}/partitions.nex", "-m", "GTR+F+G4",
                "-B", "1000", "-alrt", "1000", "-T", T, "--prefix", f"{O}/ml", "-redo",
                "-o", "Astyanax_mexicanus,Ictalurus_furcatus"], check=False)
print("FIN arbol_ext", flush=True)
