#!/usr/bin/env python3
"""Controles de atraccion de ramas largas (LBA) para la posicion de Paedocypris.
Nuestro ML en codones lo pone hermano de todos los Cypriniformes con rama troncal 0,2031 (2x la mayor
del grupo interno) = patron tipico de LBA. Si la colocacion basal es artefacto, el rechazo de la
topologia D tampoco es evidencia independiente. Tres controles:
  (a) AA  : aminoacidos (las sustituciones sinonimas saturadas son la fuente principal de LBA en codones)
  (b) NOOG: sin el grupo externo Astyanax (quita el atractor)
  (c) SLOW: eliminando el 25% de genes de evolucion mas rapida
En cada uno: posicion de Paedocypris + prueba AU frente a la topologia D.
"""
import glob, os, subprocess, collections, re
from Bio import SeqIO
from Bio.Seq import Seq
D = "${DATA}/p5/relax"
O = "${DATA}/arbol_ctrl"; os.makedirs(O, exist_ok=True)
T = os.environ.get("THREADS", "5"); IQ = "/usr/bin/iqtree2"
seqs, longs, usados = collections.defaultdict(dict), {}, []
for f in sorted(glob.glob(f"{D}/*.codon.fas")):
    g = os.path.basename(f).replace(".codon.fas", "")
    d = {r.id: str(r.seq) for r in SeqIO.parse(f, "fasta")}
    if len(d) < 12: continue
    L = len(next(iter(d.values())))
    if any(len(v) != L for v in d.values()) or L % 3: continue
    usados.append(g); longs[g] = L
    for sp, s in d.items(): seqs[sp][g] = s
taxones = sorted(seqs)
print(f"genes={len(usados)} taxones={len(taxones)}", flush=True)
# divergencia media por gen (proxy de tasa) para el control SLOW
def diverg(g):
    ss = [seqs[sp][g] for sp in taxones if g in seqs[sp]]
    if len(ss) < 4: return 9.9
    tot = n = 0
    for i in range(len(ss)):
        for j in range(i + 1, len(ss)):
            a, b = ss[i], ss[j]
            pares = [(x, y) for x, y in zip(a, b) if x not in "-N" and y not in "-N"]
            if len(pares) < 60: continue
            tot += sum(x != y for x, y in pares) / len(pares); n += 1
    return tot / n if n else 9.9
dv = {g: diverg(g) for g in usados}
lentos = sorted(usados, key=lambda g: dv[g])[:int(len(usados) * 0.75)]
print(f"SLOW: {len(lentos)} genes mas lentos (divergencia {dv[lentos[0]]:.3f}-{dv[lentos[-1]]:.3f}); "
      f"rapidos excluidos hasta {max(dv.values()):.3f}", flush=True)
def escribir(nombre, genes, taxa, aa):
    with open(f"{O}/{nombre}.fas", "w") as fo:
        for sp in taxa:
            partes = []
            for g in genes:
                s = seqs[sp].get(g)
                if s is None: partes.append(("-" * (longs[g] // 3)) if aa else ("-" * longs[g])); continue
                if aa:
                    p = str(Seq(s.replace("-", "N")).translate()).replace("*", "X")
                    p = "".join(("-" if s[i*3:i*3+3] == "---" else c) for i, c in enumerate(p))
                    partes.append(p)
                else: partes.append(s)
            fo.write(f">{sp}\n" + "".join(partes) + "\n")
    ini = 1
    with open(f"{O}/{nombre}.nex", "w") as fo:
        fo.write("#nexus\nbegin sets;\n"); mod = "LG+G" if aa else "GTR+G"
        for g in genes:
            L = longs[g] // 3 if aa else longs[g]
            fo.write(f"  charset {g} = {ini}-{ini+L-1};\n"); ini += L
        fo.write("  charpartition mine = " + ", ".join(f"{mod}:{g}" for g in genes) + ";\nend;\n")
    return ini - 1
PAE = [t for t in taxones if t.startswith("Paedocypris")]
DAN = [t for t in taxones if t.startswith("Danio_")]
DNL = [t for t in taxones if t.startswith("Danionella")]
def cl(xs): return "(" + ",".join(xs) + ")" if len(xs) > 1 else xs[0]
res = []
for nombre, genes, taxa, aa, modelo in [
        ("a_AA",   usados, taxones, True,  "LG+F+G4"),
        ("b_NOOG", usados, [t for t in taxones if t != "Astyanax_mexicanus"], False, "GTR+F+G4"),
        ("c_SLOW", lentos, taxones, False, "GTR+F+G4")]:
    n = escribir(nombre, genes, taxa, aa)
    print(f"\n### {nombre}: {n} posiciones, {len(taxa)} taxones, modelo {modelo}", flush=True)
    resto = [t for t in taxa if t not in PAE + DAN + DNL]
    topoD = f"(({cl(PAE)},{cl(DAN)}),{cl(DNL)},{cl(resto)});"
    open(f"{O}/{nombre}_topoD.nwk", "w").write(topoD + "\n")
    base = [IQ, "-s", f"{O}/{nombre}.fas", "-p", f"{O}/{nombre}.nex", "-m", modelo, "-T", T, "-redo"]
    raiz = ["-o", "Astyanax_mexicanus"] if "Astyanax_mexicanus" in taxa else []
    subprocess.run(base + ["-B", "1000", "--prefix", f"{O}/{nombre}_ml"] + raiz, check=False)
    subprocess.run(base + ["-g", f"{O}/{nombre}_topoD.nwk", "--prefix", f"{O}/{nombre}_mlD"] + raiz, check=False)
    tf, tfd = f"{O}/{nombre}_ml.treefile", f"{O}/{nombre}_mlD.treefile"
    if os.path.exists(tf) and os.path.exists(tfd):
        with open(f"{O}/{nombre}_cands.trees", "w") as fo:
            fo.write(open(tf).read()); fo.write(open(tfd).read())
        subprocess.run(base + ["-z", f"{O}/{nombre}_cands.trees", "-n", "0", "-zb", "10000", "-zw", "-au",
                               "--prefix", f"{O}/{nombre}_test"], check=False)
    # leer resultado
    au = "-"
    rep = f"{O}/{nombre}_test.iqtree"
    if os.path.exists(rep):
        filas = [l.split() for l in open(rep) if re.match(r"^\s+[12] -?\d", l)]
        if len(filas) >= 2: au = f"ML p-AU={filas[0][-2]}  topoD p-AU={filas[1][-2]}"
    t = open(tf).read() if os.path.exists(tf) else ""
    m = re.search(r"Paedocypris\w*:([\d.]+)", t)
    basal = bool(re.search(r"\)[^(),;]*;?\s*$", t)) and ("Paedocypris" in t.split(",")[-1] or t.rstrip().endswith(");"))
    # posicion: hermana del resto si el clado Paedocypris cuelga del nodo mas externo tras la raiz
    res.append((nombre, n, len(taxa), au, t))
    print(f"   {nombre}: {au}", flush=True)
with open(f"{O}/controles_results.txt", "w") as fo:
    for nombre, n, nt, au, t in res:
        fo.write(f"== {nombre}  posiciones={n} taxones={nt}\n   {au}\n   arbol: {t}\n")
print("FIN controles", flush=True)
