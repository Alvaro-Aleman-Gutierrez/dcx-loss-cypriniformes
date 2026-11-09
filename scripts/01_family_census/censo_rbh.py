#!/usr/bin/env python3
"""Censo v2 de la familia DCX/DCLK con ORTOLOGIA VERIFICADA, para resolver dos problemas del censo v1:
 (1) dclk2b sale en mosaico ("ausente" en 3 Danio y en las 2 Nemacheilidae, "pseudo?" en otras 4):
     un patron asi dentro de un mismo genero suele indicar fallo del metodo, no 6 perdidas reales.
 (2) 44 entradas "pseudo?/parcial" nunca revisadas.
Diferencias clave frente al v1: busqueda tblastn en el GENOMA COMPLETO (el v1 usaba miniprot sobre
regiones extraidas, que puede perder el gen si la extraccion falla) y verificacion por mejor hit
reciproco contra los proteomas de Pimephales/Gobio/Danio + las 16 proteinas de la familia.
"""
import os, subprocess, collections, re, glob
P = "${PROJECT}"
S = "${DATA}/syn2"
O = "${DATA}/censo2"; os.makedirs(O, exist_ok=True)
B2 = "${DATA}/p1b/genomes"
T = os.environ.get("THREADS", "4")
QF = f"{P}/census/queries.faa"
GEN = sorted({os.path.basename(f)[:-4] for f in glob.glob(f"{P}/genomes/*.fna")} |
             {os.path.basename(f)[:-4] for f in glob.glob(f"{B2}/*.fna")})
def run(c, **k): return subprocess.run(c, check=True, capture_output=True, text=True, **k).stdout
def dbof(sp):
    for c in (f"{P}/genomes/db_{sp}", f"{B2}/db_{sp}"):
        if os.path.exists(c + ".ndb") or os.path.exists(c + ".nal"): return c
    return None
cache = {}; cf = f"{O}/rbh_cache.tsv"
if os.path.exists(cf):
    for l in open(cf):
        a, b = l.rstrip("\n").split("\t"); cache[a] = b
def vuelta(ids):
    falta = [i for i in ids if i not in cache]
    if falta:
        open(f"{O}/_s.ids", "w").write("\n".join(falta) + "\n")
        open(f"{O}/_s.faa", "w").write(run(["seqkit", "grep", "-f", f"{O}/_s.ids", f"{S}/refprot.faa"]))
        got = {}
        for l in run(["blastp", "-query", f"{O}/_s.faa", "-db", f"{S}/qdb", "-outfmt",
                      "6 qseqid sseqid bitscore", "-max_target_seqs", "1", "-max_hsps", "1",
                      "-num_threads", T]).splitlines():
            a, b, _ = l.split("\t")
            if a not in got: got[a] = b
        with open(cf, "a") as fo:
            for i in falta: cache[i] = got.get(i, "NINGUNO"); fo.write(f"{i}\t{cache[i]}\n")
    return {i: cache[i] for i in ids}
genes = [l[1:].split()[0] for l in open(QF) if l.startswith(">")]
qlen = {}
for l in run(["seqkit", "fx2tab", "-nl", QF]).splitlines():
    p = l.split("\t"); qlen[p[0].split()[0]] = int(p[-1])
tabla = {}
det = open(f"{O}/censo2_detalle.txt", "w")
for sp in GEN:
    db = dbof(sp)
    if not db: continue
    tf = f"{O}/{sp}.tbn.tsv"
    if not os.path.exists(tf) or not os.path.getsize(tf):
        open(tf, "w").write(run(["tblastn", "-query", QF, "-db", db, "-evalue", "1e-8", "-seg", "no",
            "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore sseq",
            "-num_threads", T, "-max_target_seqs", "20", "-max_hsps", "5"]))
    por = collections.defaultdict(list)
    for l in open(tf):
        f_ = l.rstrip("\n").split("\t")
        if len(f_) < 10: continue
        por[(f_[0], f_[1])].append(f_)
    loci = []
    for (q, c), hs in por.items():
        hs.sort(key=lambda x: min(int(x[6]), int(x[7])))
        grp = []
        for h in hs:
            if grp and min(int(h[6]), int(h[7])) - max(int(grp[-1][6]), int(grp[-1][7])) < 100000: grp.append(h)
            else:
                if grp: loci.append((q, c, grp))
                grp = [h]
        if grp: loci.append((q, c, grp))
    faa = f"{O}/{sp}.loci.faa"; meta = []
    with open(faa, "w") as fo:
        for i, (q, c, grp) in enumerate(loci):
            mejor = max(grp, key=lambda x: float(x[8]))
            cov = len({x for h in grp for x in range(int(h[4]), int(h[5]) + 1)}) / qlen[q]
            fo.write(f">L{i}\n{mejor[9].replace('-','').replace('*','X')}\n")
            meta.append(dict(i=i, q=q, c=c, bs=sum(float(x[8]) for x in grp), pid=float(mejor[2]),
                             cov=cov, n=len(grp), pos=min(int(h[6]) for h in grp)))
    if not meta: continue
    mr = {}
    for l in run(["blastp", "-query", faa, "-db", f"{S}/refdb", "-outfmt", "6 qseqid sseqid bitscore",
                  "-max_target_seqs", "1", "-max_hsps", "1", "-num_threads", T]).splitlines():
        a, b, _ = l.split("\t")
        if a not in mr: mr[a] = b
    vb = vuelta(sorted(set(mr.values())))
    mejorporgen = {}
    for m in meta:
        R = mr.get(f"L{m['i']}"); back = vb.get(R, "NINGUNO") if R else "NINGUNO"
        m["R"] = R or "-"; m["back"] = back; m["ok"] = (back == m["q"])
        if m["ok"] and (m["q"] not in mejorporgen or m["bs"] > mejorporgen[m["q"]]["bs"]):
            mejorporgen[m["q"]] = m
    fila = {}
    for g in genes:
        m = mejorporgen.get(g)
        if not m: fila[g] = "ausente"
        else:
            est = "presente" if m["cov"] >= 0.5 else "parcial"
            fila[g] = f'{est}({m["cov"]:.2f};{m["pid"]:.0f})'
    tabla[sp] = fila
    det.write(f"== {sp}\n")
    for g in genes:
        m = mejorporgen.get(g)
        if m: det.write(f"   {g:12s} {fila[g]:22s} {m['c']}:{m['pos']} HSPs={m['n']} bit={m['bs']:.0f} RBH={m['R']}\n")
        else:
            rech = [x for x in meta if x["q"] == g]
            if rech:
                b = max(rech, key=lambda x: x["bs"])
                det.write(f"   {g:12s} AUSENTE tras RBH; mejor candidato bit={b['bs']:.0f} pid={b['pid']:.0f} -> RBH dio {b['back']}\n")
            else: det.write(f"   {g:12s} AUSENTE (sin ningun HSP)\n")
    det.flush(); print("hecho:", sp, flush=True)
det.close()
with open(f"{O}/censo2_table.tsv", "w") as fo:
    fo.write("especie\t" + "\t".join(g[2:] if g.startswith("Q_") else g for g in genes) + "\n")
    for sp in sorted(tabla): fo.write(sp + "\t" + "\t".join(tabla[sp][g] for g in genes) + "\n")
print("FIN censo2", flush=True)
