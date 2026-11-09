#!/usr/bin/env python3
"""Estratifica la nula de compactacion por numero de genes de Gobio intercalados entre las dos anclas
(0 = intervalo puramente intergenico; >=1 = intervalo con genes sin ancla, perdidos o ambiguos)."""
import re, sys, os, collections, statistics, bisect
D = "${DATA}"; O = f"{D}/compactacion"
pos = collections.defaultdict(list); gpos = {}
for l in open(f"{D}/syn2/gff/Gobio.gff"):
    c = l.split("\t")
    if len(c) < 9 or c[2] != "gene" or "gene_biotype=protein_coding" not in c[8]: continue
    g = re.search(r"gene=([^;\n]+)", c[8]).group(1); pos[c[0]].append(int(c[3])); gpos[g] = (c[0], int(c[3]), int(c[4]))
for ch in pos: pos[ch].sort()
def entre(ch, a, b):   # genes de Gobio que empiezan estrictamente entre a y b
    p = pos[ch]; return bisect.bisect_left(p, b) - bisect.bisect_right(p, a)
for sp in sys.argv[1:]:
    filas = [l.rstrip().split("\t") for l in open(f"{O}/{sp}_intervalos.tsv")][1:]
    est = collections.defaultdict(list)
    for ch, g1, g2, gg, tg, r in filas:
        a = gpos[g1][2] if g1 in gpos else None; b = gpos[g2][1] if g2 in gpos else None
        if a is None or b is None: continue
        gg = float(gg); n = entre(ch, gpos[g1][1], gpos[g2][1]) - 0
        if not (128e3 <= gg <= 515e3): continue
        est["0" if n <= 0 else ("1" if n == 1 else "2+")].append(float(r))
    r0 = 5.732 / 257.264 if "micromegethes" in sp else (5.842 / 257.264 if "Singkep" in sp else None)
    print(f"== {sp}  (intervalos de 128-515 kb en Gobio; intervalo DCX: 257 kb, 2 genes intercalados)")
    for k in ("0", "1", "2+"):
        v = sorted(est[k])
        if not v: continue
        linea = f"  {k:3s} genes intercalados: n={len(v):4d} mediana={statistics.median(v):.3f} p05={v[int(.05*len(v))]:.3f} p01={v[int(.01*len(v))]:.3f}"
        if r0: linea += f"  percentil DCX={100*sum(x<=r0 for x in v)/len(v):.1f}"
        print(linea)

# ---- tabla para la Fig. 3B: cocientes de intervalos de 128-515 kb en Gobio, por numero de genes intercalados
def exportar(especies, salida=f"{O}/nula_estratificada.tsv"):
    with open(salida, "w") as fo:
        fo.write("especie\tn_intercalados\thueco_gobio\thueco_diana\tcociente\n")
        for sp in especies:
            for ch, g1, g2, gg, tg, r in [l.rstrip().split("\t") for l in open(f"{O}/{sp}_intervalos.tsv")][1:]:
                if g1 not in gpos or g2 not in gpos or not (128e3 <= float(gg) <= 515e3): continue
                n = entre(ch, gpos[g1][1], gpos[g2][1])
                fo.write(f"{sp}\t{max(n,0)}\t{gg}\t{tg}\t{r}\n")
    print("escrito", salida)
if os.environ.get("EXPORTAR"): exportar(sys.argv[1:])
