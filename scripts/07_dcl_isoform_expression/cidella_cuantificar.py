#!/usr/bin/env python3
"""Cuenta lecturas (R1, primarias, MAPQ>=0, alineamiento >=60 nt) por region en los PAF de cidella_expresion.sh.
Regiones (coordenadas 0-based sobre los mRNA RefSeq; ver blastn X3 vs X1):
  DC compartido   X1 0-1187 / X3 66-1253  (isoforma larga + DCL)
  kinasa          X1 >=1188               (larga y solo-kinasa)
  DCL especifico  X3 >=1254               (exon terminal alternativo + 3'UTR)
  union DCL       lecturas en X3 que cruzan 1253 (unen exon DC compartido con exon DCL)
  union larga     lecturas en X1 que cruzan 1187
Salida: cid_expr/expresion_tejidos.tsv con conteos y RPKM (por millon de R1)."""
import glob, os, collections
O = "${DATA}/rnaseq/cid_expr"
L = {"DCX": 1071, "dclk2a": 2935, "actb1": 1781, "elavl3": 6151, "gfap": 2562}
LX1, LX3 = 3578, 5337
filas = []
for ok in sorted(glob.glob(f"{O}/*.paf.ok")):
    paf = ok[:-3]; t, r = os.path.basename(paf).split(".")[:2]
    n = int(open(ok).read()); c = collections.Counter()
    for l in open(paf):
        x = l.split("\t"); tn = x[5]; ts, te = int(x[7]), int(x[8])
        if int(x[10]) < 60: continue
        if tn == "dclk1a_X1_long":
            if te <= 1188: c["DC_compartido"] += 1
            else: c["kinasa"] += 1
            if ts < 1180 and te > 1196: c["union_larga"] += 1
        elif tn == "dclk1a_X3_DCL":
            if te <= 1254: c["DC_compartido"] += 1
            else: c["DCL_especifico"] += 1
            if ts < 1246 and te > 1262: c["union_DCL"] += 1
        else: c[tn] += 1
    f = dict(tejido=t, corrida=r, lecturas_R1=n, **{k: c[k] for k in ["DCX", "DC_compartido", "kinasa", "DCL_especifico", "union_larga", "union_DCL", "dclk2a", "elavl3", "gfap", "actb1"]})
    rp = lambda k, ln: round(c[k] / (ln / 1e3) / (n / 1e6), 2)
    for k, ln in L.items(): f[f"RPKM_{k}"] = rp(k, ln)
    f["RPKM_kinasa"] = rp("kinasa", LX1 - 1188); f["RPKM_DCL"] = rp("DCL_especifico", LX3 - 1254)
    u = c["union_larga"] + c["union_DCL"]; f["frac_DCL_uniones"] = round(c["union_DCL"] / u, 3) if u else ""
    filas.append(f)
import pandas as pd
df = pd.DataFrame(filas); df.to_csv(f"{O}/expresion_tejidos.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(df.to_string(index=False))
