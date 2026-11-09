#!/usr/bin/env python3
"""Cuenta uniones de empalme de dclk1a en pez cebra (dre_R1.paf). Gen en hebra menos de la ventana:
exones DC 642-586 kb, exon DCL 545.494-545.747 kb, kinasa 530.7-514.6 kb (coordenadas de ventana, tBLASTn/dc-megablast).
Una union = operacion N del CIGAR; extremo bajo del intron = aceptor (hebra menos)."""
import re, collections, sys
O = "${DATA}/rnaseq/dre_expr"
DCL = (545494, 545747); KIN1 = (530526, 530702)
j = collections.Counter(); exon = collections.Counter(); n = int(open(f"{O}/n_R1.txt").read())
for l in open(f"{O}/dre_R1.paf"):
    x = l.rstrip("\n").split("\t")
    if int(x[11]) < 1: continue
    t = int(x[7]); cg = next((c[5:] for c in x[12:] if c.startswith("cg:Z:")), "")
    blocks = []
    for k, op in re.findall(r"(\d+)([MIDN])", cg):
        k = int(k)
        if op == "N": j[(t, t + k)] += 1; blocks.append(("N", t, t + k)); t += k
        elif op in "MD": blocks.append(("M", t, t + k)); t += k
    for b in blocks:
        if b[0] == "M":
            if b[1] < DCL[1] and b[2] > DCL[0]: exon["DCL"] += 1; break
    for b in blocks:
        if b[0] == "M" and b[1] < KIN1[1] and b[2] > KIN1[0]: exon["kin1"] += 1; break
print("lecturas R1", n)
print("lecturas en exon DCL", exon["DCL"], " en primer exon kinasa", exon["kin1"])
cerca = lambda p, r: r[0] - 30 <= p <= r[1] + 30
print("\nuniones con aceptor (extremo bajo) en exon DCL:")
for (a, b), c in sorted(j.items()):
    if cerca(b, DCL) or cerca(a, DCL): print(f"  intron {a}-{b} ({b-a} nt)  lecturas={c}")
print("\nuniones con el primer exon de kinasa:")
for (a, b), c in sorted(j.items()):
    if cerca(b, KIN1) or cerca(a, KIN1): print(f"  intron {a}-{b} ({b-a} nt)  lecturas={c}")
print("\nuniones mas frecuentes del locus (>=5):")
for (a, b), c in sorted(j.items(), key=lambda z: -z[0][1]):
    if c >= 5 and 514000 < a < 643000: print(f"  {a}-{b}  {c}")
