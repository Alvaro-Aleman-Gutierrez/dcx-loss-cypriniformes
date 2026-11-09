#!/usr/bin/env python3
"""Cuenta, en el Iso-Seq completo de P. micromegethes (SRR11971539), las lecturas que salen del ultimo exon compartido de
dclk1a (donador 205548 en la region) hacia el exon DCL (aceptor ~205685) o hacia el primer exon de cinasa (aceptor 208526).
Usa las coordenadas de rnaseq_dcl/isoseq_conteos.json y actualiza sus conteos. Tolerancia de +-15 nt por errores PacBio."""
import re, json, sys, collections
D = "${DATA}/rnaseq_dcl/"
J = json.load(open(D + "isoseq_conteos.json"))
DON = J["exones_compartidos"][-1][1]; ACC_DCL = J["exon_dcl"][0]; ACC_KIN = J["exones_cinasa"][0][0]; TOL = 15
def intrones(pos, cig):
    x, out = pos, []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op == "N": out.append((x - 1, x + n)); x += n   # (ultima base del exon, primera base del siguiente)
        elif op in "M=XD": x += n
    return out
c = collections.Counter(); utr = 0
for l in open(D + "pmic_dclk1a.sam"):
    if l[0] == "@": continue
    f = l.split("\t"); fl = int(f[1])
    if fl & 0x904 or int(f[4]) < 20: continue
    p, cig = int(f[3]), f[5]
    for a, b in intrones(p, cig):
        if abs(a - DON) <= TOL:
            if abs(b - ACC_DCL) <= TOL:
                c["dcl"] += 1
                if not any(abs(b2 - ACC_KIN) <= TOL for _, b2 in intrones(p, cig)): utr += 1
            elif abs(b - ACC_KIN) <= TOL: c["largo"] += 1
            else: c["otro:%d" % b] += 1
print("donador", DON, "aceptores", ACC_DCL, ACC_KIN); print(dict(c), "dcl sin cinasa:", utr)
if "--escribir" in sys.argv:
    J.update(largo=c["largo"], dcl=c["dcl"], dcl_termina_en_utr=utr, fraccion_procesada=1.0,
             nota="Mapeo completo (9.962.155 lecturas, minimap2 -ax splice, servidor). Conteo de uniones desde el ultimo exon compartido.")
    json.dump(J, open(D + "isoseq_conteos.json", "w")); print("json actualizado")
