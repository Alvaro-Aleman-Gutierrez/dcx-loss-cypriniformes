#!/usr/bin/env python3
"""Expresion de la isoforma DCL de dclk1a en Paedocypris micromegethes (sin DCX), con Iso-Seq PacBio de cuerpo entero
(SRR11971539) mapeado con minimap2 -ax splice a la region CM051494.1:250.000-750.000.
Coordenadas en la region (posicion genomica - 249.999):
  exones compartidos (DC):    194.295-202.125 (aa 1-303 por tBLASTn)
  exon terminal DCL:          205.533-206.120 (cola EMDKVIFNSTA* + 3'UTR)
  exones de la cinasa (largo): 208.526-212.421
Clasificacion de cada lectura empalmada que cubre el locus:
  DCL   = tiene un bloque alineado dentro del exon DCL, procede de un intron que empieza en los exones DC
          y NO continua hacia los exones de cinasa
  largo = salta del tramo DC a la cinasa sin pasar por el exon DCL
"""
import re, sys, collections
SAM = sys.argv[1] if len(sys.argv) > 1 else "${DATA}/rnaseq_dcl/pmic_dclk1a.sam"
DC = (194000, 204500); DCL = (205533, 206120); KIN = (208400, 212500)
def bloques(pos, cig):
    b, x = [], pos; ini = pos
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
        n = int(n)
        if op in "M=XD": x += n
        elif op == "N": b.append((ini, x - 1)); x += n; ini = x
    b.append((ini, x - 1)); return b
cuenta = collections.Counter(); ejemplos = collections.defaultdict(list)
for l in open(SAM):
    if l.startswith("@"): continue
    c = l.split("\t"); fl = int(c[1])
    if fl & 4 or fl & 256 or fl & 2048 or int(c[4]) < 20: continue
    bl = bloques(int(c[3]), c[5])
    en = lambda r: any(a <= r[1] and b >= r[0] for a, b in bl)
    dc, dcl, kin = en(DC), en(DCL), en(KIN)
    if not (dc or dcl or kin): continue
    empalme_a_dcl = any(b0[1] < DC[1] and DCL[0] - 50 <= b1[0] <= DCL[1] for b0, b1 in zip(bl, bl[1:]))
    if dcl and empalme_a_dcl and not kin: k = "DCL (DC -> exon DCL, termina)"
    elif dcl and not kin and not dc: k = "solo exon DCL (sin empalme visible)"
    elif dc and kin and not dcl: k = "largo (DC -> cinasa)"
    elif dc and kin and dcl: k = "cruza DCL y cinasa (retencion/otro)"
    elif kin and not dc: k = "solo cinasa (corto/CPG16 o parcial)"
    else: k = "solo DC (parcial)"
    cuenta[k] += 1
    if len(ejemplos[k]) < 3: ejemplos[k].append((c[0], c[3], c[5][:80]))
print(f"Lecturas en dclk1a (MAPQ>=20): {sum(cuenta.values())}")
for k, n in cuenta.most_common(): print(f"  {n:6d}  {k}")
for k, ej in ejemplos.items():
    for e in ej: print(f"     ej. {k[:20]:20s} {e[0][:40]} pos={e[1]} {e[2]}")
