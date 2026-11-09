#!/usr/bin/env python3
"""Control de compactacion genomica para la delecion de DCX en Paedocypris.

Pregunta: el intervalo hip1rb-pak3loc mide 20-23 kb en Paedocypris y 188-326 kb en cipriniformes con DCX.
El genoma de Paedocypris es ~3,5 veces menor que el de Gobio. Se acorto ese intervalo mas de lo que
explica la compactacion general?

Metodo:
 1. Proteina mas larga de cada gen anotado de Gobio gobio (RefSeq) -> miniprot contra cada genoma diana.
 2. Ancla = proteina con hit unico (2.o hit < 0.9 x AS del 1.o) y cobertura >= 50 %.
 3. Para pares de anclas CONSECUTIVAS en el orden de Gobio que caen en la misma secuencia diana (<5 Mb):
    hueco = distancia entre extensiones CDS (fin de A -> inicio de B), en Gobio y en la diana.
    cociente = hueco_diana / hueco_Gobio. Distribucion nula = todos esos intervalos.
 4. Intervalo de interes: hip1rb -> pak3loc (en Gobio contiene LOC307125817 y DCX).
    Se calcula su cociente y su percentil frente a intervalos de Gobio de longitud comparable.
 Pimephales (con DCX) es la calibracion: alli el intervalo deberia caer en el centro de la distribucion.
"""
import os, re, sys, subprocess, collections, statistics, json
D = "${DATA}"; O = f"{D}/compactacion"; os.makedirs(O, exist_ok=True)
MP = os.path.expanduser("~/bin_extra/miniprot-0.13_x64-linux/miniprot"); T = os.environ.get("THREADS", "6")
DIANAS = sys.argv[1:] or ["Paedocypris_micromegethes", "Paedocypris_sp_Pulau_Singkep", "Pimephales_promelas"]
HIP, PAK, DCX, TM = "hip1rb", "LOC307124699", "LOC307125267", "LOC307125817"

# ---- 1. anotacion de Gobio: extension CDS por proteina; proteina mas larga por gen
cds = {}
for l in open(f"{D}/syn2/gff/Gobio.gff"):
    c = l.split("\t")
    if len(c) < 9 or c[2] != "CDS": continue
    p = re.search(r"protein_id=([^;\n]+)", c[8]); g = re.search(r"gene=([^;\n]+)", c[8])
    if not p or not g: continue
    k = p.group(1); a, b = int(c[3]), int(c[4])
    if k in cds: x = cds[k]; x[2] = min(x[2], a); x[3] = max(x[3], b); x[4] += b - a + 1
    else: cds[k] = [g.group(1), c[0], a, b, b - a + 1]
porgen = {}
for k, (g, ch, a, b, L) in cds.items():
    if g not in porgen or L > cds[porgen[g]][4]: porgen[g] = k
GEN = {porgen[g]: g for g in porgen}
faa = f"{O}/gobio_longest.faa"
if not os.path.exists(faa):
    seqs, cur = {}, None
    for l in open(f"{D}/syn2/prot/Gobio.faa"):
        if l.startswith(">"): cur = l[1:].split()[0].split("|")[-1]; seqs[cur] = []
        elif cur: seqs[cur].append(l.strip())
    with open(faa, "w") as fo:
        for k in GEN:
            if k in seqs: fo.write(f">{k}\n{''.join(seqs[k])}\n")
print(f"Gobio: {len(porgen)} genes con CDS", flush=True)
# orden de Gobio por cromosoma
orden = collections.defaultdict(list)
for k, g in GEN.items():
    _, ch, a, b, _ = cds[k]; orden[ch].append((a, b, g, k))
for ch in orden: orden[ch].sort()
G = {g: (cds[k][1], cds[k][2], cds[k][3]) for k, g in GEN.items()}
for x in (HIP, PAK, DCX): assert x in G, x

def anclas(sp):
    paf = f"{O}/{sp}.paf"
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        print(f"  miniprot -> {sp}", flush=True)
        with open(paf, "w") as fo:
            subprocess.run([MP, "-t", T, "-N", "3", "--outs=0.5", f"{D}/genomes/{sp}.fna", faa], stdout=fo,
                           stderr=subprocess.DEVNULL, check=True)
    hits = collections.defaultdict(list)
    for l in open(paf):
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t")
        if len(c) < 12: continue
        AS = next((int(x[5:]) for x in c[12:] if x.startswith("AS:i:")), 0)
        hits[c[0]].append((AS, c[5], int(c[7]), int(c[8]), (int(c[3]) - int(c[2])) / int(c[1])))
    A = {}
    for q, hs in hits.items():
        hs.sort(reverse=True)
        if hs[0][4] < 0.5: continue
        if len(hs) > 1 and hs[1][0] >= 0.9 * hs[0][0]: continue
        if q in GEN: A[GEN[q]] = hs[0][1:4]
    return A

def analizar(sp):
    A = anclas(sp)
    filas = []
    for ch, gs in orden.items():
        m = [(a, b, g) for a, b, g, k in gs if g in A]
        for (a1, b1, g1), (a2, b2, g2) in zip(m, m[1:]):
            s1, s2 = A[g1], A[g2]
            if s1[0] != s2[0]: continue
            gg = a2 - b1
            lo, hi = sorted([s1, s2], key=lambda s: s[1])
            tg = hi[1] - lo[2]
            if gg < 5000 or tg <= 0 or abs(tg) > 5e6: continue
            filas.append((ch, g1, g2, gg, tg, tg / gg))
    # intervalo de interes (puede no ser consecutivo entre anclas si TM/DCX no anclan)
    res = {"sp": sp, "n_anclas": len(A), "n_intervalos": len(filas)}
    if HIP in A and PAK in A and A[HIP][0] == A[PAK][0]:
        gch, ga1, gb1 = G[HIP]; _, ga2, gb2 = G[PAK]
        gg = max(ga1, ga2) - min(gb1, gb2)
        lo, hi = sorted([A[HIP], A[PAK]], key=lambda s: s[1]); tg = hi[1] - lo[2]
        r = tg / gg
        # nula: intervalos de Gobio de longitud comparable (0,5x-2x)
        comp = [f[5] for f in filas if 0.5 * gg <= f[3] <= 2 * gg]
        todos = [f[5] for f in filas]
        pct = 100 * sum(x <= r for x in comp) / len(comp)
        res.update(dict(gobio_kb=gg / 1e3, diana_kb=tg / 1e3, cociente=r, n_comp=len(comp),
                        mediana_comp=statistics.median(comp), p05_comp=sorted(comp)[int(0.05 * len(comp))],
                        percentil=pct, mediana_todos=statistics.median(todos), dcx_anclado=DCX in A))
    with open(f"{O}/{sp}_intervalos.tsv", "w") as fo:
        fo.write("cromosoma_gobio\tgen_a\tgen_b\thueco_gobio\thueco_diana\tcociente\n")
        for f in filas: fo.write("\t".join(map(str, f)) + "\n")
    return res

if __name__ == "__main__":
    tam = {s: sum(int(l.split()[1]) for l in open(f"{D}/genomes/{s}.fna.fai")) for s in ["Gobio_gobio"] + DIANAS}
    out = []
    for sp in DIANAS:
        r = analizar(sp); r["genoma_ratio"] = tam[sp] / tam["Gobio_gobio"]; out.append(r)
        print(json.dumps(r, ensure_ascii=False), flush=True)
    json.dump(out, open(f"{O}/compactacion_results.json", "w"), indent=1, ensure_ascii=False)
    print("FIN compactacion")
