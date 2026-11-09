#!/usr/bin/env python3
"""Restos de pseudogen de DCX, version 2 (corregida).
Calibrado en los controles con DCX: el orden real del bloque es
   TMEM47 -- DCX -- PAK3loc -- ZFPL1 -- RWDD2B
(Barbatula: 5.961 | 6.024 | 6.113 | 6.138 | 6.139 Mb;  Triplophysa: 8.485 | 8.533 | 8.727 | 8.754 | 8.757 Mb)
=> el intervalo correcto es el que queda ENTRE TMEM47 y PAK3loc, con margen de 50 kb.
Correcciones frente a v1: (a) parseo correcto del id (ref|X| -> X), (b) recorte a la longitud real
del contig leida del .fai, (c) el intervalo ya no excluye el DCX real, (d) se exige que los dos
anclajes esten en el MISMO contig; si no, se informa 'bloque roto / contexto perdido'."""
import os, subprocess, collections, re
P = "${PROJECT}"
D = "${DATA}/pseudo2"; os.makedirs(D, exist_ok=True)
B2 = "${DATA}/p1b/genomes"
T = os.environ.get("THREADS", "4"); PAD = 50000
A1, A2 = "TMEM47", "PAK3loc_Pimephales_LOC120471754"
AUS = [l.strip() for l in open(f"{P}/relax/dcx_absent.txt") if l.strip()] + \
      ["Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus"]
CTL = ["Gobio_gobio", "Barbatula_barbatula", "Triplophysa_rosa", "Danionella_cerebrum",
       "Pimephales_promelas", "Rasbora_lateristriata", "Ctenopharyngodon_idella", "Labeo_rohita"]
def pick(sp, *cands):
    for c in cands:
        if os.path.exists(c) or os.path.exists(c + ".ndb") or os.path.exists(c + ".nal"): return c
    return None
def run(c, **k): return subprocess.run(c, check=True, capture_output=True, text=True, **k).stdout
def clean(s): return s.split("|")[1] if "|" in s else s
out = open(f"{D}/pseudo2_results.txt", "w")
def say(*a): print(*a, flush=True); out.write(" ".join(map(str, a)) + "\n"); out.flush()
def lens(fna):
    fai = fna + ".fai"
    if not os.path.exists(fai): run(["seqkit", "faidx", fna])
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
for sp in AUS + CTL:
    grp = "AUSENTE" if sp in AUS else "CONTROL(con DCX)"
    db = pick(sp, f"{P}/genomes/db_{sp}", f"{B2}/db_{sp}")
    fna = pick(sp, f"{P}/genomes/{sp}.fna", f"{B2}/{sp}.fna")
    if not db or not fna: say(f"== {sp} [{grp}]: SIN genoma/base"); continue
    f = f"{D}/{sp}.flank.tsv"
    if not os.path.exists(f) or not os.path.getsize(f):
        open(f, "w").write(run(["tblastn", "-query", f"{P}/syn/block_q.faa", "-db", db, "-evalue", "1e-10",
            "-outfmt", "6 qseqid sseqid pident length sstart send bitscore", "-num_threads", T,
            "-max_target_seqs", "10", "-max_hsps", "5"]))
    # mejor locus por gen ancla, en cada contig
    hits = collections.defaultdict(dict)   # contig -> gen -> (bitscore, pos)
    for l in open(f):
        q, s, pid, ln, ss, se, bs = l.split("\t"); bs = float(bs); s = clean(s)
        if q not in (A1, A2): continue
        cur = hits[s].get(q)
        if not cur or bs > cur[0]: hits[s][q] = (bs, int(ss))
    juntos = {c: g for c, g in hits.items() if A1 in g and A2 in g}
    if not juntos:
        c1 = [c for c, g in hits.items() if A1 in g]; c2 = [c for c, g in hits.items() if A2 in g]
        say(f"== {sp} [{grp}]: BLOQUE ROTO -> TMEM47 en {c1[:2]}, PAK3loc en {c2[:2]}; "
            f"contexto perdido, no se puede delimitar el intervalo de DCX")
        continue
    contig = max(juntos, key=lambda c: juntos[c][A1][0] + juntos[c][A2][0])
    p1, p2 = juntos[contig][A1][1], juntos[contig][A2][1]
    L = lens(fna).get(contig)
    ini, fin = max(1, min(p1, p2) - PAD), min(L, max(p1, p2) + PAD) if L else max(p1, p2) + PAD
    say(f"== {sp} [{grp}]  {contig} (len={L})  TMEM47={p1}  PAK3loc={p2}  intervalo={ini}-{fin} ({(fin-ini)/1000:.0f} kb)")
    reg = f"{D}/{sp}.region.fna"
    if not os.path.exists(reg) or os.path.getsize(reg) < 1000:
        open(reg, "w").write(run(["seqkit", "faidx", fna, f"{contig}:{ini}-{fin}"]))
    if os.path.getsize(reg) < 1000: say("   fallo al extraer la region"); continue
    h = f"{D}/{sp}.dcxperm.tsv"
    if not os.path.exists(h) or not os.path.getsize(h):
        open(h, "w").write(run(["tblastn", "-query", f"{P}/tbn/q_dcx.faa", "-subject", reg, "-evalue", "10",
            "-seg", "no", "-comp_based_stats", "0", "-max_hsps", "20",
            "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq"]))
    hs = [l.split("\t") for l in open(h) if l.strip()]
    hs = [x for x in hs if int(x[3]) >= 25]
    if not hs: say("   sin HSPs de DCX (>=25 aa) -> NINGUN RESTO DETECTABLE"); continue
    faa = f"{D}/{sp}.hsp.faa"
    with open(faa, "w") as fo:
        for i, x in enumerate(hs):
            fo.write(f">h{i}\n{x[10].replace('-','').replace('*','X')}\n")
    rbh = {}
    for l in run(["blastp", "-query", faa, "-db", f"{P}/tbn/famdb", "-outfmt",
                  "6 qseqid sseqid pident bitscore", "-max_target_seqs", "1", "-max_hsps", "1"]).splitlines():
        q, s, pid, bs = l.split("\t")
        if q not in rbh: rbh[q] = (s, pid)
    ndcx = fuertes = 0
    for i, x in enumerate(hs):
        tgt, pid = rbh.get(f"h{i}", ("sin_hit", "-"))
        es = bool(re.search(r"_dcx_|_DCX_|LOC113084161|LOC103035707|LOC102692282", tgt, re.I))
        E = float(x[8]); fuerte = E < 1e-3
        if es: ndcx += 1
        if es and fuerte: fuertes += 1
        if es or fuerte:
            say(f"   HSP{i}: q{x[4]}-{x[5]} len={x[3]} pid={x[2]}% E={x[8]} pos={x[6]}-{x[7]} -> {tgt} ({pid}%)"
                f" {'<<< DCX' if es else ''}{' [SIGNIFICATIVO]' if fuerte else ' [ruido, E>1e-3]'}")
    say(f"   resumen {sp}: {len(hs)} HSPs; {ndcx} con RBH=DCX; {fuertes} de ellos con E<1e-3 (significativos)")
out.write("FIN pseudo2\n"); out.close(); print("FIN")
