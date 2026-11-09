#!/usr/bin/env python3
"""Busqueda de restos de pseudogen de DCX en el intervalo sintenico de las especies sin DCX.
Estrategia: 1) ubicar los genes flanqueantes (TMEM47, PAK3loc, ZFPL1, RWDD2B) por tblastn;
2) definir el intervalo del bloque en el cromosoma donde co-ocurren; 3) tblastn PERMISIVO de DCX
solo en ese intervalo; 4) clasificar cada HSP por mejor hit reciproco contra la familia.
Controles: especies con DCX, donde debe recuperarse el gen intacto."""
import os, subprocess, collections, sys
P = "${PROJECT}"
D = "${DATA}/pseudo"; os.makedirs(D, exist_ok=True)
B2 = "${DATA}/p1b/genomes"
T = os.environ.get("THREADS", "4")
FLANK = {"TMEM47", "PAK3loc_Pimephales_LOC120471754", "ZFPL1", "RWDD2B"}
PAD = 150000
AUSENTES = [l.strip() for l in open(f"{P}/relax/dcx_absent.txt") if l.strip()] + \
           ["Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus"]
CONTROLES = ["Gobio_gobio", "Barbatula_barbatula", "Triplophysa_rosa", "Danionella_cerebrum",
             "Pimephales_promelas", "Rasbora_lateristriata"]
def db_of(sp):
    for c in (f"{P}/genomes/db_{sp}", f"{B2}/db_{sp}"):
        if os.path.exists(c + ".ndb") or os.path.exists(c + ".nal"): return c
    return None
def fna_of(sp):
    for c in (f"{P}/genomes/{sp}.fna", f"{B2}/{sp}.fna"):
        if os.path.exists(c): return c
    return None
def run(cmd, **kw): return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout
out = open(f"{D}/pseudo_results.txt", "w")
def say(*a):
    print(*a, flush=True); out.write(" ".join(map(str, a)) + "\n"); out.flush()
for sp in AUSENTES + CONTROLES:
    db, fna = db_of(sp), fna_of(sp)
    if not db or not fna: say(f"== {sp}: SIN genoma/base"); continue
    grupo = "AUSENTE" if sp in AUSENTES else "CONTROL(con DCX)"
    # 1) flancos
    f = f"{D}/{sp}.flank.tsv"
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        r = run(["tblastn", "-query", f"{P}/syn/block_q.faa", "-db", db, "-evalue", "1e-10",
                 "-outfmt", "6 qseqid sseqid pident length sstart send bitscore",
                 "-num_threads", T, "-max_target_seqs", "5", "-max_hsps", "3"])
        open(f, "w").write(r)
    best = {}
    for l in open(f):
        q, s, pid, ln, ss, se, bs = l.split()
        if q not in FLANK: continue
        bs = float(bs)
        if q not in best or bs > best[q][0]: best[q] = (bs, s, int(ss))
    porcrom = collections.defaultdict(list)
    for q, (bs, s, ss) in best.items(): porcrom[s].append((q, ss))
    if not porcrom: say(f"== {sp} [{grupo}]: sin flancos"); continue
    crom, genes = max(porcrom.items(), key=lambda kv: len(kv[1]))
    pos = [p for _, p in genes]
    ini, fin = max(1, min(pos) - PAD), max(pos) + PAD
    say(f"== {sp} [{grupo}]  bloque en {crom}  flancos={sorted(g for g,_ in genes)}  intervalo={ini}-{fin} ({(fin-ini)/1000:.0f} kb)")
    if len(genes) < 2:
        say("   AVISO: solo 1 flanco en este cromosoma; el bloque esta disperso, intervalo poco informativo")
    # 2) extraer intervalo
    reg = f"{D}/{sp}.region.fna"
    if not os.path.exists(reg) or os.path.getsize(reg) == 0:
        crom_id = crom.replace("ref|", "").replace("|", "")
        open(f"{D}/{sp}.list", "w").write(f"{crom_id}:{ini}-{fin}\n")
        r = run(["seqkit", "faidx", "-l", f"{D}/{sp}.list", fna])
        open(reg, "w").write(r)
    if os.path.getsize(reg) < 1000: say("   fallo al extraer la region"); continue
    # 3) tblastn permisivo de DCX
    h = f"{D}/{sp}.dcxperm.tsv"
    if not os.path.exists(h) or os.path.getsize(h) == 0:
        r = run(["tblastn", "-query", f"{P}/tbn/q_dcx.faa", "-subject", reg, "-evalue", "10",
                 "-seg", "no", "-comp_based_stats", "0", "-max_hsps", "20",
                 "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq"])
        open(h, "w").write(r)
    hsps = [l.split("\t") for l in open(h) if l.strip()]
    hsps = [x for x in hsps if int(x[3]) >= 25]
    if not hsps: say("   sin HSPs de DCX (>=25 aa) en el intervalo  -> NINGUN RESTO"); continue
    # 4) RBH contra la familia
    faa = f"{D}/{sp}.hsp.faa"
    with open(faa, "w") as fo:
        for i, x in enumerate(hsps):
            s = x[10].replace("-", "").replace("*", "X")
            fo.write(f">h{i}_{x[0]}_q{x[4]}-{x[5]}_s{x[6]}-{x[7]}_pid{x[2]}_E{x[8]}\n{s}\n")
    rbh = run(["blastp", "-query", faa, "-db", f"{P}/tbn/famdb", "-outfmt", "6 qseqid sseqid pident bitscore",
               "-max_target_seqs", "1", "-max_hsps", "1"])
    bestrbh = {}
    for l in rbh.splitlines():
        q, s, pid, bs = l.split("\t")
        if q not in bestrbh: bestrbh[q] = (s, pid)
    ndcx = 0
    for i, x in enumerate(hsps):
        k = [kk for kk in bestrbh if kk.startswith(f"h{i}_")]
        tgt, pid = bestrbh[k[0]] if k else ("sin_hit", "-")
        esdcx = "_dcx_" in tgt.lower() or "_DCX_" in tgt
        if esdcx: ndcx += 1
        say(f"   HSP{i}: q{x[4]}-{x[5]} len={x[3]} pid={x[2]}% E={x[8]} pos={x[6]}-{x[7]} -> RBH={tgt} ({pid}%) {'<< DCX' if esdcx else ''}")
    say(f"   resumen {sp}: {len(hsps)} HSPs, {ndcx} con mejor hit = DCX")
out.write("FIN pseudo\n"); out.close()
print("FIN")
