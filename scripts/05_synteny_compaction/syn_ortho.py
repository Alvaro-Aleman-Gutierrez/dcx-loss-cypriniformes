#!/usr/bin/env python3
"""Sintenia del bloque de DCX con ORTOLOGIA VERIFICADA por mejor hit reciproco (RBH).
Motivo: TMEM47 y PAK3 son multicopia; en v1 los anclajes cayeron en parálogos (demostrado en
Triplophysa_rosa, donde los anclajes dieron NC_079894.1 mientras el DCX real esta en NC_079891.1).
Esto invalidaba la afirmacion de que Paedocypris "conserva el bloque intacto".

Prueba RBH, sin tabla de simbolos:
  consulta q --(tblastn, genoma completo)--> locus L --> proteina del HSP
  proteina --(blastp)--> proteomas de referencia (Pimephales, Gobio, Danio: anotados y con DCX salvo Danio)
  mejor proteina de referencia R --(blastp)--> base de consultas
  se acepta el locus SOLO si vuelve a q.  Asi un exon de pak1/pak2 no puede pasar como PAK3loc.
"""
import os, subprocess, collections, re, sys
P = "${PROJECT}"
O = "${DATA}/syn2"; os.makedirs(O, exist_ok=True)
B2 = "${DATA}/p1b/genomes"
T = os.environ.get("THREADS", "8")
GEN = [g.strip() for g in os.environ.get("GENOMES", "").split(",") if g.strip()] or [
    "Danio_rerio", "Danionella_dracula", "Paedocypris_micromegethes", "Paedocypris_sp_Pulau_Singkep",
    "Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus",
    "Gobio_gobio", "Barbatula_barbatula", "Pimephales_promelas", "Rasbora_lateristriata",
    "Ctenopharyngodon_idella", "Triplophysa_rosa", "Danionella_cerebrum"]
SIN_DCX = {"Danio_rerio", "Danionella_dracula", "Paedocypris_micromegethes",
           "Paedocypris_sp_Pulau_Singkep", "Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus"}
def run(c, **k): return subprocess.run(c, check=True, capture_output=True, text=True, **k).stdout
def pick(*c):
    for x in c:
        if os.path.exists(x) or os.path.exists(x + ".ndb") or os.path.exists(x + ".nal") or os.path.exists(x + ".pdb"): return x
    return None
# --- bases: consultas (flancos + DCX) y proteomas de referencia
QF = f"{O}/queries_syn.faa"
if not os.path.exists(QF):
    with open(QF, "w") as fo:
        for f in (f"{P}/syn/block_q.faa", f"{P}/tbn/q_dcx.faa"):
            for l in open(f): fo.write(l)
# base DISCRIMINANTE del paso reciproco: flancos + DCX + TODA la familia DCX/DCLK.
# Imprescindible: sin los parálogos dclk/dcdc dentro, un locus de dclk1a reciproca a DCX y se
# acepta como DCX. Ese error hizo que la v1 "detectara" DCX en Danio rerio (RBH->DCLK1a XP_087125674.1).
DF = f"{O}/queries_disc.faa"
if not os.path.exists(DF) or os.path.getsize(DF) < 1000:
    with open(DF, "w") as fo:
        for f_ in (QF, f"{P}/census/queries.faa"):
            for l in open(f_): fo.write(l)
if not os.path.exists(f"{O}/qdb.pdb"): run(["makeblastdb", "-in", DF, "-dbtype", "prot", "-out", f"{O}/qdb"])
REFS = [f"{O}/prot/{s}.faa" for s in ("Pimephales", "Gobio", "Danrer")]
REFS = [r for r in REFS if os.path.exists(r) and os.path.getsize(r) > 1000]
if not REFS: sys.exit("faltan proteomas de referencia")
cat = f"{O}/refprot.faa"
if not os.path.exists(cat) or os.path.getmtime(cat) < max(os.path.getmtime(r) for r in REFS):
    with open(cat, "w") as fo:
        for r in REFS:
            for l in open(r): fo.write(l)
    run(["makeblastdb", "-in", cat, "-dbtype", "prot", "-out", f"{O}/refdb"])
elif not os.path.exists(f"{O}/refdb.pdb"):
    run(["makeblastdb", "-in", cat, "-dbtype", "prot", "-out", f"{O}/refdb"])
print("referencias:", [os.path.basename(r) for r in REFS], flush=True)
# --- paso 2 del RBH, con cache: proteina de referencia -> mejor consulta
cache = {}
cf = f"{O}/rbh_cache.tsv"
if os.path.exists(cf):
    for l in open(cf):
        a, b = l.rstrip("\n").split("\t"); cache[a] = b
def ref_a_consulta(ids):
    falta = [i for i in ids if i not in cache]
    if falta:
        sub = f"{O}/_sub.faa"; lst = f"{O}/_sub.ids"
        open(lst, "w").write("\n".join(falta) + "\n")
        open(sub, "w").write(run(["seqkit", "grep", "-f", lst, cat]))
        r = run(["blastp", "-query", sub, "-db", f"{O}/qdb", "-outfmt", "6 qseqid sseqid bitscore",
                 "-max_target_seqs", "1", "-max_hsps", "1", "-num_threads", T])
        got = {}
        for l in r.splitlines():
            a, b, bs = l.split("\t")
            if a not in got: got[a] = b
        with open(cf, "a") as fo:
            for i in falta:
                cache[i] = got.get(i, "NINGUNO"); fo.write(f"{i}\t{cache[i]}\n")
    return {i: cache[i] for i in ids}
out = open(f"{O}/syn_ortho_results.txt", "w")
def say(*a): print(*a, flush=True); out.write(" ".join(map(str, a)) + "\n"); out.flush()
def clean(s): return s.split("|")[1] if "|" in s and s.count("|") >= 2 else s.strip("|")
DCXQ = {l[1:].split()[0] for l in open(f"{P}/tbn/q_dcx.faa") if l.startswith(">")} | {"Q_DCX"}
for sp in GEN:
    db = pick(f"{P}/genomes/db_{sp}", f"{B2}/db_{sp}")
    if not db: say(f"\n== {sp}: sin base BLAST"); continue
    say(f"\n== {sp} [{'SIN DCX' if sp in SIN_DCX else 'con DCX'}]")
    tf = f"{O}/{sp}.tbn.tsv"
    if not os.path.exists(tf) or not os.path.getsize(tf):
        open(tf, "w").write(run(["tblastn", "-query", QF, "-db", db, "-evalue", "1e-8",
            "-outfmt", "6 qseqid sseqid pident length sstart send bitscore sseq", "-num_threads", T,
            "-max_target_seqs", "20", "-max_hsps", "5", "-seg", "no"]))
    # agrupar HSPs en loci (mismo contig, misma consulta, <100 kb)
    por = collections.defaultdict(list)
    for l in open(tf):
        f_ = l.rstrip("\n").split("\t")
        if len(f_) < 8: continue
        por[(f_[0], clean(f_[1]))].append((int(f_[4]), int(f_[5]), float(f_[6]), float(f_[2]), f_[7]))
    loci = []
    for (q, c), hs in por.items():
        hs.sort(key=lambda x: min(x[0], x[1]))
        grp = []
        for h in hs:
            if grp and min(h[0], h[1]) - max(grp[-1][0], grp[-1][1]) < 100000: grp.append(h)
            else:
                if grp: loci.append((q, c, grp))
                grp = [h]
        if grp: loci.append((q, c, grp))
    # proteina representativa por locus = HSP de mayor bitscore
    faa = f"{O}/{sp}.loci.faa"; meta = []
    with open(faa, "w") as fo:
        for i, (q, c, grp) in enumerate(loci):
            mejor = max(grp, key=lambda x: x[2])
            bs = sum(x[2] for x in grp); pos = min(min(x[0], x[1]) for x in grp)
            fo.write(f">L{i}\n{mejor[4].replace('-','').replace('*','X')}\n")
            meta.append(dict(i=i, q=q, c=c, pos=pos, bs=bs, pid=mejor[3], n=len(grp)))
    if not meta: say("   sin HSPs"); continue
    r = run(["blastp", "-query", faa, "-db", f"{O}/refdb", "-outfmt", "6 qseqid sseqid bitscore",
             "-max_target_seqs", "1", "-max_hsps", "1", "-num_threads", T])
    mejorref = {}
    for l in r.splitlines():
        a, b, bs = l.split("\t")
        if a not in mejorref: mejorref[a] = b
    vuelta = ref_a_consulta(sorted(set(mejorref.values())))
    ok = []
    for m in meta:
        R = mejorref.get(f"L{m['i']}")
        back = vuelta.get(R, "NINGUNO") if R else "NINGUNO"
        esDCX_q = m["q"] in DCXQ
        pasa = ((back == m["q"]) or (esDCX_q and back in DCXQ)) and m["bs"] >= 100
        m.update(R=R or "-", back=back, pasa=pasa)
        if pasa: ok.append(m)
    # un locus verificado por gen: el de mayor bitscore acumulado
    gen = {}
    for m in sorted(ok, key=lambda m: -m["bs"]):
        g = "DCX" if m["q"] in DCXQ else m["q"]
        if g not in gen: gen[g] = m
    for g, m in sorted(gen.items(), key=lambda kv: (kv[1]["c"], kv[1]["pos"])):
        say(f"   {g:34s} {m['c']:22s} pos={m['pos']:>12,}  HSPs={m['n']:2d} bit={m['bs']:7.0f} pid={m['pid']:.1f}%  RBH->{m['R']}")
    rech = [m for m in meta if not m["pasa"]]
    say(f"   -- loci candidatos={len(meta)}  verificados={len(ok)}  genes con ortologo unico={len(gen)}  rechazados por RBH={len(rech)}")
    for m in sorted(rech, key=lambda m: -m["bs"])[:6]:
        say(f"      RECHAZADO {m['q']:30s} {m['c']:20s} pos={m['pos']:>12,} bit={m['bs']:6.0f} -> RBH dio {m['back']}")
    # bloque: contigs y orden de los genes verificados
    cc = collections.Counter(m["c"] for m in gen.values())
    say(f"   >> reparto por contig: {dict(cc)}")
    if "DCX" in gen:
        d = gen["DCX"]; vec = sorted([m for m in gen.values() if m["c"] == d["c"] and m is not d], key=lambda m: m["pos"])
        say(f"   >> DCX en {d['c']}:{d['pos']:,}  vecinos en el mismo contig: "
            + ", ".join(f"{('DCX' if m['q'] in DCXQ else m['q'])}@{m['pos']:,}" for m in vec) if vec else "   >> DCX sin vecinos verificados en su contig")
    else:
        say("   >> DCX NO detectado (coherente con la perdida)" if sp in SIN_DCX else "   >> AVISO: DCX no detectado en una especie que deberia tenerlo")
out.write("\nFIN syn_ortho\n"); out.close(); print("FIN")
