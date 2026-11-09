#!/usr/bin/env python3
"""Extiende la matriz de 48/50 genes (punto 5) a Cobitidae (Misgurnus, Paramisgurnus) e
Ictalurus furcatus, para que la Figura 1 pueda mostrar las CUATRO perdidas de DCX.
Reutiliza los CDS ya calculados en p5 (p5/relax/cds/*.fna) y solo añade los taxones nuevos:
  tblastn de las 50 consultas -> regiones +-150 kb -> miniprot -> mejor locus por gen -> CDS
  -> se anexa a los CDS de p5 -> MAFFT -> retrotraduccion -> supermatriz -> IQ-TREE
No toca los datos de p5 ni vuelve a ejecutar RELAX.
"""
import os, subprocess, collections, glob
from Bio import SeqIO
from Bio.Seq import Seq
P = "${PROJECT}"
P5 = "${DATA}/p5"
O = "${DATA}/matriz_ext"; os.makedirs(f"{O}/regions", exist_ok=True)
GB = "${DATA}/p1b/genomes"
MINIPROT = "${HOME}/anaconda3/envs/dcx_env/bin/miniprot"
T = os.environ.get("THREADS", "6")
QF = f"{P5}/random_queries.faa"
NUEVOS = ["Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus", "Ictalurus_furcatus"]
PAD = 150000
def run(c, **k): return subprocess.run(c, check=True, capture_output=True, text=True, **k).stdout
qlen = {r.id: len(r.seq) for r in SeqIO.parse(QF, "fasta")}
# ---------- 1) regiones candidatas + miniprot ----------
for sp in NUEVOS:
    gff = f"{O}/{sp}.gff"
    if os.path.exists(gff) and os.path.getsize(gff): continue
    reg = f"{O}/regions/{sp}.fna"
    if not (os.path.exists(reg) and os.path.getsize(reg) > 1000):
        hits = run(["tblastn", "-query", QF, "-db", f"{GB}/db_{sp}", "-evalue", "1e-10", "-task", "tblastn-fast",
                    "-outfmt", "6 sseqid sstart send", "-max_target_seqs", "5", "-num_threads", T])
        iv = collections.defaultdict(list)
        for l in hits.splitlines():
            s, a, b = l.split("\t"); a, b = int(a), int(b)
            iv[s.split("|")[1] if s.count("|") >= 2 else s].append((min(a, b) - PAD, max(a, b) + PAD))
        lst = []
        for c, vs in iv.items():
            vs.sort(); cur = list(vs[0])
            for a, b in vs[1:]:
                if a <= cur[1]: cur[1] = max(cur[1], b)
                else: lst.append((c, cur)); cur = [a, b]
            lst.append((c, cur))
        with open(f"{O}/{sp}.list", "w") as fo:
            for c, (a, b) in lst: fo.write(f"{c}:{max(1,a)}-{b}\n")
        open(reg, "w").write(run(["seqkit", "faidx", "-l", f"{O}/{sp}.list", f"{GB}/{sp}.fna"]))
        run(["sed", "-i", "-E", r"s/^>([^: ]+):([0-9]+)-[0-9]+.*/>\1__\2/", reg])
    print(f"{sp}: region {os.path.getsize(reg)/1e6:.0f} Mb -> miniprot", flush=True)
    open(gff, "w").write(run([MINIPROT, "-t", T, "--gff", "-P", sp, "-N", "30", reg, QF]))
# ---------- 2) parsear GFF y elegir el mejor locus por gen ----------
def parse_gff(fn):
    mr = {}
    for l in open(fn):
        if l.startswith("#") or not l.strip(): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            mr[attr["ID"]] = dict(seq=f[0], strand=f[6], q=attr.get("Target", " ").split()[0],
                                  ident=float(attr.get("Identity", 0)),
                                  qs=int(attr.get("Target", "x 1 1").split()[1]),
                                  qe=int(attr.get("Target", "x 1 1").split()[2]),
                                  score=int(attr.get("Rank", 99)), cds=[])
        elif f[2] == "CDS" and attr.get("Parent") in mr:
            mr[attr["Parent"]]["cds"].append((int(f[3]), int(f[4])))
    return list(mr.values())
_IDX = {}
def extract_cds(sp, m):
    if sp not in _IDX: _IDX[sp] = SeqIO.index(f"{O}/regions/{sp}.fna", "fasta")
    rec = _IDX[sp][m["seq"]]
    parts = sorted(m["cds"])
    s = "".join(str(rec.seq[a - 1:b]) for a, b in parts)
    if m["strand"] == "-": s = str(Seq(s).reverse_complement())
    return s
nuevos_cds = collections.defaultdict(dict)   # gen -> sp -> secuencia
for sp in NUEVOS:
    ms = parse_gff(f"{O}/{sp}.gff")
    mejor = {}
    for m in ms:
        cov = (m["qe"] - m["qs"] + 1) / qlen.get(m["q"], 1)
        m["cov"] = cov
        k = m["q"]
        if cov < 0.5: continue
        if k not in mejor or m["ident"] * cov > mejor[k]["ident"] * mejor[k]["cov"]: mejor[k] = m
    n = 0
    for q, m in mejor.items():
        s = extract_cds(sp, m)
        if len(s) % 3: s = s[:len(s) - len(s) % 3]
        prot = str(Seq(s).translate())
        if prot.rstrip("*").count("*") > 0 or len(s) < 150: continue
        nuevos_cds[q][sp] = s; n += 1
    print(f"{sp}: {n} genes con CDS utilizable de {len(qlen)}", flush=True)
# ---------- 3) anexar a los CDS de p5, realinear ----------
os.makedirs(f"{O}/cds", exist_ok=True); os.makedirs(f"{O}/aln", exist_ok=True)
hechos = 0
for fn in sorted(glob.glob(f"{P5}/relax/cds/*.fna")):
    g = os.path.basename(fn)[:-4]
    q = [k for k in qlen if k[2:] == g or k == g]
    q = q[0] if q else None
    recs = list(SeqIO.parse(fn, "fasta"))
    extra = nuevos_cds.get(q, {}) if q else {}
    out = f"{O}/cds/{g}.fna"
    with open(out, "w") as fo:
        for r in recs: fo.write(f">{r.id}\n{r.seq}\n")
        for sp, s in extra.items(): fo.write(f">{sp}\n{s}\n")
    # MAFFT sobre proteinas + retrotraduccion
    faa = f"{O}/cds/{g}.faa"
    with open(faa, "w") as fa:
        for r in SeqIO.parse(out, "fasta"):
            fa.write(f">{r.id}\n{Seq(str(r.seq)).translate()}\n")
    aln = f"{O}/aln/{g}.aln.faa"
    with open(aln, "w") as fo: fo.write(run(["mafft", "--auto", "--quiet", faa]))
    nt = {r.id: str(r.seq) for r in SeqIO.parse(out, "fasta")}
    cod = {}
    for r in SeqIO.parse(aln, "fasta"):
        s = nt[r.id]; i = 0; buf = []
        for c in str(r.seq):
            if c == "-": buf.append("---")
            else: buf.append(s[i*3:i*3+3]); i += 1
        cod[r.id] = "".join(buf)
    L = len(next(iter(cod.values())))
    keep = [j for j in range(0, L, 3)
            if sum(cod[k][j:j+3] == "---" for k in cod) <= 0.5 * len(cod)]
    with open(f"{O}/aln/{g}.codon.fas", "w") as fo:
        for k, v in cod.items(): fo.write(f">{k}\n" + "".join(v[j:j+3] for j in keep) + "\n")
    hechos += 1
print(f"alineamientos reconstruidos: {hechos}", flush=True)
