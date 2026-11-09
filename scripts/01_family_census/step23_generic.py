#!/usr/bin/env python3
"""Pasos 2-3 del piloto DCX.
Paso 2: censo de la familia DCX (+ controles) en todos los genomas con miniprot.
Paso 3: extracción de CDS, alineamiento codónico, árbol de especies y RELAX (HyPhy)
        linajes sin DCX (Test) vs linajes con DCX (Reference).
Uso: python3 -I step23_pipeline.py   (rutas fijas abajo)
"""
import os, re, subprocess, json, collections, sys, glob
from Bio import SeqIO, Phylo
from Bio.Seq import Seq

P = "${PROJECT}"
G = f"{P}/genomes"; C = os.environ.get("C_DIR", f"{P}/census"); R = os.environ.get("R_DIR", f"{P}/relax")
QFILE = os.environ.get("QFILE", f"{C}/queries.faa")
MINIPROT = "${HOME}/anaconda3/envs/dcx_env/bin/miniprot"
HYPHY = "${HOME}/anaconda3/envs/hyphy/bin/hyphy"
THREADS = os.environ.get("THREADS", "6")
HYCPU = os.environ.get("HYCPU", THREADS)
FAST = os.environ.get("FAST") == "1"
os.makedirs(f"{C}/gff", exist_ok=True); os.makedirs(R, exist_ok=True)
log = open(os.environ.get("LOGF", f"{P}/step23.log"), "a")
def say(*a):
    print(*a, file=log, flush=True); print(*a, flush=True)
def run(cmd, **kw):
    say("$", cmd if isinstance(cmd, str) else " ".join(cmd))
    return subprocess.run(cmd, shell=isinstance(cmd, str), check=True, **kw)

species = []
for f in ["list.txt", "list2.txt", "list3.txt"]:
    for l in open(f"{G}/{f}"):
        if l.strip():
            acc, sp = l.split()[:2]
            if sp not in species: species.append(sp)
# usar solo un ensamblaje de D. cerebrum para RELAX/árbol (el v1); v2 solo para el censo
qlen = {r.id: len(r.seq) for r in SeqIO.parse(QFILE, "fasta")}

# ---------- Paso 2: miniprot ----------
for sp in species:
    out = f"{C}/gff/{sp}.gff"
    if os.path.exists(out) and os.path.getsize(out) > 0: continue
    if not os.path.exists(f"{G}/{sp}.fna"): say("FALTA genoma", sp); continue
    # 1) regiones candidatas con tBLASTn (ahorra memoria): ±150 kb alrededor de cada hit
    os.makedirs(f"{C}/regions", exist_ok=True)
    hits = subprocess.run(["tblastn", "-query", QFILE, "-db", f"{G}/db_{sp}", "-evalue", "1e-10",
                           "-outfmt", "6 sseqid sstart send", "-max_target_seqs", "5" if FAST else "50", "-num_threads", THREADS] + (["-task", "tblastn-fast"] if FAST else []),
                          capture_output=True, text=True).stdout
    iv = collections.defaultdict(list)
    for l in hits.splitlines():
        s, a, b = l.split(); s = s.split("|")[1] if "|" in s else s
        a, b = sorted((int(a), int(b))); iv[s].append((max(1, a - 150000), b + 150000))
    # extracción con seqkit faidx (no carga cromosomas enteros en memoria)
    reglist = f"{C}/regions/{sp}.list"
    with open(reglist, "w") as fo:
        for s, L in iv.items():
            L.sort(); merged = [list(L[0])]
            for a, b in L[1:]:
                if a <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], b)
                else: merged.append([a, b])
            for a, b in merged: fo.write(f"{s}:{a}-{b}\n")
    with open(f"{C}/regions/{sp}.fna", "w") as fo:
        run(["seqkit", "faidx", "-l", reglist, f"{G}/{sp}.fna"], stdout=fo, stderr=subprocess.DEVNULL)
    # cabeceras simples para miniprot (seq:ini-fin -> seq__ini)
    run(["sed", "-i", "-E", r"s/^>([^: ]+):([0-9]+)-[0-9]+.*/>\1__\2/", f"{C}/regions/{sp}.fna"])
    # 2) miniprot sobre las regiones
    with open(out + ".tmp", "w") as fo:
        run([MINIPROT, "-t", THREADS, "--gff", "-P", sp, "-N", "30", f"{C}/regions/{sp}.fna", QFILE], stdout=fo, stderr=subprocess.DEVNULL)
    os.rename(out + ".tmp", out)

def parse_gff(fn):
    """devuelve lista de mRNA con atributos y sus CDS"""
    mr = {}; order = []
    for l in open(fn):
        if l.startswith("#") or "\t" not in l: continue
        f = l.rstrip("\n").split("\t")
        attr = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "mRNA":
            q = attr["Target"].split()[0]; qs, qe = map(int, attr["Target"].split()[1:3])
            m = dict(id=attr["ID"], seq=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                     score=float(f[5]), q=q, qs=qs, qe=qe, ident=float(attr.get("Identity", 0)),
                     fs=int(attr.get("Frameshift", 0)), stop=int(attr.get("StopCodon", 0)), cds=[])
            mr[m["id"]] = m; order.append(m["id"])
        elif f[2] == "CDS":
            mr[attr["Parent"]]["cds"].append((int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0))
    return [mr[i] for i in order]

census = {}; chosen = {}
for sp in species:
    fn = f"{C}/gff/{sp}.gff"
    if not os.path.exists(fn): continue
    ms = parse_gff(fn)
    # agrupar en loci por solapamiento
    ms.sort(key=lambda m: (m["seq"], m["start"]))
    loci = []
    for m in ms:
        if loci and loci[-1][0]["seq"] == m["seq"] and m["start"] <= max(x["end"] for x in loci[-1]):
            loci[-1].append(m)
        else: loci.append([m])
    best_per_q = {}
    for L in loci:
        top = max(L, key=lambda m: m["score"])  # la consulta que mejor explica el locus
        cov = (top["qe"] - top["qs"] + 1) / qlen[top["q"]]
        rec = dict(top, cov=cov)
        q = top["q"]
        if q not in best_per_q or rec["score"] > best_per_q[q]["score"]: best_per_q[q] = rec
    census[sp] = {}
    for q in qlen:
        b = best_per_q.get(q)
        if b and b["cov"] >= 0.5:
            status = "presente" if (b["fs"] == 0 and b["stop"] == 0) else "pseudo?"
            census[sp][q] = (status, round(b["cov"], 2), round(b["ident"], 3), b["fs"], b["stop"])
            chosen[(sp, q)] = b
        elif b:
            census[sp][q] = ("parcial", round(b["cov"], 2), round(b["ident"], 3), b["fs"], b["stop"])
        else:
            census[sp][q] = ("ausente", 0, 0, 0, 0)

qs = list(qlen)
with open(f"{C}/census_table.tsv", "w") as fo:
    fo.write("especie\t" + "\t".join(q[2:] for q in qs) + "\n")
    for sp in census:
        fo.write(sp + "\t" + "\t".join(f"{census[sp][q][0]}({census[sp][q][1]};{census[sp][q][2]})" for q in qs) + "\n")
say("censo escrito:", f"{C}/census_table.tsv")

# ---------- Paso 3: CDS ----------
_IDX = {}
def extract_cds(sp, m):
    if sp not in _IDX: _IDX[sp] = SeqIO.index(f"{C}/regions/{sp}.fna", "fasta")
    rec = _IDX[sp][m["seq"]]
    parts = sorted(m["cds"])
    s = "".join(str(rec.seq[a - 1:b]) for a, b, ph in parts)
    if m["strand"] == "-": s = str(Seq(s).reverse_complement())
    s = s[: len(s) - len(s) % 3]
    return s

relax_species = [s for s in species if s != "Danionella_cerebrum_v2"]
genes = [q for q in qs if q != "Q_DCX"]
os.makedirs(f"{R}/cds", exist_ok=True)
for q in qs:
    fo = f"{R}/cds/{q[2:]}.fna"
    if os.path.exists(fo): continue
    with open(fo, "w") as out:
        for sp in relax_species:
            m = chosen.get((sp, q))
            if not m or census[sp][q][0] != "presente": continue
            s = extract_cds(sp, m)
            prot = str(Seq(s).translate())
            if prot.rstrip("*").count("*") > 0: continue
            out.write(f">{sp}\n{s}\n")

def codon_align(q):
    fn = f"{R}/cds/{q}.fna"; recs = list(SeqIO.parse(fn, "fasta"))
    if len(recs) < 8: return None
    with open(f"{R}/cds/{q}.faa", "w") as fa:
        for r in recs:
            p = str(r.seq.translate()).rstrip("*"); fa.write(f">{r.id}\n{p}\n")
    with open(f"{R}/cds/{q}.aln.faa", "w") as fo:
        run(["mafft", "--auto", "--quiet", f"{R}/cds/{q}.faa"], stdout=fo)
    nt = {r.id: str(r.seq) for r in recs}
    out = {}
    for r in SeqIO.parse(f"{R}/cds/{q}.aln.faa", "fasta"):
        s = nt[r.id]; i = 0; cod = []
        for aa in str(r.seq):
            if aa == "-": cod.append("---")
            else: cod.append(s[i:i + 3]); i += 3
        out[r.id] = "".join(cod)
    # quitar columnas codónicas con >50% gaps
    L = len(next(iter(out.values()))) // 3
    keep = [k for k in range(L) if sum(v[3*k:3*k+3] == "---" for v in out.values()) <= 0.5 * len(out)]
    with open(f"{R}/{q}.codon.fas", "w") as fo:
        for k, v in out.items(): fo.write(f">{k}\n" + "".join(v[3*j:3*j+3] for j in keep) + "\n")
    return f"{R}/{q}.codon.fas"

alns = {}
for q in qs:
    a = codon_align(q[2:])
    if a: alns[q[2:]] = a
say("alineamientos:", list(alns))

# ---------- árbol de especies: concatenado de proteínas ----------
sup = f"{R}/species_concat.faa"
if not os.environ.get("SPTREE") and not os.path.exists(f"{R}/species.treefile"):
    cat = collections.defaultdict(str)
    for q in alns:
        seqs = {r.id: str(r.seq.translate()) for r in SeqIO.parse(alns[q], "fasta")}
        L = len(next(iter(seqs.values())))
        for sp in relax_species: cat[sp] += seqs.get(sp, "-" * L)
    with open(sup, "w") as fo:
        for sp, s in cat.items():
            if s.replace("-", "").replace("X", ""): fo.write(f">{sp}\n{s}\n")
    run(["iqtree2", "-s", sup, "-m", "LG+F+R4", "-B", "1000", "-T", "6", "-o", "Astyanax_mexicanus",
         "--prefix", f"{R}/species", "-redo", "--quiet"])
if os.environ.get("SPTREE"):
    import shutil; shutil.copyfile(os.environ["SPTREE"], f"{R}/species.treefile")
sptree = Phylo.read(f"{R}/species.treefile", "newick")
sptree.root_with_outgroup("Astyanax_mexicanus")
Phylo.write(sptree, f"{R}/species_rooted.nwk", "newick")

dcx_absent = set(open(os.environ["DCX_ABSENT"]).read().split()) if os.environ.get("DCX_ABSENT") else {sp for sp in relax_species if census.get(sp, {}).get("Q_DCX", ("ausente",))[0] not in ("presente", "pseudo?")}
say("especies sin DCX:", sorted(dcx_absent))

def label_tree(taxa, test_taxa):
    t = Phylo.read(f"{R}/species_rooted.nwk", "newick")
    for tip in [x.name for x in t.get_terminals()]:
        if tip not in taxa: t.prune(tip)
    for cl in t.find_clades():
        tips = {x.name for x in cl.get_terminals()}
        cl._lab = "Test" if tips and tips <= test_taxa else "Reference"
        cl.confidence = None
    def nw(c):
        if c.is_terminal(): return f"{c.name}{{{c._lab}}}:{c.branch_length or 0.01}"
        return "(" + ",".join(nw(x) for x in c.clades) + f"){{{c._lab}}}:{c.branch_length or 0.01}"
    root = t.root
    return "(" + ",".join(nw(x) for x in root.clades) + ");"

results = []
for q, a in alns.items():
    taxa = [r.id for r in SeqIO.parse(a, "fasta")]
    for setname, test in [("todas_las_perdidas", dcx_absent), ("solo_Danio", {s for s in dcx_absent if s.startswith("Danio_")})]:
        use = taxa if setname == "todas_las_perdidas" else [t for t in taxa if t in test or t not in dcx_absent]
        tset = set(use) & test
        if len(tset) < 2 or len(set(use) - tset) < 3: continue
        tag = f"{q}.{setname}"
        tre = f"{R}/{tag}.nwk"; open(tre, "w").write(label_tree(set(use), tset))
        alnf = f"{R}/{tag}.fas"
        with open(alnf, "w") as fo:
            for r in SeqIO.parse(a, "fasta"):
                if r.id in use: fo.write(f">{r.id}\n{r.seq}\n")
        js = f"{R}/{tag}.RELAX.json"
        if not os.path.exists(js):
            try:
                run([HYPHY, "relax", "--alignment", alnf, "--tree", tre, "--test", "Test", "--reference", "Reference",
                     "--models", "Minimal", "--output", js, f"CPU={HYCPU}"], stdout=open(f"{R}/{tag}.log", "w"), stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                say("RELAX falló", tag); continue
        d = json.load(open(js))
        tr = d.get("test results", {})
        results.append((q, setname, len(use), len(tset), tr.get("relaxation or intensification parameter"), tr.get("p-value")))
        say("RELAX", results[-1])

with open(f"{R}/relax_summary.tsv", "w") as fo:
    fo.write("gen\tconjunto\tn_taxa\tn_test\tK\tp\n")
    for r in results: fo.write("\t".join(map(str, r)) + "\n")
say("PASOS 2-3 COMPLETOS")
