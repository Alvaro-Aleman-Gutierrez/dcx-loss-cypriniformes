#!/usr/bin/env python3
"""Exon terminal propio de la isoforma DCL de dclk1a: cola codificante (11 aa + stop = 36 nt) + 600 nt de 3'UTR,
tomada del ARNm DCL de una especie de referencia (coordenadas CDS del registro XM). BLASTN (dc-megablast) en cada
genoma, aceptando solo hits a <= 300 kb del locus dclk1a (definido por tBLASTn de los aa 1-335 de DCL de
Ctenopharyngodon, compartidos con la isoforma larga). Se traduce la cola en el genoma para ver si el ORF esta intacto."""
import subprocess, sys, re, urllib.request, urllib.parse, time
from Bio import SeqIO
from Bio.Seq import Seq
D = "${DATA}"; O = f"{D}/isoformas/cola"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u, **p):
    for i in range(6):
        try: return urllib.request.urlopen(E + u, data=urllib.parse.urlencode(p).encode(), timeout=90).read().decode()
        except Exception: time.sleep(3 + 3 * i)
XP = sys.argv[1]
gp = get("efetch.fcgi", db="protein", id=XP, rettype="gp", retmode="text")
m = re.search(r'coded_by="(complement\()?([^:]+):<?(\d+)\.\.>?(\d+)', gp); acc, a, b = m.group(2), int(m.group(3)), int(m.group(4))
s = "".join(get("efetch.fcgi", db="nuccore", id=acc, rettype="fasta", retmode="text").split("\n")[1:]).upper()
q = s[b - 36:b + 600]; open(f"{O}/q_nt.fna", "w").write(f">{XP}_DCLexon\n{q}\n")
print(f"consulta {XP} ({acc}, CDS {a}-{b}); cola codificada: {Seq(q[:36]).translate()}  ({len(q)} nt)")
full = next(str(r.seq) for r in SeqIO.parse(f"{D}/isoformas/dclk_isoformas.faa", "fasta") if "XP_051766222.1" in r.id)
open(f"{O}/q_aa.faa", "w").write(f">dcl_1_335\n{full[:335]}\n")
for sp in sys.argv[2:]:
    db = f"{D}/genomes/db_{sp}"
    ta = [l.split("\t") for l in subprocess.run(["tblastn", "-query", f"{O}/q_aa.faa", "-db", db, "-evalue", "1e-20",
          "-outfmt", "6 sseqid sstart send bitscore", "-max_target_seqs", "5"], capture_output=True, text=True).stdout.splitlines()]
    if not ta: print(f"{sp:30s} sin locus dclk1"); continue
    # locus por puntuacion acumulada (dclk1a y dclk1b compiten; se listan los dos mejores)
    loci = {}
    for h in ta: loci.setdefault(h[0], []).append(h)
    nb = subprocess.run(["blastn", "-task", "dc-megablast", "-query", f"{O}/q_nt.fna", "-db", db, "-evalue", "1e-5",
                         "-outfmt", "6 sseqid qstart qend sstart send pident bitscore", "-max_target_seqs", "50"],
                        capture_output=True, text=True).stdout.splitlines()
    nb = [l.split("\t") for l in nb]
    res = []
    for sid, hs in sorted(loci.items(), key=lambda x: -sum(float(h[3]) for h in x[1]))[:2]:
        c0 = min(min(int(h[1]), int(h[2])) for h in hs); c1 = max(max(int(h[1]), int(h[2])) for h in hs)
        cerca = [h for h in nb if h[0] == sid and c0 - 3e5 <= int(h[3]) <= c1 + 3e5]
        cov = set()
        for h in cerca: cov.update(range(int(h[1]), int(h[2]) + 1))
        cola = [h for h in cerca if int(h[1]) <= 4 and int(h[2]) >= 36]
        trad = "-"
        if cola:
            h = cola[0]; x, y = int(h[3]), int(h[4]); qs = int(h[1]); rev = y < x
            ini = x - (qs - 1) if not rev else x + (qs - 1)
            rng = f"{ini}-{ini+35}" if not rev else f"{ini-35}-{ini}"
            nt = "".join(subprocess.run(["blastdbcmd", "-db", db, "-entry", sid, "-range", rng], capture_output=True,
                                        text=True).stdout.split("\n")[1:]).upper()
            if rev: nt = str(Seq(nt).reverse_complement())
            trad = str(Seq(nt[:36]).translate())
        pid = max((float(h[5]) for h in cerca), default=0)
        res.append(f"{sid.split('|')[-2] if '|' in sid else sid}:{c0}  exonDCL cobertura={len(cov)}/{len(q)} id_max={pid:.0f}% cola={trad}")
    print(f"{sp:30s} " + "  ||  ".join(res), flush=True)
