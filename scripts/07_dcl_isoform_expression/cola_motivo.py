#!/usr/bin/env python3
"""Busqueda del motivo de la cola DCL ([EK]-x-D-K-x-x-[IV]-D-x-[ST]-A-stop y variantes laxas) en las 6 fases,
dentro de +-200 kb de cada locus dclk1 (tBLASTn de los aa 1-335). Complementa a exon_dcl.py en genomas divergentes."""
import subprocess, sys, re
from Bio.Seq import Seq
D = "${DATA}"; O = f"{D}/isoformas/cola"
MOT = re.compile(r"E.D[KR]..[IVL]D.[ST]A\*")
for sp in sys.argv[1:]:
    db = f"{D}/genomes/db_{sp}"
    ta = [l.split("\t") for l in subprocess.run(["tblastn", "-query", f"{O}/q_aa.faa", "-db", db, "-evalue", "1e-20",
          "-outfmt", "6 sseqid sstart send bitscore slen", "-max_target_seqs", "5"], capture_output=True, text=True).stdout.splitlines()]
    loci = {}
    for h in ta: loci.setdefault(h[0], []).append(h)
    out = []
    for sid, hs in sorted(loci.items(), key=lambda x: -sum(float(h[3]) for h in x[1]))[:2]:
        c0 = min(min(int(h[1]), int(h[2])) for h in hs); c1 = max(max(int(h[1]), int(h[2])) for h in hs); L = int(hs[0][4])
        a, b = max(1, c0 - 200000), min(L, c1 + 200000)
        nt = "".join(subprocess.run(["blastdbcmd", "-db", db, "-entry", sid, "-range", f"{a}-{b}"], capture_output=True,
                                    text=True).stdout.split("\n")[1:]).upper()
        hits = []
        for hebra, s in (("+", nt), ("-", str(Seq(nt).reverse_complement()))):
            for f in range(3):
                p = str(Seq(s[f:len(s) - ((len(s) - f) % 3)]).translate())
                hits += [(hebra, m.group(0)) for m in MOT.finditer(p)]
        bits = sum(float(h[3]) for h in hs)
        out.append(f"{sid}:{c0} (bit {bits:.0f}) motivos={[h[1] for h in hits] or 'ninguno'}")
    print(f"{sp:30s} " + "  ||  ".join(out), flush=True)
