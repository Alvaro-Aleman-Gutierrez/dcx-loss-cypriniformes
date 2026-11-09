#!/usr/bin/env python3
"""Busca en cada genoma el extremo 3' de la isoforma DCL de dclk1a (aa 280-346 de Ctenopharyngodon XP_051766222.1:
ultimo tramo compartido + cola propia EMDKAIIDSTA + stop). Si el HSP cruza la union 335/336 sin interrupcion,
la cola esta codificada inmediatamente despues del exon (retencion de intron), y se verifica el codon de parada."""
import subprocess, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
D = "${DATA}"; O = f"{D}/isoformas/cola"
s = next(str(r.seq) for r in SeqIO.parse(f"{D}/isoformas/dclk_isoformas.faa", "fasta") if "XP_051766222.1" in r.id)
q = s[279:]; open(f"{O}/q.faa", "w").write(f">dcl_tail\n{q}\n")   # 67 aa; cola = posiciones 57-67 del query
SP = sys.argv[1:]
for sp in SP:
    out = subprocess.run(["tblastn", "-query", f"{O}/q.faa", "-db", f"{D}/genomes/db_{sp}", "-evalue", "1e-3", "-seg", "no",
                          "-outfmt", "6 sseqid sstart send qstart qend pident evalue bitscore sseq", "-max_target_seqs", "20"],
                         capture_output=True, text=True).stdout.splitlines()
    hits = [l.split("\t") for l in out]
    hits.sort(key=lambda h: -float(h[7]))
    mejor = None
    for h in hits:
        if int(h[3]) <= 45 and int(h[4]) >= 64: mejor = h; break
    if not mejor:
        b = hits[0] if hits else None
        print(f"{sp:30s} SIN cola contigua; mejor HSP: " + (f"q{b[3]}-{b[4]} {b[5]}% {b[0]}:{b[1]}-{b[2]}" if b else "ninguno"))
        continue
    sid, a, b = mejor[0], int(mejor[1]), int(mejor[2])
    # codon siguiente al HSP (stop esperado)
    rec = subprocess.run(["blastdbcmd", "-db", f"{D}/genomes/db_{sp}", "-entry", sid, "-range",
                          f"{min(a,b)-3}-{max(a,b)+3}"], capture_output=True, text=True).stdout
    nt = "".join(rec.split("\n")[1:]).upper()
    if b < a: nt = str(Seq(nt).reverse_complement())
    sig = nt[-3:]
    print(f"{sp:30s} cola contigua q{mejor[3]}-{mejor[4]} id={mejor[5]}% {sid}:{a}-{b}  traduccion: ...{mejor[8][-14:]}  codon siguiente={sig} {'(STOP)' if sig in ('TAA','TAG','TGA') else ''}")
