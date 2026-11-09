#!/usr/bin/env python3
"""Analisis 4: existe la isoforma corta de DCLK1/DCLK2 con solo dominios DC (tipo "DCL", sin cinasa)
en los linajes que perdieron DCX? En raton, DCL es funcionalmente redundante con DCX
(Dcx/Dclk1 doble KO; Deuel et al. 2006, Koizumi et al. 2006).

1. Para cada especie anotada (RefSeq), todos los genes DCLK1/DCLK2 (incluidos LOC "DCLK1-like") y TODAS
   sus proteinas RefSeq (isoformas).
2. Clasificacion por dominios con blastp contra dos segmentos de DCLK1 humano (NP_004725.1):
   DC = aa 1-330 (tandem DC + region SP inicial), KIN = aa 370-729 (dominio cinasa + cola).
   DC presente si el HSP cubre >= 150 aa de DC; cinasa si cubre >= 150 aa de KIN.
   Clases: largo (DC+KIN), DC-solo (tipo DCL), KIN-solo (tipo CPG16/corto).
"""
import urllib.request, urllib.parse, json, time, re, os, subprocess, sys, collections
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
O = "${DATA}/isoformas"; os.makedirs(O, exist_ok=True)
def get(u, **p):
    for i in range(6):
        try: return urllib.request.urlopen(E + u, data=urllib.parse.urlencode(p).encode(), timeout=90).read().decode()
        except Exception: time.sleep(3 + 3 * i)
    raise RuntimeError(u)
ESP = {  # especie: tiene DCX?
    "Homo sapiens": "ref", "Mus musculus": "ref",
    "Astyanax mexicanus": True, "Ictalurus punctatus": True,
    "Gobio gobio": True, "Pimephales promelas": True, "Ctenopharyngodon idella": True, "Labeo rohita": True,
    "Carassius auratus": True, "Cyprinus carpio": True, "Triplophysa rosa": True,
    "Danio rerio": False, "Misgurnus anguillicaudatus": False, "Paramisgurnus dabryanus": False,
}
fa = f"{O}/dclk_isoformas.faa"; meta = f"{O}/dclk_isoformas.tsv"
if not os.path.exists(meta):
    fo, mo = open(fa, "w"), open(meta, "w")
    mo.write("especie\tgen_id\tsimbolo\tdescripcion\tproteina\tlongitud\n")
    for sp in ESP:
        term = (f'(dclk1*[sym] OR dclk2*[sym] OR "doublecortin like kinase"[title] OR "DCLK1"[title] OR "DCLK2"[title]) '
                f'AND "{sp}"[orgn] AND alive[prop] AND "genetype protein coding"[prop]')
        ids = re.findall(r"<Id>(\d+)</Id>", get("esearch.fcgi", db="gene", term=term, retmax=60))
        if not ids: print("sin genes", sp, flush=True); continue
        sm = json.loads(get("esummary.fcgi", db="gene", id=",".join(ids), retmode="json"))["result"]
        for g in ids:
            r = sm[g]; desc = r["description"]; sym = r["name"]
            if not re.search(r"(?i)dclk[12]|doublecortin.like kinase [12]|DCLK[12]", sym + " " + desc): continue
            pids = re.findall(r"<Id>(\d+)</Id>", get("elink.fcgi", dbfrom="gene", db="protein",
                                                       linkname="gene_protein_refseq", id=g))
            if not pids: continue
            txt = get("efetch.fcgi", db="protein", id=",".join(pids), rettype="fasta", retmode="text")
            for bloque in txt.split(">")[1:]:
                cab, *sq = bloque.split("\n"); acc = cab.split()[0]; s = "".join(sq)
                tag = f"{sp.replace(' ', '_')}|{sym}|{acc}"
                fo.write(f">{tag}\n{s}\n"); mo.write(f"{sp}\t{g}\t{sym}\t{desc}\t{acc}\t{len(s)}\n")
            print(f"  {sp:28s} {sym:14s} {len(pids):3d} isoformas  {desc}", flush=True)
            time.sleep(0.4)
    fo.close(); mo.close()
# referencias: DCLK1 humano
ref = f"{O}/ref.faa"
if not os.path.exists(ref):
    s = "".join(get("efetch.fcgi", db="protein", id="NP_004725.1", rettype="fasta", retmode="text").split("\n")[1:])
    open(ref, "w").write(f">DC\n{s[:330]}\n>KIN\n{s[369:]}\n")
subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "prot", "-out", f"{O}/refdb"], check=True, capture_output=True)
out = subprocess.run(["blastp", "-query", fa, "-db", f"{O}/refdb", "-evalue", "1e-5", "-outfmt",
                      "6 qseqid sseqid sstart send", "-max_hsps", "5"], capture_output=True, text=True, check=True).stdout
cob = collections.defaultdict(lambda: {"DC": set(), "KIN": set()})
for l in out.splitlines():
    q, s, a, b = l.split("\t"); cob[q][s].update(range(int(a), int(b) + 1))
filas = []
for l in open(meta).read().splitlines()[1:]:
    sp, g, sym, desc, acc, L = l.split("\t"); tag = f"{sp.replace(' ', '_')}|{sym}|{acc}"
    dc, kin = len(cob[tag]["DC"]) >= 150, len(cob[tag]["KIN"]) >= 150
    clase = "largo" if dc and kin else ("DC-solo" if dc else ("KIN-solo" if kin else "otro"))
    filas.append((sp, sym, desc, acc, int(L), clase))
with open(f"{O}/isoformas_clasificadas.tsv", "w") as fo:
    fo.write("especie\tsimbolo\tdescripcion\tproteina\tlongitud\tclase\n")
    for f in filas: fo.write("\t".join(map(str, f)) + "\n")
# resumen por especie x gen
res = collections.defaultdict(collections.Counter)
for sp, sym, desc, acc, L, clase in filas: res[(sp, sym)][clase] += 1
print("\nespecie\tDCX\tgen\tlargo\tDC-solo\tKIN-solo\totro")
for (sp, sym), c in sorted(res.items(), key=lambda x: (list(ESP).index(x[0][0]), x[0][1])):
    print(f"{sp}\t{ESP[sp]}\t{sym}\t{c['largo']}\t{c['DC-solo']}\t{c['KIN-solo']}\t{c['otro']}")
print("FIN isoformas")
