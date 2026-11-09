#!/usr/bin/env python3
"""Punto 1: presencia de DCX en todos los vertebrados con genoma anotado (RefSeq).
Para cada especie anotada busca en NCBI Gene registros 'doublecortin' (DCX/LOC) y registra
si hay un gen DCX (descripción 'doublecortin' o 'neuronal migration protein doublecortin', sin 'like kinase'/'domain').
Salida: p1/vertebrate_dcx_presence.tsv
"""
import urllib.request, urllib.parse, json, time, re, os, sys
OUT = "${DATA}/p1"; os.makedirs(OUT, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(url, tries=8):
    for i in range(tries):
        try: return urllib.request.urlopen(url, timeout=90).read().decode()
        except Exception as e: time.sleep(20 * (i + 1))
    raise RuntimeError(url)
# 1) especies vertebradas con anotación RefSeq
spp = {}; tok = ""
while True:
    u = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/7742/dataset_report?page_size=1000"
         "&filters.has_annotation=true&filters.reference_only=true" + (f"&page_token={tok}" if tok else ""))
    d = json.loads(get(u))
    for r in d.get("reports", []):
        o = r["organism"]; spp[o["tax_id"]] = (o["organism_name"], r["accession"])
    tok = d.get("next_page_token")
    if not tok: break
print("especies anotadas:", len(spp), flush=True)
done = set()
fn = f"{OUT}/vertebrate_dcx_presence.tsv"
if os.path.exists(fn): done = {l.split("\t")[0] for l in open(fn)}
fo = open(fn, "a")
if not done: fo.write("taxid\tespecie\taccesion\tn_DCX\tgenes_DCX\tn_doublecortin_total\n")
for tx, (name, acc) in spp.items():
    if str(tx) in done: continue
    term = f'doublecortin[title] AND txid{tx}[orgn:noexp] AND alive[prop]'
    ids = re.findall(r"<Id>(\d+)</Id>", get(E + "esearch.fcgi?" + urllib.parse.urlencode(dict(db="gene", term=term, retmax=100))))
    dcx = []
    if ids:
        sm = json.loads(get(E + "esummary.fcgi?" + urllib.parse.urlencode(dict(db="gene", id=",".join(ids), retmode="json"))))["result"]
        for u in sm["uids"]:
            desc = sm[u]["description"].lower()
            if ("like kinase" in desc or "domain" in desc or "kinase" in desc): continue
            if desc in ("doublecortin", "neuronal migration protein doublecortin") or desc.startswith("neuronal migration protein doublecortin"):
                dcx.append(f'{sm[u]["name"]}:{u}')
    fo.write(f"{tx}\t{name}\t{acc}\t{len(dcx)}\t{','.join(dcx)}\t{len(ids)}\n"); fo.flush()
    time.sleep(0.4)
print("FIN punto 1", flush=True)
