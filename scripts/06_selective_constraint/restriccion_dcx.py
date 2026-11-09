#!/usr/bin/env python3
"""Analisis 5: restriccion selectiva de DCX en los linajes que lo conservan.
Todos los genes (DCX, 15 miembros de la familia, 5 controles de sintenia y 49 genes al azar) se podan a las MISMAS
9 especies con DCX, con el arbol de especies final (mf10) podado, para que los omegas sean comparables.
 (a) FEL -> omega global (MG94xREV) y fraccion de sitios con seleccion purificadora (p<0,1).
 (b) RELAX: Test = Danionella cerebrum + Rasbora (danionidos que conservan DCX; hermanos de linajes que lo perdieron)
     vs Reference = resto. Pregunta: se relajo DCX antes de perderse? Se calibra con los genes al azar.
"""
import glob, os, subprocess, json, sys, io
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO, Phylo
R1 = "${PROJECT}/relax"; R2 = "${DATA}/p5/relax"
O = "${DATA}/restriccion"; os.makedirs(O, exist_ok=True)
ARBOL = "${DATA}/arbol_ext/mf10.treefile"
TAX = ["Danionella_cerebrum", "Rasbora_lateristriata", "Pimephales_promelas", "Ctenopharyngodon_idella", "Labeo_rohita",
       "Gobio_gobio", "Barbatula_barbatula", "Triplophysa_rosa", "Astyanax_mexicanus"]
TEST = {"Danionella_cerebrum", "Rasbora_lateristriata"}
HY = "${HOME}/anaconda3/envs/hyphy/bin/hyphy"; P = int(os.environ.get("PAR", "4"))
genes = {os.path.basename(f)[:-10]: f for f in glob.glob(f"{R1}/*.codon.fas")}
genes.update({os.path.basename(f)[:-10]: f for f in glob.glob(f"{R2}/*.codon.fas")})
STOP = {"TAA", "TAG", "TGA"}
def preparar(g):
    d = {r.id: str(r.seq).upper() for r in SeqIO.parse(genes[g], "fasta") if r.id in TAX}
    if len(d) < 7 or "Astyanax_mexicanus" not in d: return None
    L = len(next(iter(d.values())))
    keep = []
    for i in range(0, L - 2, 3):
        cods = [s[i:i + 3] for s in d.values()]
        if all(c == "---" for c in cods): continue
        if any(c in STOP for c in cods): continue
        keep.append(i)
    fas = f"{O}/{g}.fas"
    with open(fas, "w") as fo:
        for k, s in d.items():
            fo.write(f">{k}\n" + "".join(s[i:i + 3].replace("N", "-") if "-" in s[i:i + 3] else s[i:i + 3] for i in keep) + "\n")
    t = Phylo.read(ARBOL, "newick")
    for x in list(t.get_terminals()):
        if x.name not in d: t.prune(x)
    for c in t.find_clades(): c.confidence = None; c.comment = None
    h = io.StringIO(); Phylo.write(t, h, "newick", plain=False); nw = h.getvalue().strip()
    open(f"{O}/{g}.nwk", "w").write(nw + "\n")
    for x in TEST & set(d): nw = nw.replace(x + ":", x + "{Test}:")
    open(f"{O}/{g}.relax.nwk", "w").write(nw + "\n")
    return fas
def correr(g):
    fas = preparar(g)
    if not fas: return g, None
    for an, extra, out in [("fel", ["--tree", f"{O}/{g}.nwk"], f"{O}/{g}.FEL.json"),
                           ("relax", ["--tree", f"{O}/{g}.relax.nwk", "--test", "Test"],
                            f"{O}/{g}.RELAX.json")]:
        if os.path.exists(out) and os.path.getsize(out) > 0: continue
        with open(out.replace(".json", ".log"), "w") as lg:
            subprocess.run([HY, "CPU=2", an, "--alignment", fas, *extra, "--output", out], stdout=lg, stderr=lg)
    return g, fas
def leer(g):
    r = {"gen": g}
    try:
        F = json.load(open(f"{O}/{g}.FEL.json"))
        rd = F["fits"]["Global MG94xREV"]["Rate Distributions"]
        r["omega"] = next(v[0][0] for k, v in rd.items() if "non-synonymous/synonymous" in k)
        filas = F["MLE"]["content"]["0"]; n = len(filas)
        r["sitios"] = n; r["purif"] = sum(1 for x in filas if x[1] < x[0] and x[4] < 0.1) / n
    except Exception as e: r["err_fel"] = str(e)[:60]
    try:
        X = json.load(open(f"{O}/{g}.RELAX.json"))["test results"]; r["K"] = X["relaxation or intensification parameter"]; r["pK"] = X["p-value"]
    except Exception as e: r["err_relax"] = str(e)[:60]
    return r
if __name__ == "__main__":
    lista = sorted(genes, key=lambda g: (g != "DCX", g))
    with ThreadPoolExecutor(P) as ex:
        for g, ok in ex.map(correr, lista): print(f"  {g}: {'ok' if ok else 'omitido'}", flush=True)
    res = [leer(g) for g in lista if os.path.exists(f"{O}/{g}.fas")]
    json.dump(res, open(f"{O}/restriccion_results.json", "w"), indent=1)
    azar = [r for r in res if r["gen"].startswith("R") and "omega" in r]
    for r in res:
        if r["gen"].startswith("R") or "omega" not in r: continue
        pct = 100 * sum(a["omega"] <= r["omega"] for a in azar) / len(azar)
        pk = ""
        if "K" in r:
            ka = [a["K"] for a in azar if "K" in a]; pk = f"K={r['K']:.2f} (p={r['pK']:.3g}; {100*sum(k<=r['K'] for k in ka)/len(ka):.0f}% de azar con K<=)"
        print(f"{r['gen']:14s} omega={r['omega']:.4f} (percentil {pct:.0f} entre {len(azar)} al azar) purif={r.get('purif',0):.2f} {pk}")
    om = sorted(a["omega"] for a in azar); print(f"azar: omega mediana={om[len(om)//2]:.4f} min={om[0]:.4f} max={om[-1]:.4f}")
    print("FIN restriccion")
