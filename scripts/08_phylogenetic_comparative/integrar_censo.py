#!/usr/bin/env python3
"""Integra el estado de DCX (genomas + lecturas SRA) con el arbol calibrado de Sudasinghe et al. 2026 (PNAS;
Zenodo 10.5281/zenodo.23008328, CC-BY; 335 taxones, BEAST) y su Dataset S1 (miniaturizacion, SL maxima, genoma, ploidia).

Emparejamiento con las puntas del arbol (etiquetas Suborden_Familia_Genero_especie_VOUCHER|ACCESION):
  1. corridas SRA -> Voucher_reads del Dataset S1 (columna Accession) -> punta que termina en ese voucher;
  2. si no, por Genero_especie (sinonimias ortograficas en SINON); varias puntas de la misma especie: se usa la primera
     y se anota.
Regla de estado a partir del cociente obs/esp de DCX (lecturas; parse_dmnd.py): <= 0,10 ausente; >= 0,30 presente;
entre ambos 'ambiguo'; esperado < 10 lecturas -> 'cobertura_baja'. Genomas: estado verificado del censo por RBH
(censo2), con la correccion de I. furcatus (presente en lecturas).
Salida: ${DATA}/integracion/tabla_dcx_arbol.tsv y especies_sin_punta.tsv
"""
import re, glob, os, collections
import pandas as pd
D = "${DATA}"; S = f"{D}/sudasinghe2026"; O = f"{D}/integracion"; os.makedirs(O, exist_ok=True)
SINON = {"Danio_choprae": "Danio_choprai"}
tre = open(f"{S}/Sudasinghe-etal.-BEAST-time.calibrated.phylogeny-335taxa.tre").read()
tips = re.search(r"Taxlabels(.*?);", tre, re.S).group(1).split()
def gsp(tip):
    p = tip.split("_"); return f"{p[2]}_{p[3]}"
por_sp = collections.defaultdict(list)
for t in tips: por_sp[gsp(t)].append(t)
s1 = pd.read_excel(f"{S}/Sudasinghe-etal.-Dataset_S1.xlsx")
acc2vou = dict(zip(s1["Accession"].astype(str), s1["Voucher_reads"].astype(str)))
rasgo_vou = {str(r.Voucher_reads): r for r in s1.itertuples()}
rasgo_sp = {f"{r.Genus}_{r.Species}": r for r in s1.itertuples()}

def estado_lecturas(ratio, esp):
    if esp < 10: return "cobertura_baja"
    return "ausente" if ratio <= 0.10 else ("presente" if ratio >= 0.30 else "ambiguo")

filas = []
# --- lecturas SRA (cola prioritaria + censo completo + validaciones con corridas propias)
res = []
for f in [f"{D}/sra_ext/sra_ext_results.tsv"] + sorted(glob.glob(f"{D}/servidor_censo/resultados/*.res")):
    for l in open(f):
        if "\t" not in l: continue
        c = l.rstrip("\n").split("\t"); kv = dict(x.split("=", 1) for x in c[1:] if "=" in x)
        res.append((c[0], kv))
corrida_val = {"VAL1furcatus": "SRR18966194", "VAL2paramisgurnus": "SRR23104325"}
for etq, kv in res:
    nombre, run = etq.split("|"); run = corrida_val.get(run, run)
    sp = re.sub(r"^CONTROL[+-]_", "", nombre); sp = "_".join(sp.split("_")[:2]).replace("sp.", "sp")
    ratio, esp = float(kv["obs/esp"]) if kv["obs/esp"] != "nan" else float("nan"), float(kv["esperado"])
    vou = acc2vou.get(run)
    tip = next((t for t in tips if vou and t.endswith("_" + vou)), None)
    if tip is None and por_sp.get(SINON.get(sp, sp)): tip = por_sp[SINON.get(sp, sp)][0]
    filas.append(dict(especie=sp, fuente="lecturas", corrida=run, voucher=vou, punta=tip,
                      DCX_obs=int(kv["DCX"]), DCX_esp=esp, cociente=ratio, estado=estado_lecturas(ratio, esp)))
# --- genomas (censo2, ortologia por RBH)
cz = pd.read_csv(f"{D}/censo2/censo2_table.tsv", sep="\t")
for r in cz.itertuples():
    sp = r.especie
    if sp == "Danionella_cerebrum_v2": continue
    est = "presente" if str(r.DCX).startswith(("presente", "parcial")) else "ausente"
    if sp == "Ictalurus_furcatus": est = "presente"     # hueco de ensamblado; DCX presente en lecturas (Fig. 2)
    tips_sp = por_sp.get(SINON.get(sp, sp), [])
    tip = next((t for t in tips_sp if "_GC" in t), tips_sp[0] if tips_sp else None)
    filas.append(dict(especie=sp, fuente="genoma", corrida="", voucher="", punta=tip, DCX_obs=None, DCX_esp=None,
                      cociente=None, estado=est))
T = pd.DataFrame(filas)
s1i = s1.reset_index(drop=True)
idx_vou = {str(v): i for i, v in enumerate(s1i["Voucher_reads"])}
idx_sp = {}
for i, r in s1i.iterrows(): idx_sp.setdefault(f"{r.Genus}_{r.Species}", i)
def fila_s1(row):
    i = idx_vou.get(str(row["voucher"]))
    if i is None: i = idx_sp.get(SINON.get(row["especie"], row["especie"]))
    return s1i.loc[i] if i is not None else None
CAMPOS = {"Family": "familia", "Miniature": "miniatura", "SL.max.cm": "SL_max_cm",
          "mean_k31_genomesize_MB": "genoma_Mb", "Ploidy": "ploidia"}
for k, v in CAMPOS.items():
    T[v] = [(lambda f: f[k] if f is not None else None)(fila_s1(r)) for _, r in T.iterrows()]
# Codificacion de miniaturizacion de los 309 taxones (Tabla S5 del suplemento; extraida a tablaS5_miniaturizacion.tsv).
# Tiene prioridad sobre el Dataset S1 porque cubre tambien los genomas de NCBI. Busqueda por voucher/accesion y luego especie.
s5 = pd.read_csv(f"{S}/tablaS5_miniaturizacion.tsv", sep="\t")
s5["sp"] = s5.genero + "_" + s5.especie.str.replace(r"\d+$", "", regex=True)
s5_vou = {str(r.voucher): r for r in s5.itertuples()}; s5_sp = {}
for r in s5.itertuples(): s5_sp.setdefault(r.sp, r)
def fila_s5(row):
    tip = row["punta"] if isinstance(row["punta"], str) else ""
    for v in [str(row["voucher"])] + [p for p in tip.split("_") if p.startswith(("LR", "GC")) or p[:3].isdigit()]:
        if v in s5_vou: return s5_vou[v]
    gc = re.search(r"(GC[AF]_\d+\.\d+)", tip)
    if gc and gc.group(1) in s5_vou: return s5_vou[gc.group(1)]
    return s5_sp.get(SINON.get(row["especie"], row["especie"]))
f5 = [fila_s5(r) for _, r in T.iterrows()]
T["miniatura"] = [f.miniatura if f is not None else m for f, m in zip(f5, T["miniatura"])]
T["SL_max_cm"] = [f.SL_max_cm if f is not None else m for f, m in zip(f5, T["SL_max_cm"])]
T["ploidia"] = [f.ploidia if f is not None else m for f, m in zip(f5, T["ploidia"])]
# conflictos dentro de una misma punta (p. ej. genoma y lecturas)
conf = T.dropna(subset=["punta"]).groupby("punta")["estado"].agg(lambda x: "|".join(sorted(set(x))))
T["estado_punta"] = T["punta"].map(conf)
T.to_csv(f"{O}/tabla_dcx_arbol.tsv", sep="\t", index=False)
T[T.punta.isna()].to_csv(f"{O}/especies_sin_punta.tsv", sep="\t", index=False)
print(f"registros: {len(T)} | con punta en el arbol: {T.punta.notna().sum()} | puntas distintas: {T.punta.nunique()}")
print(T.groupby(["fuente", "estado"]).size().to_string())
print("\nsin punta:", ", ".join(T[T.punta.isna()].especie))
c = T.dropna(subset=["punta"]).drop_duplicates("punta")
print("\nconflictos genoma/lecturas en la misma punta:", (c.estado_punta.str.contains(r"\|")).sum())
print("\nausentes con punta:"); print(c[c.estado_punta.str.contains("ausente")][["especie", "fuente", "estado_punta", "miniatura"]].to_string(index=False))
