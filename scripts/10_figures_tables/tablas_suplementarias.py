#!/usr/bin/env python3
"""Tablas suplementarias del manuscrito (Excel con una hoja por tabla + TSV por tabla).
  S1 Genome assemblies analysed (accession, level, size, DCX status, miniaturization)
  S2 Raw-read tests for the species of Fig. 2 (run, spots, read counts, observed/expected, 95% CI, P)
  S3 Read-based DCX census across Cypriniformes (DIAMOND; one row per run) -- se regenera al terminar el censo
  S4 Transcriptomic data (Iso-Seq / RNA-seq) used for the DCL isoform and DCX expression
  S5 Genes used for the species tree and the selection analyses (omega, purifying codons, RELAX K)
Uso: python3 tablas_suplementarias.py   (salida en manuscrito/suplemento/)
"""
import json, re, os, glob, math
import pandas as pd
from scipy.stats import chi2, poisson
P = "${PROJECT}"; D = "${DATA}"
O = f"{P}/manuscrito/suplemento"; os.makedirs(O, exist_ok=True)
S5T = pd.read_csv(f"{D}/sudasinghe2026/tablaS5_miniaturizacion.tsv", sep="\t")
S5T["sp"] = S5T.genero + "_" + S5T.especie.str.replace(r"\d+$", "", regex=True)
MINI = dict(zip(S5T.sp, S5T.miniatura)); MINI["Danio_choprae"] = MINI.get("Danio_choprai", MINI.get("Danio_choprae"))
def mini(sp): return MINI.get(sp, "")
def sp_txt(s): return s.replace("_sp_Pulau_Singkep", " sp. Pulau Singkep").replace("_", " ")

# ---------------------------------------------------------------- S1
acc = json.load(open(f"{D}/integracion/genomas_accesiones.json"))
WGS = {"Danio_albolineatus": ("GCA_903798035.1", "Scaffold", ""), "Danio_choprae": ("GCA_903798125.1", "Scaffold", ""),
       "Danio_jaintianensis": ("GCA_903798115.1", "Scaffold", ""), "Danio_tinwini": ("GCA_903798205.1", "Scaffold", ""),
       "Danio_nigrofasciatus": ("GCA_053539255.1", "Scaffold", "ASM5353925v1"),
       "Danionella_cerebrum": ("GCA_007224835.1", "Scaffold", "ASM722483v1"),
       "Danionella_cerebrum_v2": ("GCA_903798025.1", "Scaffold", "fDanTra1.1"),
       "Paedocypris_carbunculus": ("GCA_028128575.1", "Scaffold", "ASM2812857v1")}
cz = pd.read_csv(f"{D}/censo2/censo2_table.tsv", sep="\t").set_index("especie")
filas = []
for sp, v in acc.items():
    a = v["ensamblado"] if v["ensamblado"] else ("", "", "")
    a = WGS.get(sp, tuple(a)) if not a[0] else a
    if sp in cz.index:
        st = "present" if str(cz.loc[sp, "DCX"]).startswith(("presente", "parcial")) else "absent"
    else: st = ""
    nota = ""
    if sp == "Ictalurus_furcatus": st, nota = "present", "DCX missing from the assembly, present in raw reads (assembly gap)"
    if sp == "Danionella_cerebrum_v2": nota = "second assembly of D. cerebrum (replicate; excluded from counts)"
    grupo = "outgroup" if sp in ("Astyanax_mexicanus", "Ictalurus_furcatus") else "Cypriniformes"
    filas.append(dict(species=sp_txt(sp.replace("_v2", "")), assembly_accession=a[0], assembly_name=a[2], assembly_level=a[1],
                      size_Mb=v["tamano_Mb"], n_sequences=v["n_secuencias"], group=grupo, DCX=st,
                      miniaturization=mini(sp.replace("_v2", "")), note=nota))
S1 = pd.DataFrame(filas).sort_values(["group", "species"], ascending=[False, True])

# ---------------------------------------------------------------- S2 (counts from the raw-read test logs)
L = {"DCX": 357, "dclk1a": 770, "dclk2a": 826, "dclk3": 721}
LECT = [  # species, run, BioProject source, spots, DCX, dclk1a, dclk2a, dclk3, DCX present
    ("Danionella cerebrum", "ERR3284980", 15_000_000, 39, 76, 101, 104, True),
    ("Ictalurus furcatus", "SRR18966194", 40_000_000, 52, 84, 112, 174, True),
    ("Danionella dracula", "ERR2639754", 15_000_000, 1, 88, 75, 87, False),
    ("Danio tinwini", "ERR3284969", 40_000_000, 2, 102, 120, 150, False),
    ("Paedocypris micromegethes", "SRR11971489", 12_000_000, 2, 173, 113, 246, False),
    ("Paedocypris carbunculus", "SRR11971490", 40_000_000, 4, 698, 564, 867, False),
    ("Misgurnus anguillicaudatus", "SRR36982626", 35_000_000, 0, 186, 171, 208, False),
    ("Paramisgurnus dabryanus", "SRR23104325", 35_000_000, 1, 46, 52, 49, False)]
filas = []
for sp, run, spots, k, a1, a2, a3, pres in LECT:
    esp = (a1 + a2 + a3) / (L["dclk1a"] + L["dclk2a"] + L["dclk3"]) * L["DCX"]
    lo = chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0; hi = chi2.ppf(0.975, 2 * (k + 1)) / 2
    filas.append(dict(species=sp, run=run, spots_analysed=spots, reads_DCX=k, reads_dclk1a=a1, reads_dclk2a=a2, reads_dclk3=a3,
                      expected_DCX=round(esp, 1), obs_exp=round(k / esp, 3), CI95_low=round(lo / esp, 3), CI95_high=round(hi / esp, 3),
                      P_le_obs_if_present=f"{poisson.cdf(k, esp):.2g}", DCX_status="present" if pres else "absent"))
S2 = pd.DataFrame(filas)

# ---------------------------------------------------------------- S3 (census)
T = pd.read_csv(f"{D}/integracion/tabla_dcx_arbol.tsv", sep="\t")
T3 = T[T.fuente == "lecturas"].copy()
runinfo = {}
for f in glob.glob("/tmp/claude-1000/*/*/scratchpad/sra_cyp.tsv"):
    for l in open(f):
        c = l.rstrip("\n").split("\t")
        if len(c) >= 6: runinfo[c[2].split(",")[0]] = c[4]
T3["BioProject"] = T3.corrida.map(lambda r: runinfo.get(r, ""))
S3 = T3[["especie", "corrida", "BioProject", "voucher", "DCX_obs", "DCX_esp", "cociente", "estado", "punta", "miniatura"]].rename(
    columns=dict(especie="species", corrida="run", DCX_obs="reads_DCX", DCX_esp="expected_DCX", cociente="obs_exp",
                 estado="DCX_status", punta="tip_in_Sudasinghe_2026_tree", miniatura="miniaturization"))
spots = {}
for f in glob.glob(f"{D}/**/*.spots", recursive=True):
    m = re.search(r"Written (\d+) spots", open(f).read())
    if m: spots[os.path.basename(f)[:-6]] = max(int(m.group(1)), spots.get(os.path.basename(f)[:-6], 0))
S3.insert(4, "spots_analysed", S3.run.map(spots))
# runs between thresholds after 15 M spots were repeated once with 40 M; both results are reported
R15 = {}
for f in glob.glob(f"{D}/servidor_censo/resultados_15M_ambiguos/*.res"):
    t = open(f).read()
    R15[os.path.basename(f)[:-4]] = (int(re.search(r"DCX=(\d+)", t).group(1)), float(re.search(r"esperado=([\d.]+)", t).group(1)),
                                     float(re.search(r"obs/esp=([\d.]+)", t).group(1)))
S3["reads_DCX_15M"] = S3.run.map(lambda r: R15[r][0] if r in R15 else None)
S3["expected_DCX_15M"] = S3.run.map(lambda r: R15[r][1] if r in R15 else None)
S3["obs_exp_15M"] = S3.run.map(lambda r: R15[r][2] if r in R15 else None)
S3["note"] = S3.run.map(lambda r: "between thresholds with 15 M spots (columns *_15M); repeated with 40 M (main columns)" if r in R15 else "")
S3["species"] = S3.species.map(sp_txt); S3 = S3.sort_values("species")
S3["DCX_status"] = S3.DCX_status.map({"presente": "present", "ausente": "absent", "ambiguo": "ambiguous",
                                      "cobertura_baja": "low coverage"}).fillna(S3.DCX_status)

# ---------------------------------------------------------------- S4
S4 = pd.DataFrame([
    ("Paedocypris micromegethes", "SRR11971539", "PRJNA638472", "PacBio Iso-Seq (Sequel)", "whole body", "no", "DCL isoform structure (Fig. 4A)"),
    ("Danio rerio", "SRR37646263", "PRJNA1438185", "Illumina NovaSeq X, paired (read 1, 20 M)", "adult whole brain, wild type", "no", "dclk1a junction reads (Supplementary Fig. S4B)")]
    + [("Ctenopharyngodon idella", r, "PRJNA1306682", "Illumina NovaSeq X, paired", t, "yes", "DCX and dclk1a expression by tissue (Supplementary Fig. S4A, Table S6)")
       for r, t in [("SRR35031399", "brain"), ("SRR35031400", "brain"), ("SRR35029734", "gills"), ("SRR35029735", "gills"),
                    ("SRR35030603", "intestine"), ("SRR35030604", "intestine"), ("SRR35032171", "kidney"), ("SRR35032172", "kidney"),
                    ("SRR35030834", "liver"), ("SRR35030835", "liver"), ("SRR35031011", "spleen"), ("SRR35031012", "spleen")]],
    columns=["species", "run", "BioProject", "platform", "tissue", "DCX_in_genome", "use"])

# ---------------------------------------------------------------- S5
parts = re.findall(r"charset (\S+) = (\d+)-(\d+);", open(f"{D}/arbol_ext/partitions.nex").read())
arbol = {n: int(b) - int(a) + 1 for n, a, b in parts}
R = {r["gen"]: r for r in json.load(open(f"{D}/restriccion/restriccion_results.json"))}
genes = sorted(set(arbol) | {g for g in R if g.startswith("R")})
filas = []
for g in genes:
    r = R.get(g, {})
    filas.append(dict(gene=g, zebrafish_query=g.split("_", 1)[1] if "_" in g else "", in_species_tree="yes" if g in arbol else "no",
                      alignment_length_nt=arbol.get(g, ""), omega=round(r["omega"], 4) if "omega" in r else "",
                      purifying_codons_pct=round(100 * r["purif"], 1) if "purif" in r else "", RELAX_K=round(r["K"], 3) if "K" in r else "",
                      RELAX_p=f"{r['pK']:.3g}" if "pK" in r else "", category="random"))
for g in [g for g in R if not g.startswith("R")]:
    r = R[g]
    filas.append(dict(gene=g, zebrafish_query="", in_species_tree="no", alignment_length_nt="", omega=round(r.get("omega", float("nan")), 4),
                      purifying_codons_pct=round(100 * r.get("purif", float("nan")), 1), RELAX_K=round(r.get("K", float("nan")), 3),
                      RELAX_p=f"{r.get('pK', float('nan')):.3g}", category="DCX" if g == "DCX" else
                      ("synteny control" if g in ("hip1rb", "kntc1", "rwdd2b", "spata22", "zfpl1") else "DCX/DCLK family")))
S5 = pd.DataFrame(filas)

E6 = pd.read_csv(f"{D}/rnaseq/cid_expr/expresion_tejidos.tsv", sep="\t")
TEJ = {"cerebro": "brain", "intestino": "intestine", "branquias": "gills", "rinon": "kidney", "higado": "liver", "bazo": "spleen"}
E6["tejido"] = E6.tejido.map(TEJ)
S6 = E6.rename(columns={"tejido": "tissue", "corrida": "run", "lecturas_R1": "reads_R1_analysed", "DC_compartido": "reads_dclk1a_shared_DC",
    "kinasa": "reads_dclk1a_kinase", "DCL_especifico": "reads_dclk1a_DConly_exon", "union_larga": "junction_reads_to_kinase",
    "union_DCL": "junction_reads_to_DConly", "RPKM_kinasa": "RPKM_dclk1a_kinase", "RPKM_DCL": "RPKM_dclk1a_DConly_exon",
    "frac_DCL_uniones": "fraction_DConly_junctions", "DCX": "reads_DCX", "dclk2a": "reads_dclk2a", "elavl3": "reads_elavl3",
    "gfap": "reads_gfap", "actb1": "reads_actb1"})
S6.insert(0, "species", "Ctenopharyngodon idella")
S6 = pd.concat([S6, pd.DataFrame([dict(species="Danio rerio", tissue="brain", run="SRR37646263", reads_R1_analysed=20_000_000,
    junction_reads_to_kinase=67, junction_reads_to_DConly=0, fraction_DConly_junctions=0.0)])], ignore_index=True)
TITULOS = {"S1": "Supplementary Table S1. Genome assemblies analysed.",
           "S2": "Supplementary Table S2. Raw-read tests of DCX absence (Fig. 2).",
           "S3": "Supplementary Table S3. Read-based DCX census across Cypriniformes.",
           "S4": "Supplementary Table S4. Transcriptomic data.",
           "S5": "Supplementary Table S5. Genes used for the species tree and selection analyses.",
           "S6": "Supplementary Table S6. Expression of DCX and dclk1a isoforms in adult tissues."}
with pd.ExcelWriter(f"{O}/Supplementary_Tables_S1-S6.xlsx", engine="openpyxl") as xw:
    for k, df in [("S1", S1), ("S2", S2), ("S3", S3), ("S4", S4), ("S5", S5), ("S6", S6)]:
        df.to_excel(xw, sheet_name=k, index=False, startrow=1)
        xw.sheets[k].cell(row=1, column=1, value=TITULOS[k])
        df.to_csv(f"{O}/Table_{k}.tsv", sep="\t", index=False)
        print(f"{k}: {len(df)} filas")
print("escrito:", O)
