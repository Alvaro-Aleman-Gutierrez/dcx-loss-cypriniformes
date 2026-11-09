#!/usr/bin/env python3
"""Manuscript figures, English version, renumbered 2026-11-06.
  Fig. 1  phylogeny + DCX/DCLK family census (former Figs. 1 and 2, unified)
  Fig. 2  raw-read support: observed/expected DCX reads
  Fig. 3  microsynteny of the DCX-pak3loc locus
  Fig. S1 DCX annotation across RefSeq vertebrates | Fig. S2 topology tests | Fig. S3 RELAX null distribution
No titles or captions inside the images: legends are written to figuras/figure_legends.md.
Palette derived from Leon et al. 2024 (Mol Biol Evol 41:msae166, CC BY 4.0); HEX values measured on its figures
(scripts/paleta_leon.py). Fixed meaning in every figure:
  #664BBC DCX present (filled diamond) | #FA41CA DCX absent / loss (open diamond, x) | #686868 node support (circle) | #66FEFF syntenic neighbour
  #71D226 outgroup | #E5C0DA lost region / lineage without DCX | #CAD9F6 gene present (census) | #FEDC00 assembly gap
Never colour alone. Accessibility check: scripts/verificar_figuras.py ; palette check: scripts/validar_paleta.py
Spanish version of this script (former numbering): scripts/figuras_es.py.bak
Usage: python3 -I figuras.py [fig1 fig2 fig3 figS1 figS2 figS3]   (no arguments: all)
"""
import os, sys, csv, json, glob, re, math, collections, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.text
import matplotlib.patheffects as pe
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon
from matplotlib.lines import Line2D
from Bio import Phylo
P = "${PROJECT}"; D = "${DATA}"
F = f"{P}/figuras"; os.makedirs(F, exist_ok=True)
# ---------------------------------------------------------------- palette (Leon et al. 2024, measured)
PRES, AUS, PARC, EXT, SINT = "#664BBC", "#FA41CA", "#FEDC00", "#71D226", "#66FEFF"
SINORT, INK, INK2, GRID, CROM, HUECO, BG = "#C4BAB9", "#000000", "#525252", "#EBEBEB", "#EBEBEB", "#E5C0DA", "#FFFFFF"
PRES_T, PARC_E, BANDA_AUS, BANDA_EXT, CINTA = "#CAD9F6", "#B4A7D5", "#E5C0DA", "#D0E9B5", "#C4BAB9"
def _Lstar(h):
    c = [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    c = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in c]
    return 116 * (0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]) ** (1 / 3) - 16
def borde(col):
    """Very light fill (L* > 85) gets a dark outline so it reads on white and in greyscale."""
    return dict(edgecolor=INK2, linewidth=0.5) if _Lstar(col) > 85 else dict(linewidth=0)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5, "axes.labelsize": 7.5, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "axes.facecolor": BG, "figure.facecolor": BG, "savefig.facecolor": BG, "text.color": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "legend.frameon": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
    "mathtext.fontset": "custom", "mathtext.rm": "DejaVu Sans", "mathtext.it": "DejaVu Sans:italic", "mathtext.bf": "DejaVu Sans:bold",
    "mathtext.sf": "DejaVu Sans:italic"})   # it() uses \mathsf so digits are italic too
NO_CAP = re.compile(r"^(n = |p = |p=|dN/dS|kb\b|Mb\b|Ma\b|K\b|ω)")
def _cap(t):
    """Sentence case for figure text: first letter upper-case, unless the text is a gene/species name (italic or
    mathtext) or starts with a symbol/unit (n, p, dN/dS, kb, Mb, Ma)."""
    s_ = t.get_text()
    if not s_ or s_.startswith("$") or t.get_fontstyle() == "italic" or NO_CAP.match(s_): return
    i = next((k for k, ch in enumerate(s_) if ch.isalnum()), None)
    if i is None or not s_[i].isalpha() or s_[i].isupper(): return
    t.set_text(s_[:i] + s_[i].upper() + s_[i + 1:])
def mayusculas(fig):
    fig.canvas.draw()
    for ax in fig.axes:   # freeze tick labels so the change survives the next draw
        for eje, get, set_ in (("x", ax.get_xticklabels, ax.set_xticklabels), ("y", ax.get_yticklabels, ax.set_yticklabels)):
            tl = get()
            if any(t.get_text() and t.get_text()[0].isalpha() and t.get_fontstyle() != "italic" and t.get_text()[0].islower() for t in tl):
                loc = ax.get_xticks() if eje == "x" else ax.get_yticks()
                for t in tl: _cap(t)
                props = dict(rotation=tl[0].get_rotation(), fontsize=tl[0].get_fontsize(), fontstyle=tl[0].get_fontstyle())
                (ax.set_xticks if eje == "x" else ax.set_yticks)(loc); set_([t.get_text() for t in tl], **props)
    for t in fig.findobj(matplotlib.text.Text): _cap(t)
def guardar(fig, nombre):
    mayusculas(fig)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{F}/{nombre}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig); print("written:", nombre)
AUSENTE = {"Danio_rerio", "Danio_aesculapii", "Danio_kyathit", "Danio_albolineatus", "Danio_tinwini", "Danio_choprae",
           "Danio_jaintianensis", "Danio_nigrofasciatus", "Danionella_dracula", "Paedocypris_micromegethes",
           "Paedocypris_carbunculus", "Paedocypris_sp_Pulau_Singkep", "Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus"}
EXTERNO = {"Astyanax_mexicanus", "Ictalurus_furcatus"}
LINAJES = [("Danio", lambda s: s.startswith("Danio_")), ("D. dracula", lambda s: s == "Danionella_dracula"),
           ("Cobitidae", lambda s: s in ("Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus")),
           ("Paedocypris", lambda s: s.startswith("Paedocypris"))]
# final species tree: ModelFinder+MERGE (147 -> 32 partitions), 10 independent ML searches, RF = 0
ARBOL = os.environ.get("ARBOL_FIG1", f"{D}/arbol_ext/mf10.treefile")
def cient(p):
    m, e = f"{p:.1e}".split("e"); return rf"${m}\times10^{{{int(e)}}}$"
def bonito(s):
    s = s.replace("_sp_Pulau_Singkep", " sp. Singkep").replace("_", " ")
    p = s.split(); return (p[0][0] + ". " + " ".join(p[1:])) if len(p) > 1 else s
def marca_estado(ax, x, y, presente, ms=4.6, **kw):
    """DCX status: filled diamond = present, open diamond = absent (circles are reserved for node support)."""
    if presente: ax.plot([x], [y], "D", ms=ms, color=PRES, mec=PRES, mew=0.8, **kw)
    else:        ax.plot([x], [y], "D", ms=ms, mfc=BG, mec=AUS, mew=1.4, **kw)
def it(t):
    """Italic run for mixed labels (mathtext, same sans font)."""
    rep = {". ": r"{.}\,", ".": "{.}", ", ": r"{,}\ ", ",": "{,}", " ": r"\ "}
    return r"$\mathsf{" + re.sub(r"\. |, |[., ]", lambda m: rep[m.group(0)], t) + "}$"
def limpiar(ax, left=True):
    for s_ in ("top", "right") + (("left",) if left else ()): ax.spines[s_].set_visible(False)

# =============================================================== FIG 1: phylogeny (A) + family census (B)
GRUPOS_FAM = [("DCX", ["DCX"]), ("DCLK", ["dclk1a", "dclk1b", "dclk2a", "dclk2b", "dclk3"]),
              ("DCDC2", ["dcdc2", "dcdc2b", "dcdc2c"]), ("RP1L1", ["rp1l1a", "rp1l1b"]),
              ("synteny controls", ["kntc1", "hip1rb", "zfpl1", "rwdd2b", "spata22"])]
# node support: one dot per node, size AND opacity scale with UFBoot (70 -> 100); dark red #9A0000 measured in
# Leon et al. 2024 (GO bars, Fig. 3a); validated against every other colour of the figure, incl. CVD simulations
SOP_COL = os.environ.get("SOP_COL", "#686868"); SOP_MIN, SOP_MAX = 70.0, 100.0   # medium grey measured in Leon et al. 2024 (map land, Fig. 6)
def sop_estilo(uf):
    f = min(max((uf - SOP_MIN) / (SOP_MAX - SOP_MIN), 0.0), 1.0)
    return dict(ms=2.2 + 2.6 * f, alpha=0.35 + 0.65 * f)
def fig1():
    rows = list(csv.reader(open(f"{D}/censo2/censo2_table.tsv"), delimiter="\t")); hdr = rows[0][1:]
    dat = {r[0]: {g: v.split("(")[0] for g, v in zip(hdr, r[1:])} for r in rows[1:] if r[0] != "Danionella_cerebrum_v2"}
    t = Phylo.read(ARBOL, "newick"); t.root_at_midpoint(); t.ladderize(reverse=True)
    lados = [set(x.name for x in c.get_terminals()) for c in t.root.clades]
    assert any(s == EXTERNO for s in lados), "midpoint root does not separate the outgroup"
    tips = t.get_terminals(); orden = [x.name for x in tips]; Y = {n: i for i, n in enumerate(orden)}
    assert set(orden) == set(dat)
    prof = t.depths(); X = max(prof.values()); XEND = X * 1.04
    def yp(cl): return Y[cl.name] if cl.is_terminal() else sum(yp(h) for h in cl.clades) / len(cl.clades)
    def soporte(cl):
        """(SH-aLRT, UFBoot) from IQ-TREE node labels, or None."""
        if cl.is_terminal() or not cl.name or "/" not in str(cl.name): return None
        a, b = str(cl.name).split("/")[:2]; return float(a), float(b)
    es_ext = lambda cl: all(z.name in EXTERNO for z in cl.get_terminals())
    cols, xcol, xg, x = [], [], [], 0
    for nom, gs_ in GRUPOS_FAM:
        xs = []
        for g in gs_: cols.append(g); xcol.append(x); xs.append(x); x += 1
        xg.append((nom, min(xs), max(xs))); x += 0.6
    XM = x - 0.6
    fig = plt.figure(figsize=(7.2, 5.0))
    gsp = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.46, 1.38], wspace=0.0, top=0.86, bottom=0.13, left=0.01, right=0.88)
    at = fig.add_subplot(gsp[0]); al = fig.add_subplot(gsp[1], sharey=at); am = fig.add_subplot(gsp[2], sharey=at)
    bandas = [(nom, f, BANDA_AUS) for nom, f in LINAJES] + [("outgroup", lambda s: s in EXTERNO, BANDA_EXT)]
    for nom, f, col in bandas:          # A: lineage bands from each lineage stem
        sel = [c for c in tips if f(c.name)]
        mrca = sel[0] if len(sel) == 1 else t.common_ancestor(sel)
        cam = t.get_path(mrca); padre = cam[-2] if len(cam) > 1 else t.root
        x0 = prof[padre] + 0.004; ys = [Y[c.name] for c in sel]
        at.add_patch(FancyBboxPatch((x0, min(ys) - 0.44), XEND - x0, max(ys) - min(ys) + 0.88,
                     boxstyle="round,pad=0,rounding_size=0.006", mutation_aspect=60, facecolor=col, edgecolor="none",
                     alpha=0.55, zorder=0))
    n_perd = 0
    def rec(cl):
        nonlocal n_perd
        if cl.is_terminal(): return
        x0 = prof[cl]; ys = [yp(h) for h in cl.clades]
        at.plot([x0, x0], [min(ys), max(ys)], color=EXT if es_ext(cl) else INK, lw=1.8 if es_ext(cl) else 0.9,
                solid_capstyle="round", zorder=2)
        for h, yh in zip(cl.clades, ys):
            perd = all(z.name in AUSENTE for z in h.get_terminals()) and not all(z.name in AUSENTE for z in cl.get_terminals())
            at.plot([x0, prof[h]], [yh, yh], color=AUS if perd else (EXT if es_ext(h) else INK),
                    lw=2.2 if perd else (1.8 if es_ext(h) else 0.9), solid_capstyle="round", zorder=3)
            if perd:
                n_perd += 1; at.plot([(x0 + prof[h]) / 2], [yh], marker="X", ms=7, color=AUS, mec=BG, mew=1.0, zorder=5)
            s = soporte(h)
            if s is not None:
                at.plot(prof[h], yh, "o", color=SOP_COL, mec="none", zorder=6, **sop_estilo(s[1]))
            rec(h)
    rec(t.root)
    PUNT = dict(color="#BBBBBB", lw=0.5, ls=(0, (1, 2)), zorder=1)
    for c in tips:
        yy = Y[c.name]; xx = prof[c]
        marca_estado(at, xx + 0.004, yy, c.name not in AUSENTE, ms=4.0, zorder=6)
        at.plot([xx + 0.009, XEND], [yy, yy], **PUNT); al.plot([0, 1], [yy, yy], **PUNT); am.plot([-0.7, -0.47], [yy, yy], **PUNT)
        al.text(0.5, yy, bonito(c.name), ha="center", va="center", fontsize=6.4, style="italic", zorder=3,
                bbox=dict(boxstyle="square,pad=0.15", facecolor=BG, edgecolor="none"))
    yb = len(orden) + 0.3
    at.plot([0, 0.05], [yb, yb], color=INK, lw=1.0, solid_capstyle="butt", clip_on=False)
    at.text(0.025, yb + 0.3, "0.05 subst./site", ha="center", va="top", fontsize=6.0, clip_on=False)
    at.set_xlim(-0.004, XEND); at.axis("off"); al.set_xlim(0, 1); al.axis("off")
    for g, xc in zip(cols, xcol):       # B: census tiles
        for nn in orden:
            est = dat[nn][g]; yy = Y[nn]
            if est == "ausente":
                am.add_patch(Rectangle((xc - 0.44, yy - 0.44), 0.88, 0.88, facecolor=AUS, lw=0, zorder=2))
                am.plot(xc, yy, "x", ms=4.4, color=BG, mew=1.5, zorder=4)
            elif est == "presente":
                am.add_patch(Rectangle((xc - 0.44, yy - 0.44), 0.88, 0.88, facecolor=PRES_T, lw=0, zorder=2))
            else:
                am.add_patch(Rectangle((xc - 0.40, yy - 0.40), 0.80, 0.80, facecolor=BG, edgecolor=PARC_E, lw=0.8, zorder=2))
    am.add_patch(Rectangle((-0.5, -0.5), 1, len(orden), fill=False, ec=INK, lw=1.0, zorder=5))
    for nom, f, col in bandas:
        ys = [Y[nn] for nn in orden if f(nn)]
        am.add_patch(Rectangle((XM + 0.75, min(ys) - 0.44), 0.22, max(ys) - min(ys) + 0.88,
                               facecolor=AUS if col == BANDA_AUS else EXT, lw=0, clip_on=False))
        am.text(XM + 1.2, (min(ys) + max(ys)) / 2, nom, va="center", fontsize=6.5, clip_on=False,
                style="normal" if nom in ("outgroup", "Cobitidae") else "italic", weight="normal" if nom == "outgroup" else "bold")
    am.set_xlim(-0.7, XM + 0.7); am.set_ylim(len(orden) - 0.4, -0.6)
    am.set_xticks(xcol); am.set_xticklabels(cols, rotation=90, fontsize=6.3, style="italic")
    for tl, g in zip(am.get_xticklabels(), cols):
        if g == "DCX": tl.set_weight("bold")
        if g in GRUPOS_FAM[-1][1]: tl.set_color(INK2)
    am.xaxis.set_ticks_position("top"); am.tick_params(length=0); am.set_yticks([])
    for s_ in am.spines.values(): s_.set_visible(False)
    for nom, a, b in xg:
        ctrl = nom.startswith("synteny")
        am.plot([a - 0.3, b + 0.3], [len(orden) - 0.15] * 2, color=INK2 if ctrl else INK, lw=0.8, clip_on=False)
        am.text((a + b) / 2, len(orden) + 0.25, nom, ha="center", va="top", fontsize=6.0, clip_on=False,
                color=INK2 if ctrl else INK, style="normal")
    fig.text(0.01, 0.985, "A", fontsize=10, weight="bold", va="top"); fig.text(0.47, 0.985, "B", fontsize=10, weight="bold", va="top")
    hs = [Line2D([], [], marker="D", ms=4.2, color=PRES, ls="none", label="DCX present"),
          Line2D([], [], marker="D", ms=4.2, mfc=BG, mec=AUS, mew=1.4, ls="none", label="DCX absent"),
          Line2D([0, 1], [0, 0], color=AUS, lw=2.2, marker="X", ms=6.5, mec=BG, label="inferred loss"),
          Rectangle((0, 0), 1, 1, facecolor=BANDA_AUS, alpha=0.6, label="lineage without DCX"),
          Rectangle((0, 0), 1, 1, facecolor=BANDA_EXT, alpha=0.6, label="outgroup"),
          Line2D([], [], marker="o", color=SOP_COL, mec="none", ls="none", label="UFBoot 100", **sop_estilo(100)),
          Line2D([], [], marker="o", color=SOP_COL, mec="none", ls="none", label="UFBoot 90", **sop_estilo(90)),
          Line2D([], [], marker="o", color=SOP_COL, mec="none", ls="none", label="UFBoot ≤ 75", **sop_estilo(75)),
          Rectangle((0, 0), 1, 1, facecolor=AUS, label="gene absent  ×"),
          Rectangle((0, 0), 1, 1, facecolor=PRES_T, label="gene present"),
          Rectangle((0, 0), 1, 1, facecolor=BG, edgecolor=PARC_E, lw=0.8, label="partial detection")]
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.47, 0.065), ncol=3, fontsize=6.2, handlelength=1.5,
               columnspacing=1.2, frameon=False)
    nsin = {g: sum(1 for n in orden if dat[n][g] == "ausente") for g in cols}
    print(f"  fig1: losses = {n_perd}; absences per gene: { {g: v for g, v in nsin.items() if v} }")
    guardar(fig, "fig1_phylogeny_census")

# =============================================================== FIG 2: raw reads, observed / expected
LECTURAS = [  # species, DCX, dclk1a, dclk2a, dclk3 (reads per 100 aa), has_DCX
    ("Danionella cerebrum", 10.92, 9.87, 12.23, 14.42, True), ("Ictalurus furcatus", 14.57, 10.91, 13.56, 24.13, True),
    ("Danionella dracula", 0.28, 11.43, 9.08, 12.07, False), ("Danio tinwini", 0.56, 13.25, 14.53, 20.80, False),
    ("Paedocypris micromegethes", 0.56, 22.47, 13.68, 34.12, False),
    ("Paedocypris carbunculus", 1.12, 90.65, 68.28, 120.25, False),
    ("Misgurnus anguillicaudatus", 0.00, 24.16, 20.70, 28.85, False), ("Paramisgurnus dabryanus", 0.28, 5.97, 6.30, 6.80, False)]   # 1 ambiguous DC-domain read counted as DCX (conservative, same rule as all species)
LONG_AA = {"DCX": 357, "dclk1a": 770, "dclk2a": 826, "dclk3": 721}
def fig2():
    from scipy import stats
    gs_ = ["DCX", "dclk1a", "dclk2a", "dclk3"]; filas = []
    for nom, *tasas, tiene in LECTURAS:
        k = {g: round(tt * LONG_AA[g] / 100) for g, tt in zip(gs_, tasas)}
        E = sum(k[g] for g in gs_[1:]) / sum(LONG_AA[g] for g in gs_[1:]) * LONG_AA["DCX"]; kd = k["DCX"]
        lo = stats.chi2.ppf(0.025, 2 * kd) / 2 / E if kd > 0 else 0.0
        hi = stats.chi2.ppf(0.975, 2 * kd + 2) / 2 / E
        filas.append(dict(nom=nom, kd=kd, E=E, r=kd / E, lo=lo, hi=hi, p=stats.poisson.cdf(kd, E), tiene=tiene))
    PISO = 0.003
    fig, ax = plt.subplots(figsize=(5.0, 3.0)); fig.subplots_adjust(top=0.9, left=0.32, right=0.86, bottom=0.25)
    ax.axvspan(PISO * 0.7, 0.2, color=HUECO, alpha=0.35, lw=0, zorder=0)
    for i, f in enumerate(filas):
        ax.plot([max(f["lo"], PISO), f["hi"]], [i, i], color=PRES if f["tiene"] else AUS, lw=1.8, solid_capstyle="round", zorder=2)
        marca_estado(ax, max(f["r"], PISO), i, f["tiene"], ms=5.6, zorder=5)
        ax.text(2.1, i, f"{f['kd']} / {f['E']:.0f}", va="center", ha="left", fontsize=6.2, clip_on=False)
    ax.text(2.1, -0.85, "obs. / exp.\nreads", va="bottom", ha="left", fontsize=6.0, color=INK2, clip_on=False)
    ax.axvline(1, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xscale("log"); ax.set_xlim(PISO * 0.7, 1.9)
    ax.set_xticks([PISO, 0.01, 0.1, 1]); ax.set_xticklabels(["0", "0.01", "0.1", "1"])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_yticks(range(len(filas))); ax.set_yticklabels([f["nom"] for f in filas], style="italic", fontsize=6.8)
    ax.set_ylim(len(filas) - 0.4, -0.7)
    ax.set_xlabel("observed / expected DCX reads (95% CI)")
    ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); limpiar(ax); ax.tick_params(axis="y", length=0)
    hs = [Line2D([0, 1], [0, 0], color=PRES, lw=1.8, marker="D", ms=5.4, label="species with DCX"),
          Line2D([0, 1], [0, 0], color=AUS, lw=1.8, marker="D", ms=5.4, mfc=BG, mec=AUS, mew=1.4, label="species without DCX")]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=2, fontsize=6.2, handlelength=1.8, columnspacing=1.4)
    for f in filas: print(f"  fig2: {f['nom']:28s} {f['kd']:3d}/{f['E']:6.1f} = {f['r']:.3f} [{f['lo']:.3f}, {f['hi']:.3f}] p={f['p']:.1e}")
    guardar(fig, "fig2_raw_reads")

# =============================================================== FIG 3: microsynteny of the DCX-pak3loc locus
BLOQUE_GENES = ["kntc1", "hip1rb", "LOC307125267", "LOC307124699", "zfpl1", "rwdd2b", "spata22", "lhfpl7"]
BLOQUE_ETQ = {"LOC307125267": "DCX", "LOC307124699": "pak3loc"}
SP3 = ["Gobio_gobio", "Ctenopharyngodon_idella", "Pimephales_promelas", "Triplophysa_rosa", "Paedocypris_micromegethes",
       "Paedocypris_sp_Pulau_Singkep", "Danio_rerio", "Misgurnus_anguillicaudatus", "Paramisgurnus_dabryanus"]
GRUPO3 = {"Gobio_gobio": "Cyprinoidei", "Ctenopharyngodon_idella": "Cyprinoidei", "Pimephales_promelas": "Cyprinoidei",
          "Danio_rerio": "Danio", "Misgurnus_anguillicaudatus": "Cobitidae", "Paramisgurnus_dabryanus": "Cobitidae",
          "Triplophysa_rosa": "Nemacheilidae", "Paedocypris_micromegethes": "Paedocypris", "Paedocypris_sp_Pulau_Singkep": "Paedocypris"}
LIN3 = ["Danio", "Cobitidae", "Paedocypris"]
PERDIDOS_EXTRA = {("Danio_rerio", "LOC307124699")}   # pak3 absent from the complete D. rerio RefSeq annotation
VENTANA = 1_500_000
def fig3():
    DCXG = "LOC307125267"
    tabla = json.load(open(f"{D}/syn2/bloque/bloque_tabla.json"))
    t = Phylo.read(ARBOL, "newick")
    for x in list(t.get_terminals()):
        if x.name not in SP3: t.prune(x)
    t.root_at_midpoint(); t.ladderize(reverse=True)
    orden = [x.name for x in t.get_terminals()]; Y = {n: i for i, n in enumerate(orden)}
    def nucleo(sp):
        pos = [(g, v[0], v[1]) for g, v in tabla[sp].items() if g in BLOQUE_GENES and v]; mejor = set()
        for g0, c0, p0 in pos:
            grp = {g for g, c, p in pos if c == c0 and abs(p - p0) <= VENTANA / 2 + 1}
            if len(grp) > len(mejor) or (len(grp) == len(mejor) and "kntc1" in grp): mejor = grp
        return mejor
    NUC = {sp: nucleo(sp) for sp in orden}
    def estado(sp, g):
        v = tabla[sp].get(g)
        if v is None: return "perdido" if (g == DCXG or (sp, g) in PERDIDOS_EXTRA) else "nodet"
        return "bloque" if g in NUC[sp] else "fuera"
    def distancia(sp, g):
        v = tabla[sp][g]; c0 = tabla[sp][next(iter(NUC[sp]))][0]
        if v[0] != c0: return "other chr."
        return f"{min(abs(v[1] - tabla[sp][h][1]) for h in NUC[sp]) / 1e6:.1f} Mb"
    fig = plt.figure(figsize=(7.2, 6.0))
    gsp = fig.add_gridspec(2, 2, height_ratios=[1.3, 0.7], width_ratios=[1.55, 1.0], hspace=0.3, wspace=0.5,
                           left=0.14, right=0.98, top=0.9, bottom=0.03)
    ax = fig.add_subplot(gsp[0, 0]); axb = fig.add_subplot(gsp[0, 1]); axc = fig.add_subplot(gsp[1, :])
    def caja(x, y, est, g, w=0.72, h=0.5):
        dcx = g == DCXG; kw = dict(boxstyle="round,pad=0,rounding_size=0.08", zorder=4)
        if est == "bloque":
            c = PRES if dcx else SINT
            ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, facecolor=c, **borde(c), **kw))
        elif est == "fuera":
            ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, facecolor=BG, edgecolor=INK2, lw=1.0, **kw))
        elif est == "perdido":
            ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, facecolor=HUECO, edgecolor=AUS, lw=1.0, ls=(0, (2.2, 1.4)), **kw))
            ax.plot(x, y, "x", ms=4.5, color=AUS, mew=1.5, zorder=5)
        else:
            ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, fill=False, edgecolor=SINORT, lw=0.9, ls=(0, (2.2, 1.4)), **kw))
    n = len(BLOQUE_GENES)
    for a, b in zip(orden, orden[1:]):          # ribbons only between adjacent tracks of the same lineage
        if GRUPO3[a] != GRUPO3[b]: continue
        for j, g in enumerate(BLOQUE_GENES):
            if estado(a, g) == "bloque" and estado(b, g) == "bloque":
                dcx = g == DCXG; ya, yb = Y[a], Y[b]
                ax.add_patch(Polygon([(j - 0.25, ya + 0.25), (j + 0.25, ya + 0.25), (j + 0.25, yb - 0.25), (j - 0.25, yb - 0.25)],
                                     closed=True, facecolor=PRES if dcx else CINTA, alpha=0.55 if dcx else 0.45,
                                     edgecolor="none", zorder=1))
    for sp in orden:
        y = Y[sp]
        ax.plot([-0.45, n - 0.55], [y, y], color=CROM, lw=4.5, solid_capstyle="round", zorder=2)
        for j, g in enumerate(BLOQUE_GENES):
            e = estado(sp, g); caja(j, y, e, g)
            if e == "fuera": ax.text(j, y + 0.36, distancia(sp, g), ha="center", va="top", fontsize=5.2, color=INK2, zorder=5)
        ax.text(-0.75, y, bonito(sp), ha="right", va="center", fontsize=6.8, style="italic")
    for lin in LIN3:
        yy = [Y[s] for s in orden if GRUPO3[s] == lin]
        ax.add_patch(Rectangle((n - 0.35, min(yy) - 0.3), 0.14, max(yy) - min(yy) + 0.6, facecolor=AUS, lw=0, clip_on=False))
        ax.text(n - 0.1, (min(yy) + max(yy)) / 2, lin, va="center", fontsize=6.7, style="normal" if lin == "Cobitidae" else "italic", weight="bold", clip_on=False)
    ax.set_xticks(range(n)); ax.set_xticklabels([BLOQUE_ETQ.get(g, g) for g in BLOQUE_GENES], rotation=90, fontsize=6.6, style="italic")
    ax.get_xticklabels()[2].set_weight("bold"); ax.xaxis.set_ticks_position("top"); ax.tick_params(length=0)
    ax.set_yticks([]); [s_.set_visible(False) for s_ in ax.spines.values()]
    ax.set_xlim(-0.6, n - 0.4); ax.set_ylim(len(orden) - 0.4, -0.7)
    hs = [Rectangle((0, 0), 1, 1, facecolor=PRES, label="DCX in block"),
          Rectangle((0, 0), 1, 1, facecolor=SINT, **borde(SINT), label="neighbour in block"),
          Rectangle((0, 0), 1, 1, facecolor=BG, edgecolor=INK2, lw=1.0, label="outside block (distance)"),
          Rectangle((0, 0), 1, 1, facecolor=HUECO, edgecolor=AUS, ls=(0, (2, 1.3)), lw=1.0, label="lost  ×"),
          Rectangle((0, 0), 1, 1, facecolor=CINTA, alpha=0.5, label="conserved synteny")]
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.3, -0.02), ncol=3, fontsize=6.0, handlelength=1.2,
              handleheight=0.9, columnspacing=0.9)
    ax.text(-0.42, 1.2, "A", transform=ax.transAxes, fontsize=10, weight="bold")
    fig3_compactacion(axb); fig3_cortes(axc)
    guardar(fig, "fig3_microsynteny")

def fig3_compactacion(ax):
    """B: genome compaction control. Ratio (target / Gobio) of intergenic distances between consecutive uniquely
    anchored orthologues (miniprot, Gobio proteome) for Gobio intervals of 128-515 kb with no intervening gene,
    in P. micromegethes (no DCX) and Pimephales (DCX present); the hip1rb-pak3loc interval of each is marked."""
    import numpy as np, random
    rnd = random.Random(3)
    N = [l.rstrip().split("\t") for l in open(f"{D}/compactacion/nula_estratificada.tsv")][1:]
    res = {r["sp"]: r for r in json.load(open(f"{D}/compactacion/compactacion_results.json"))}
    LO, HI = -2.5, 1.0
    for y, sp, pres, nombre in [(1, "Pimephales_promelas", True, "P. promelas"), (0, "Paedocypris_micromegethes", False, "P. micromegethes")]:
        v = np.log10([float(x[4]) for x in N if x[0] == sp and x[1] == "0" and float(x[4]) > 0])
        ax.scatter(np.clip(v, LO, HI), [y + rnd.uniform(-0.17, 0.17) for _ in v], s=2.2, color="#A8A8A8", lw=0, zorder=2)
        q1, md, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([q1, q3], [y, y], color=INK, lw=3.4, solid_capstyle="butt", zorder=3)
        ax.plot([md], [y], "|", ms=10, mew=1.4, color=BG, zorder=4)
        ax.plot([np.log10(res[sp]["genoma_ratio"])] * 2, [y - 0.3, y + 0.3], color=INK2, lw=0.8, ls=(0, (2, 1.5)), zorder=3)
        r = res[sp]; x = np.log10(r["cociente"]); pct = 100 * np.mean(v <= x)
        marca_estado(ax, x, y, pres, ms=7, zorder=6)
        ax.text(LO - 0.12, y, nombre, ha="right", va="center", fontsize=6.6, style="italic")
        ax.text(x + (0.12 if not pres else 0.12), y + 0.27, it("hip1rb–pak3loc") + f" {r['gobio_kb']:.0f} → {r['diana_kb']:.1f} kb"
                f"\npercentile {pct:.0f}", ha="left" if not pres else "right", va="bottom", fontsize=5.8,
                transform=ax.transData) if False else ax.annotate(
                it("hip1rb–pak3loc") + f"\n{r['gobio_kb']:.0f} → {r['diana_kb']:.1f} kb (pct. {pct:.0f})", (x, y),
                xytext=(x + 0.08 if not pres else x - 0.08, y + 0.33), ha="left" if not pres else "right", va="bottom",
                fontsize=5.8, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    ax.set_xlim(LO - 0.05, HI); ax.set_ylim(-0.45, 1.85)
    ax.set_xticks([-2, -1, 0, 1]); ax.set_xticklabels(["0.01", "0.1", "1", "10"]); ax.set_yticks([])
    ax.set_xlabel("intergenic length, species / " + it("Gobio"))
    for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
    ax.text(-0.42, 1.04, "B", transform=ax.transAxes, fontsize=10, weight="bold")

def fig3_cortes(ax):
    """C: deletion in Paedocypris. Gobio NC_137100.1 (8.898-9.164 Mb) vs P. micromegethes CM051502.1 (12.545-12.562 Mb,
    drawn reversed); ribbons = dc-megablast HSPs (e < 1e-5). The same structure is found in P. sp. Singkep."""
    GL, PL = 266001, 17001
    gx = lambda p: p / GL; px = lambda p: 0.12 + 0.76 * (PL - p) / PL      # Paedocypris drawn reversed, own scale
    yg, yp = 1.0, 0.0
    ax.plot([0, 1], [yg, yg], color=INK2, lw=0.8); ax.plot([px(PL), px(1)], [yp, yp], color=INK2, lw=0.8)
    def gen(x0, x1, y, col, etq, dcx=False, ha="center"):
        ax.add_patch(Rectangle((x0, y - 0.09), max(x1 - x0, 0.004), 0.18, facecolor=col, edgecolor=INK2, lw=0.5, zorder=3))
        ax.text((x0 + x1) / 2, y + 0.16, etq, ha="center", va="bottom", fontsize=6.4, style="italic",
                weight="bold" if dcx else "normal")
    # Gobio genes (region coordinates)
    gen(gx(0), gx(5763), yg, SINT, "hip1rb"); gen(gx(31954), gx(71236), yg, SINT, "LOC307125817")
    gen(gx(176719), gx(248225), yg, PRES, "DCX", dcx=True); gen(gx(260768), gx(GL), yg, SINT, "pak3loc")
    ax.add_patch(Rectangle((gx(5763), yg - 0.2), gx(256631) - gx(5763), 0.4, facecolor=HUECO, alpha=0.5, lw=0, zorder=1))
    ax.text(gx(130000), yg - 0.26, "absent from " + it("Paedocypris") + " (≈251 kb of " + it("Gobio") + " sequence)", ha="center", va="top", fontsize=6.2,
            color=AUS)
    # Paedocypris genes
    gen(px(17001), px(11273), yp, SINT, "hip1rb"); gen(px(5541), px(1), yp, SINT, "pak3loc")
    # ribbons (Paedocypris start,end -> Gobio start,end)
    for (a, b, c, d) in [(5449, 5563, 262014, 262128), (7363, 7558, 257720, 257906), (7968, 8374, 256631, 257033),
                         (11271, 11501, 1921, 2175)]:
        ax.add_patch(Polygon([(px(a), yp + 0.09), (px(b), yp + 0.09), (gx(c), yg - 0.09), (gx(d), yg - 0.09)], closed=True,
                             facecolor=CINTA, alpha=0.75, edgecolor=INK2, lw=0.3, zorder=2))
    ax.text(px(8100), yp - 0.2, "conserved non-coding\nelements", ha="center", va="top", fontsize=5.8, color=INK2)
    ax.text(-0.01, yg, "G. gobio", ha="right", va="center", fontsize=6.8, style="italic")
    ax.text(px(PL) - 0.01, yp, "P. micromegethes", ha="right", va="center", fontsize=6.8, style="italic")
    ax.plot([gx(GL - 50000), gx(GL)], [yg + 0.42, yg + 0.42], color=INK, lw=0.9); ax.text(gx(GL - 25000), yg + 0.46, "50 kb",
            ha="center", va="bottom", fontsize=5.8)
    ax.plot([px(2000), px(0)], [yp - 0.42, yp - 0.42], color=INK, lw=0.9); ax.text(px(1000), yp - 0.46, "2 kb",
            ha="center", va="top", fontsize=5.8)
    ax.set_xlim(-0.12, 1.02); ax.set_ylim(-0.7, 1.6); ax.axis("off")
    ax.text(-0.12, 1.6, "C", fontsize=10, weight="bold", va="top")

# =============================================================== FIG S1: DCX across RefSeq vertebrates
VERIF = {"Danio rerio": "loss", "Danio aesculapii": "loss", "Misgurnus anguillicaudatus": "loss",
         "Paramisgurnus dabryanus": "loss", "Ictalurus furcatus": "gap"}
CLASES = [("Mammalia", "Mammals"), ("Aves", "Birds"), ("Lepidosauria", "Reptiles"), ("Testudines", "Reptiles"),
          ("Crocodylia", "Reptiles"), ("Amphibia", "Amphibians"), ("Cypriniformes", "Cypriniformes"),
          ("Actinopterygii", "Other ray-finned fishes"), ("Chondrichthyes", "Cartilaginous fishes"),
          ("Dipnomorpha", "Other lobe-finned fishes"), ("Coelacanthimorpha", "Other lobe-finned fishes"),
          ("Cyclostomata", "Cyclostomes")]
def marca_verif(ax, x, y, est, **kw):
    if est == "loss": ax.plot(x, y, "x", ms=5, color=AUS, mew=1.7, **kw)
    elif est == "gap": ax.plot(x, y, "o", ms=5, fillstyle="left", mfc=PARC, mfcalt=BG, mec=INK, mew=0.6, **kw)
    else: ax.plot(x, y, "o", ms=5, mfc=SINORT, mec=INK2, mew=0.8, **kw)
def figS1():
    L = json.load(open(f"{D}/p1/vertebrate_dcx_presence_lineages.json"))
    v = [r for r in csv.DictReader(open(f"{D}/p1/vertebrate_dcx_presence.tsv"), delimiter="\t") if r["accesion"].startswith("GCF")]
    def clase(tx):
        for k, nm in CLASES:
            if k in L[tx]: return nm
        return "other"
    tiene = lambda r: r["n_DCX"] != "0" or r["especie"] == "Xenopus laevis"   # dcx.L / dcx.S
    tot, sin = collections.Counter(), collections.defaultdict(list)
    for r in v:
        c = clase(r["taxid"]); tot[c] += 1
        if not tiene(r): sin[c].append(r["especie"])
    orden = ["Mammals", "Birds", "Reptiles", "Amphibians", "Cypriniformes", "Other ray-finned fishes",
             "Cartilaginous fishes", "Other lobe-finned fishes", "Cyclostomes"]
    fig = plt.figure(figsize=(6.6, 2.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.04, top=0.92, bottom=0.17)
    ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1], sharey=ax)
    for i, c in enumerate(orden):
        nn, k = tot[c], len(sin[c])
        ax.barh(i, nn - k, height=0.55, color=PRES, zorder=3)
        ax.text(nn - k + 4, i, f"{nn - k}/{nn}", va="center", fontsize=6.4)
        x, yy = 0.14, i
        for s in sin[c]:
            nom = bonito(s.replace(" ", "_")); ancho = 0.05 + 0.0165 * len(nom)
            if x + ancho > 1.02: x, yy = 0.14, yy + 0.45
            marca_verif(bx, x, yy, VERIF.get(s, "unverified"), clip_on=False)
            bx.text(x + 0.035, yy, nom, va="center", fontsize=6.2, style="italic")
            x += ancho + 0.03
    ax.set_yticks(range(len(orden))); ax.set_yticklabels(orden, fontsize=6.8)
    ax.set_ylim(len(orden) - 0.4, -0.6); ax.set_xlim(0, 300); ax.set_xlabel("genomes with annotated DCX")
    ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); limpiar(ax); ax.tick_params(axis="y", length=0)
    bx.set_xlim(0, 1); bx.axis("off"); bx.text(0.14, -0.95, "no DCX annotated", fontsize=6.6, color=INK2)
    hs = [Line2D([], [], marker="x", ms=5, color=AUS, mew=1.7, ls="none", label="verified loss"),
          Line2D([], [], marker="o", ms=5, fillstyle="left", mfc=PARC, mfcalt=BG, mec=INK, mew=0.6, ls="none",
                 label="assembly gap (DCX in reads)"),
          Line2D([], [], marker="o", ms=5, mfc=SINORT, mec=INK2, mew=0.8, ls="none", label="not verified")]
    bx.legend(handles=hs, loc="lower left", bbox_to_anchor=(0.1, 0.0), fontsize=6.2, labelspacing=0.5)
    nt = sum(tot.values()); ns = sum(len(x) for x in sin.values())
    print(f"  figS1: {nt} genomes, {nt - ns} with DCX; without: { {k: v for k, v in sin.items() if v} }")
    guardar(fig, "figS1_vertebrates")

# =============================================================== FIG S2: topology tests
def figS2():
    casos = [("arbol/test.iqtree", "Codons, GTR+F+G4"), ("arbol_ctrl/a_AA_test.iqtree", "Amino acids, LG+F+G4"),
             ("arbol_ctrl/b_NOOG_test.iqtree", "No outgroup"), ("arbol_ctrl/c_SLOW_test.iqtree", "Fastest 25% of genes removed")]
    datos = []
    for f, et in casos:
        txt = open(f"{D}/{f}").read().split("USER TREES")[1]
        c = [l.split() for l in txt.splitlines() if re.match(r"^\s+2\s+-\d", l)][0]
        datos.append((et, float(c[2]), float(c[-2])))
    fig, ax = plt.subplots(figsize=(5.0, 1.9)); fig.subplots_adjust(top=0.95, bottom=0.26)
    for i, (et, dl, pau) in enumerate(datos):
        x = -math.log10(pau)
        ax.plot([0, x], [i, i], color=GRID, lw=1.2, zorder=1); ax.plot(x, i, "o", ms=6, color=INK, zorder=3)
        ax.text(x + 2, i, f"p = {cient(pau)}", va="center", fontsize=6.3)
    ax.axvline(-math.log10(0.05), color=INK2, lw=0.9, ls=(0, (3, 2)))
    ax.text(-math.log10(0.05) + 0.8, -0.45, "p = 0.05", fontsize=6, color=INK2, va="center")
    ax.set_yticks(range(len(datos))); ax.set_yticklabels([d[0] for d in datos], fontsize=6.8)
    ax.set_ylim(len(datos) - 0.4, -0.6); ax.set_xlim(0, 135)
    ax.set_xlabel(r"$-\log_{10}$ AU p-value of the constrained topology")
    ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True); limpiar(ax); ax.tick_params(axis="y", length=0)
    print(f"  figS2: {datos}")
    guardar(fig, "figS2_topology_tests")

# =============================================================== FIG S3: RELAX null distribution
FAMILIA = ["dclk1a", "dclk1b", "dclk2a", "dclk2b", "dclk3", "dcdc2", "dcdc2b", "dcdc2c", "rp1l1a", "rp1l1b"]
VECINOS = ["hip1rb", "kntc1", "rwdd2b", "spata22", "zfpl1"]
ETIQUETAR = {"todas_las_perdidas": ["dclk2b", "dclk1a", "spata22", "dclk2a"],
             "solo_Danio": ["dclk1b", "dcdc2", "dclk1a", "dclk2a", "zfpl1"]}
def K(f):
    try: return float(json.load(open(f))["test results"]["relaxation or intensification parameter"])
    except Exception: return None
def figS3():
    import numpy as np
    fig, axs = plt.subplots(2, 1, figsize=(5.4, 4.3), sharex=True, gridspec_kw=dict(hspace=0.3))
    fig.subplots_adjust(top=0.95, bottom=0.11)
    lo, hi = -2.0, 1.5
    for ax, (conj, tit) in zip(axs, [("todas_las_perdidas", "A  Test branches: all four losses"),
                                      ("solo_Danio", "B  Test branches: " + it("Danio") + " only")]):
        nula = [k for k in (K(f) for f in glob.glob(f"{D}/p5/relax/R*.{conj}.RELAX.json")) if k is not None]
        lk = np.clip(np.log10(np.maximum(nula, 1e-9)), lo, hi - 1e-6)
        cnt, _, _ = ax.hist(lk, bins=np.linspace(lo, hi, 29), color=SINORT, edgecolor=BG, lw=0.6, zorder=2)
        ymax = max(cnt) * 1.12; ybot = -0.62 * ymax
        ax.add_patch(Rectangle((math.log10(5), 0), hi - math.log10(5), ymax, color=CROM, lw=0, zorder=0))
        ax.text(math.log10(5) + 0.04, ymax * 0.97, "K > 5\nunstable", fontsize=5.6, color=INK2, va="top")
        ax.axhline(0, color=INK2, lw=0.6, zorder=2); ax.plot([0, 0], [0, ymax], color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=1)
        ks = {}
        for genes, col in [(FAMILIA, INK), (VECINOS, SINT)]:
            for g in genes:
                k = K(f"{P}/relax/{g}.{conj}.RELAX.json")
                if k is None: continue
                x = min(max(math.log10(max(k, 1e-9)), lo), hi - 0.02); ks[g] = (x, k, col)
                if col == SINT: ax.plot([x, x], [-0.03 * ymax, -0.16 * ymax], color=INK2, lw=2.6, zorder=3, solid_capstyle="butt")
                ax.plot([x, x], [-0.03 * ymax, -0.16 * ymax], color=col, lw=1.4, zorder=3, solid_capstyle="butt")
        et = sorted(ETIQUETAR[conj], key=lambda g: ks[g][0]); sep = 0.42
        tx = [ks[g][0] for g in et]
        for j in range(1, len(tx)): tx[j] = max(tx[j], tx[j - 1] + sep)
        sh = (sum(ks[g][0] for g in et) - sum(tx)) / len(tx)
        tx = [min(max(tt + sh, lo + 0.2), hi - 0.2) for tt in tx]
        for j in range(1, len(tx)): tx[j] = max(tx[j], tx[j - 1] + sep)
        for g, tt in zip(et, tx):
            x, k, col = ks[g]; pct = 100 * sum(1 for z in nula if z >= k) / len(nula)
            ax.plot([x, tt], [-0.17 * ymax, -0.27 * ymax], color=INK2, lw=0.5, zorder=3)
            ax.text(tt, -0.29 * ymax, f"{g}\n{pct:.0f}%", ha="center", va="top", fontsize=5.7, style="italic",
                    weight="bold" if col == INK else "normal")
        ax.set_ylim(ybot, ymax); ax.set_xlim(lo, hi); ax.set_ylabel("random genes", fontsize=6.8)
        ax.set_title(f"{tit} (n = {len(nula)})", loc="left", fontsize=7.2, pad=3)
        limpiar(ax, left=False); ax.spines["left"].set_bounds(0, ymax)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=4))
        ax.set_yticks([tt for tt in ax.get_yticks() if 0 <= tt <= ymax])
        if ax is axs[0]: ax.spines["bottom"].set_visible(False)
        ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    axs[1].set_xticks([-2, -1, 0, 1]); axs[1].set_xticklabels(["0.01", "0.1", "1", "10"]); axs[0].tick_params(axis="x", length=0)
    axs[1].set_xlabel("RELAX K (log scale; K > 1 intensification, K < 1 relaxation)")
    hs = [Rectangle((0, 0), 1, 1, facecolor=SINORT, label="null: random genes"),
          Line2D([], [], color=INK, lw=1.4, label="DCX/DCLK family (bold label)"),
          Line2D([], [], color=SINT, lw=2.2, path_effects=[pe.Stroke(linewidth=3.4, foreground=INK2), pe.Normal()],
                 label="synteny controls")]
    axs[0].legend(handles=hs, loc="upper left", fontsize=6.2, handlelength=1.2)
    guardar(fig, "figS3_relax_null")

# =============================================================== FIG 4: DC-only (DCL) isoform of dclk1a
# Domain positions: human DCLK1 (UniProt O15075): doublecortin 1 = 54-140, doublecortin 2 = 180-263, kinase = 389-651.
DOM_DCL = [("DC1", 54, 140), ("DC2", 180, 263)]
def fig4():
    """A: dclk1a gene model in P. micromegethes (no DCX) with the two isoforms seen in Iso-Seq; protein cartoons.
    B: the DCL-specific C-terminal tail across Otophysi, with DCX status (diamonds as in Figs. 1-2) and evidence."""
    iso = json.load(open(f"{D}/rnaseq_dcl/isoseq_conteos.json"))
    filas = [l.rstrip("\n").split("\t") for l in open(f"{D}/isoformas/dcl_colas.tsv") if not l.startswith("#")][1:]
    x0 = iso["exones_compartidos"][0][0]; kb = lambda p: (p - x0) / 1000
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 2.25], width_ratios=[2.3, 1.0], hspace=0.02, wspace=0.06,
                          left=0.02, right=0.99, top=0.97, bottom=0.03)
    a = fig.add_subplot(gs[0, 0]); ap = fig.add_subplot(gs[0, 1]); b = fig.add_subplot(gs[1, :])
    # ---- A: gene model and isoforms
    COMP, KIN, DCLC = INK2, CINTA, PARC
    def caja(ax, s, e, y, col, h=0.34, **kw):
        w = max(kb(e) - kb(s), 0.05)
        ax.add_patch(Rectangle((kb(s), y - h / 2), w, h, facecolor=col, lw=0.5, edgecolor=INK2 if col != COMP else COMP,
                               zorder=3, **kw))
    def exon_dcl(y):   # DCL terminal exon: short coding part (tall) + long 3' UTR (half height)
        caja(a, iso["exon_dcl"][0], iso["cds_dcl_fin"], y, DCLC)
        caja(a, iso["cds_dcl_fin"], iso["exon_dcl"][1], y, DCLC, h=0.17)
    def isoforma(y, exs, etiqueta, n):
        for (s1, e1), (s2, e2) in zip(exs, exs[1:]):
            m = (kb(e1) + kb(s2)) / 2
            a.plot([kb(e1), m, kb(s2)], [y + 0.02, y + 0.2, y + 0.02], color=INK2, lw=0.6, zorder=2)
        for s, e, col in [(s, e, c) for (s, e), c in zip(exs, cols(exs))]:
            exon_dcl(y) if col == DCLC else caja(a, s, e, y, col)
        a.text(-0.4, y, etiqueta, ha="right", va="center", fontsize=6.8)
        a.text(kb(exs[-1][1]) + 0.35, y, f"n = {n}", ha="left", va="center", fontsize=6.5, color=INK2)
    def cols(exs):
        return [DCLC if (s, e) == tuple(iso["exon_dcl"]) else (KIN if [s, e] in iso["exones_cinasa"] else COMP) for s, e in exs]
    comp = [tuple(e) for e in iso["exones_compartidos"]]; dcl = tuple(iso["exon_dcl"]); kin = [tuple(e) for e in iso["exones_cinasa"]]
    y_gen, y_l, y_d = 2.0, 1.0, 0.0
    a.plot([kb(x0), kb(kin[-1][1])], [y_gen, y_gen], color=INK2, lw=0.8, zorder=1)
    for s, e in comp: caja(a, s, e, y_gen, COMP)
    exon_dcl(y_gen); [caja(a, s, e, y_gen, KIN) for s, e in kin]
    a.text(-0.4, y_gen, it("dclk1a") + " locus", ha="right", va="center", fontsize=6.8)
    isoforma(y_l, comp + kin, "long (DC + kinase)", iso["largo"])
    isoforma(y_d, comp + [dcl], "DCL (DC only)", iso["dcl"])
    # stop codon of the DCL exon
    a.plot([kb(iso["cds_dcl_fin"])], [y_d + 0.3], marker="v", ms=3.6, color=INK, zorder=5)
    a.text(kb(iso["cds_dcl_fin"]), y_d + 0.42, "stop", ha="center", va="bottom", fontsize=5.6)
    a.plot([kb(x0), kb(x0) + 5], [-0.75, -0.75], color=INK, lw=0.9, solid_capstyle="butt")
    a.text(kb(x0) + 2.5, -0.85, "5 kb", ha="center", va="top", fontsize=6)
    a.set_xlim(-6.8, kb(kin[-1][1]) + 2.2); a.set_ylim(-1.0, 2.5); a.axis("off")
    a.text(-6.8, 2.45, "A", fontsize=10, weight="bold", va="top")
    # ---- A (right): protein cartoons
    def proteina(y, L, kin_on, etq):
        ap.add_patch(Rectangle((0, y - 0.07), L, 0.14, facecolor=GRID, lw=0, zorder=1))
        for nm, s, e in DOM_DCL:
            ap.add_patch(Rectangle((s, y - 0.22), e - s, 0.44, facecolor=PRES_T, edgecolor=INK2, lw=0.5, zorder=2))
            ap.text((s + e) / 2, y, nm, ha="center", va="center", fontsize=5.6)
        if kin_on:
            ap.add_patch(Rectangle((389, y - 0.22), 651 - 389, 0.44, facecolor=KIN, edgecolor=INK2, lw=0.5, zorder=2))
            ap.text(520, y, "kinase", ha="center", va="center", fontsize=5.6)
        else:
            ap.add_patch(Rectangle((L - 11, y - 0.22), 11, 0.44, facecolor=DCLC, edgecolor=INK2, lw=0.5, zorder=2))
        ap.text(L + 12, y, etq, ha="left", va="center", fontsize=6.3, color=INK2)
    proteina(1.0, 740, True, "~740 aa"); proteina(0.0, 346, False, "346 aa")
    ap.set_xlim(-10, 900); ap.set_ylim(-1.0, 2.5); ap.axis("off")
    # ---- B: DCL tail across Otophysi
    cons = "EMDKAIIDSTA"; n = len(filas)
    xs_cola = 4.6; dx = 0.36
    grupos = []
    for i, (g, sp, dcx, cola, ev) in enumerate(filas):
        y = i
        if not grupos or grupos[-1][0] != g: grupos.append([g, i, i])
        else: grupos[-1][2] = i
        marca_estado(b, 0.35, y, dcx == "1", ms=4.6)
        b.text(0.75, y, bonito(sp), ha="left", va="center", fontsize=6.8, style="italic")
        for j, aa in enumerate(cola):
            igual = j < len(cons) and aa == cons[j]
            b.text(xs_cola + j * dx, y, aa, ha="center", va="center", fontsize=7.2, family="DejaVu Sans Mono",
                   color="#A0A0A0" if igual else INK, weight="normal" if igual else "bold")
        b.text(xs_cola + len(cola) * dx, y, "*", ha="center", va="center", fontsize=7.2, family="DejaVu Sans Mono")
        b.text(xs_cola + 12.6 * dx, y, ev, ha="left", va="center", fontsize=6.4, color=INK2)
    for g, i0, i1 in grupos:
        sin_dcx = all(f[2] == "0" for f in filas[i0:i1 + 1])
        b.add_patch(Rectangle((-0.05, i0 - 0.45), 0.12, i1 - i0 + 0.9, facecolor=AUS if sin_dcx else INK2, lw=0, clip_on=False))
        b.text(-0.15, (i0 + i1) / 2, g, ha="right", va="center", fontsize=6.4, style="italic" if g == "Paedocypris" else "normal",
               weight="bold" if sin_dcx else "normal")
        if i1 < n - 1: b.plot([-1.8, xs_cola + 17 * dx], [i1 + 0.5, i1 + 0.5], color=GRID, lw=0.6, zorder=0)
    b.text(xs_cola + 5 * dx, -1.1, "DCL-specific C-terminus", ha="center", va="bottom", fontsize=6.6)
    b.text(xs_cola + 12.6 * dx, -1.1, "evidence", ha="left", va="bottom", fontsize=6.6)
    b.set_xlim(-3.6, xs_cola + 19 * dx); b.set_ylim(n - 0.4, -1.6); b.axis("off")
    b.text(-3.6, -1.6, "B", fontsize=10, weight="bold", va="top")
    hs = [Rectangle((0, 0), 1, 1, facecolor=COMP, label="shared exons (DC domains)"),
          Rectangle((0, 0), 1, 1, facecolor=DCLC, edgecolor=INK2, lw=0.5, label="DCL terminal exon (CDS tall, 3' UTR short)"),
          Rectangle((0, 0), 1, 1, facecolor=KIN, edgecolor=INK2, lw=0.5, label="kinase exons"),
          Line2D([], [], marker="D", ms=4.2, color=PRES, ls="none", label="DCX present"),
          Line2D([], [], marker="D", ms=4.2, mfc=BG, mec=AUS, mew=1.3, ls="none", label="DCX absent"),
          Line2D([], [], marker="$\\mathbf{V}$", ms=5, color=INK, ls="none", label="residue differs from consensus (EMDKAIIDSTA)")]
    b.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=6.2, handlelength=1.3,
             columnspacing=1.2)
    guardar(fig, "fig4_dcl_isoform")

# =============================================================== FIG 5: selective constraint on DCX where it is retained
FAM5 = {"dclk1b", "dclk2a", "dclk3", "dcdc2", "dcdc2b", "dcdc2c", "rp1l1a"}
CTRL5 = {"hip1rb", "kntc1", "rwdd2b", "spata22", "zfpl1"}
def fig5():
    """A: omega (FEL global MG94xREV) vs fraction of codons under purifying selection (FEL, p < 0.1) in the nine species
    that retain DCX, for DCX, family members, synteny controls and 46 random genes (same species and pruned species tree).
    B: RELAX K for the Danionidae that retain DCX (D. cerebrum + Rasbora) vs the same random genes."""
    R = json.load(open(f"{D}/restriccion/restriccion_results.json"))
    az = [r for r in R if r["gen"].startswith("R") and "omega" in r]
    fam = [r for r in R if r["gen"] in FAM5 and "omega" in r]; ctl = [r for r in R if r["gen"] in CTRL5 and "omega" in r]
    dcx = next(r for r in R if r["gen"] == "DCX")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.32))
    OMIN = 0.005; om = lambda v: max(v, OMIN)
    a.scatter([om(r["omega"]) for r in az], [100 * r["purif"] for r in az], s=11, color="#9A9A9A", lw=0, zorder=2,
              label=f"random genes (n = {len(az)})")
    a.scatter([om(r["omega"]) for r in ctl], [100 * r["purif"] for r in ctl], s=20, marker="s", color=SINT, edgecolor=INK2,
              lw=0.5, zorder=3, label="synteny controls (" + it("hip1rb, kntc1,") + "\n" + it("rwdd2b, spata22, zfpl1") + ")")
    a.scatter([om(r["omega"]) for r in fam], [100 * r["purif"] for r in fam], s=20, marker="s", facecolor=BG, edgecolor=INK,
              lw=0.8, zorder=3, label="DCX/DCLK family")
    a.plot([dcx["omega"]], [100 * dcx["purif"]], "D", ms=7, color=PRES, mec=INK, mew=0.6, zorder=5, label="DCX")
    a.annotate("DCX", (dcx["omega"], 100 * dcx["purif"]), xytext=(-4, 6), textcoords="offset points", ha="right",
               fontsize=7.5, weight="bold")
    OFF = {"dclk1b": (-5, 4, "right"), "dclk2a": (-5, -7, "right"), "dclk3": (2, -9, "center"), "dcdc2": (6, 3, "left"),
           "dcdc2b": (6, -4, "left"), "rp1l1a": (6, 0, "left"), "dcdc2c": (6, -6, "left")}
    for r in fam:              # synteny controls are named in the legend (they overlap in the centre of the cloud)
        dx_, dy_, ha_ = OFF.get(r["gen"], (4, 0, "left"))
        a.annotate(r["gen"], (om(r["omega"]), 100 * r["purif"]), xytext=(dx_, dy_), textcoords="offset points", fontsize=5.6,
                   color=INK2, style="italic", ha=ha_, va="center")
    mo = sorted(r["omega"] for r in az)[len(az) // 2]; mp = sorted(100 * r["purif"] for r in az)[len(az) // 2]
    a.axvline(mo, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=1); a.axhline(mp, color=INK2, lw=0.6, ls=(0, (3, 2)), zorder=1)
    a.set_xscale("log"); a.set_xlim(OMIN * 0.8, 1.0); a.set_ylim(0, 80)
    a.set_xticks([0.01, 0.1, 1]); a.set_xticklabels(["0.01", "0.1", "1"])
    a.set_xlabel("dN/dS (ω), species that retain DCX"); a.set_ylabel("codons under purifying selection (%)")
    pct_o = 100 * sum(r["omega"] <= dcx["omega"] for r in az) / len(az)
    pct_p = 100 * sum(r["purif"] <= dcx["purif"] for r in az) / len(az)
    a.text(0.98, 0.97, f"DCX: ω = {dcx['omega']:.3f} (percentile {pct_o:.0f})\n{100*dcx['purif']:.0f}% purifying codons "
           f"(percentile {pct_p:.0f})", transform=a.transAxes, ha="right", va="top", fontsize=6.3,
           bbox=dict(boxstyle="square,pad=0.25", facecolor=BG, edgecolor="none"), zorder=6)
    a.legend(loc="lower left", fontsize=5.9, handletextpad=0.3, borderaxespad=0.3)
    limpiar(a, left=False); a.text(-0.16, 1.02, "A", transform=a.transAxes, fontsize=10, weight="bold")
    # ---- B: RELAX K
    KMIN, KMAX = 0.25, 4.0
    import random; rnd = random.Random(7)
    ks = [r["K"] for r in az if "K" in r]
    yj = [rnd.uniform(-0.22, 0.22) for _ in ks]
    b.scatter([min(max(k, KMIN), KMAX) for k in ks], yj, s=11, color="#9A9A9A", lw=0, zorder=2)
    nb = sum(k < KMIN for k in ks)
    if nb: b.text(KMIN * 1.02, 0.26, f"{nb} at K≈0", fontsize=5.6, color=INK2, va="bottom")
    b.plot([dcx["K"]], [0], "D", ms=7, color=PRES, mec=INK, mew=0.6, zorder=5)
    b.annotate(f"DCX\nK = {dcx['K']:.2f}, p = {dcx['pK']:.2f}", (dcx["K"], -0.03), xytext=(dcx["K"], -0.62),
               textcoords="data", ha="center", va="top", fontsize=6.4, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    b.axvline(1, color=INK, lw=0.7, zorder=1)
    b.set_xscale("log", base=2); b.set_xlim(KMIN * 0.9, KMAX * 1.1); b.set_ylim(-1.1, 0.75)
    b.set_xticks([0.25, 0.5, 1, 2, 4]); b.set_xticklabels(["0.25", "0.5", "1", "2", "4"]); b.set_yticks([])
    b.set_xlabel("RELAX K, Danionidae that retain DCX")
    b.text(0.5, 0.62, "← relaxation", ha="right", fontsize=6.2, color=INK2, transform=b.get_xaxis_transform() if False else b.transData)
    b.text(2.0, 0.62, "intensification →", ha="left", fontsize=6.2, color=INK2)
    pk = 100 * sum(k <= dcx["K"] for k in ks) / len(ks)
    b.text(0.98, 0.03, f"{pk:.0f}% of random genes\nhave K ≤ DCX", transform=b.transAxes, ha="right", va="bottom", fontsize=6.2)
    for s_ in ("top", "right", "left"): b.spines[s_].set_visible(False)
    b.text(-0.06, 1.02, "B", transform=b.transAxes, fontsize=10, weight="bold")
    guardar(fig, "fig5_constraint")

def figS4():
    """A: expression (RPKM) of DCX, the kinase region of dclk1a (long isoform) and the DC-only terminal exon of dclk1a in six
    adult tissues of C. idella (PRJNA1306682; two samples each, first 10 M read 1). B: reads joining the last shared exon of
    dclk1a to the kinase exons or to the DC-only exon (short-read junctions; Iso-Seq = full-length transcripts)."""
    import pandas as pd
    T = pd.read_csv(f"{D}/rnaseq/cid_expr/expresion_tejidos.tsv", sep="\t")
    ORD = [("cerebro", "brain"), ("intestino", "intestine"), ("branquias", "gills"), ("rinon", "kidney"), ("higado", "liver"), ("bazo", "spleen")]
    SER = [("RPKM_DCX", "DCX", "D", PRES, PRES, 5.2), ("RPKM_kinasa", it("dclk1a") + ", kinase region (long isoform)", "s", INK2, INK2, 4.6),
           ("RPKM_DCL", it("dclk1a") + ", DC-only exon", "^", BG, INK, 5.0)]
    FLO = 0.05
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw=dict(width_ratios=[1.3, 1.0], wspace=0.75))
    for i, (t, lab) in enumerate(ORD):
        sub = T[T.tejido == t]
        for j, (col, nom, mk, fc, ec, ms) in enumerate(SER):
            x = i + (j - 1) * 0.22; v = [max(z, FLO) for z in sub[col]]
            a.plot([x] * len(v), v, mk, ms=ms * 0.75, mfc=fc, mec=ec, mew=0.8, alpha=0.9, zorder=3, label=nom if i == 0 else None)
            a.plot([x - 0.08, x + 0.08], [math.exp(sum(math.log(z) for z in v) / len(v))] * 2, color=ec, lw=1.0, zorder=4)
    a.set_yscale("log"); a.set_ylim(FLO * 0.7, 2000); a.axhline(FLO, color=GRID, lw=3, zorder=1)
    a.text(-0.45, FLO * 1.15, "not detected", fontsize=5.6, color=INK2, ha="left", va="bottom")
    a.set_xticks(range(len(ORD))); a.set_xticklabels([l for _, l in ORD]); a.set_xlim(-0.5, len(ORD) - 0.5)
    a.set_ylabel("RPKM (" + it("Ctenopharyngodon idella") + ")")
    a.legend(loc="upper right", ncol=1, fontsize=5.8, handletextpad=0.3, borderaxespad=0.1)
    limpiar(a, left=False); a.text(-0.13, 1.02, "A", transform=a.transAxes, fontsize=10, weight="bold")
    cb = T[T.tejido == "cerebro"]; iso = json.load(open(f"{D}/rnaseq_dcl/isoseq_conteos.json"))
    G = [(it("C. idella") + " brain, sample 1", int(cb.union_larga.iloc[0]), int(cb.union_DCL.iloc[0]), "DCX present"),
         (it("C. idella") + " brain, sample 2", int(cb.union_larga.iloc[1]), int(cb.union_DCL.iloc[1]), "DCX present"),
         (it("D. rerio") + " brain", 67, 0, "DCX absent"),
         (it("P. micromegethes") + " whole body\n(Iso-Seq)", iso["largo"], iso["dcl"], "DCX absent")]
    for k, (n, lg, dc, est) in enumerate(G):
        y = len(G) - 1 - k
        t = lg + dc; fl, fd = 100 * lg / t, 100 * dc / t
        b.barh(y, fl, height=0.55, color=INK2, lw=0, zorder=2, label="to kinase exons" if k == 0 else None)
        b.barh(y, fd, left=fl, height=0.55, color=BG, edgecolor=INK, hatch="////", lw=0.6, zorder=2, label="to DC-only exon" if k == 0 else None)
        b.text(101.5, y, f"{dc}/{t}", va="center", fontsize=6.2)
        marca_estado(b, 124, y, est == "DCX present", ms=4.2, clip_on=False)
    b.set_yticks(range(len(G))); b.set_yticklabels([g[0] for g in G][::-1], fontsize=6.4)
    b.tick_params(axis="y", pad=3, length=0); b.set_xlim(0, 100); b.set_ylim(-0.5, len(G) + 0.3)
    b.set_xlabel("junction reads from last shared exon (%)")
    b.legend(loc="upper right", fontsize=5.9, handletextpad=0.3, borderaxespad=0.1)
    b.plot([], [], "D", ms=4, color=PRES, label="DCX present"); b.plot([], [], "D", ms=4, mfc=BG, mec=AUS, mew=1.2, label="DCX absent")
    h, l = b.get_legend_handles_labels(); b.legend(h, l, loc="upper right", fontsize=5.6, handletextpad=0.3, borderaxespad=0.0, ncol=2, columnspacing=0.8)
    limpiar(b); b.text(-0.78, 1.02, "B", transform=b.transAxes, fontsize=10, weight="bold")
    guardar(fig, "figS4_expression")

# =============================================================== FIG 6: census on the dated tree (Sudasinghe et al. 2026)
def fig6():
    import numpy as np
    """DCX status of the read-based census and genomes on the time-calibrated tree of Sudasinghe et al. (2026), pruned to
    the census (scripts/exportar_arbol_fig6.R). Branches inside clades without DCX in magenta; numbered bars = the four
    losses (Dollo), with the stem interval in Ma; outer rings: DCX status and miniaturization (Sudasinghe et al. Table S5)."""
    import pandas as pd
    tr = Phylo.read(f"{D}/integracion/arbol_censo.tre", "newick")
    Tt = pd.read_csv(f"{D}/integracion/arbol_censo_puntas.tsv", sep="\t").set_index("punta")
    perd = pd.read_csv(f"{D}/integracion/perdidas_dollo.tsv", sep="\t")
    tips = tr.get_terminals(); N = len(tips)
    dep = tr.depths(); H = max(dep[t] for t in tips)
    GAP = 0.055; ang = {t: math.pi / 2 - (2 * math.pi * (1 - GAP)) * (i + 0.5) / N - math.pi * GAP for i, t in enumerate(tips)}
    for c in tr.find_clades(order="postorder"):
        if not c.is_terminal(): ang[c] = (ang[c.clades[0]] + ang[c.clades[-1]]) / 2
    aus = {}
    for c in tr.find_clades(order="postorder"):
        aus[c] = (Tt.estado.get(c.name) == "ausente") if c.is_terminal() else all(aus[h] for h in c.clades)
    fig = plt.figure(figsize=(7.2, 7.2)); a = fig.add_axes([0, 0, 1, 1]); a.set_aspect("equal"); a.axis("off")
    xy = lambda r, t: (r * math.cos(t), r * math.sin(t))
    # time rings (every 25 Ma)
    for ma in range(25, int(H) + 1, 25):
        r = H - ma; tt = np.linspace(0, 2 * math.pi, 400)
        a.plot(r * np.cos(tt), r * np.sin(tt), color=GRID, lw=0.5, zorder=0)
        x, y = xy(r, math.pi / 2); a.text(x, y + 1.2, f"{ma}", ha="center", va="bottom", fontsize=5.2, color=INK2)
    a.text(0, H + 3.5, "Ma", ha="center", va="bottom", fontsize=5.2, color=INK2)
    for c in tr.find_clades():
        col, lw = (AUS, 1.1) if aus[c] else (INK2, 0.45)
        if not c.is_terminal():
            t0, t1 = ang[c.clades[-1]], ang[c.clades[0]]; tt = np.linspace(t0, t1, max(3, int(abs(t1 - t0) * 60)))
            a.plot(dep[c] * np.cos(tt), dep[c] * np.sin(tt), color=INK2, lw=0.45, solid_capstyle="butt", zorder=1)
            for h in c.clades:
                ch, lh = (AUS, 1.1) if aus[h] else (INK2, 0.45)
                (x0, y0), (x1, y1) = xy(dep[c], ang[h]), xy(dep[h], ang[h])
                a.plot([x0, x1], [y0, y1], color=ch, lw=lh, solid_capstyle="butt", zorder=2 if aus[h] else 1)
                if aus[h] and not aus[c]:
                    corto = lambda n: n.split("_", 2)[2]
                    k = [i for i, r_ in perd.iterrows() if corto(h.get_terminals()[0].name) in r_.puntas.split(", ")]
                    rm = (dep[c] + dep[h]) / 2; xm, ym = xy(rm, ang[h])
                    a.plot([xm], [ym], "o", ms=7.5, mfc=BG, mec=AUS, mew=1.2, zorder=5)
                    if k: a.text(xm, ym, str(k[0] + 1), ha="center", va="center", fontsize=5.6, weight="bold", color=INK, zorder=6)
                    h._perdida = k[0] + 1 if k else None
    # rings: DCX status and miniaturization
    R1, R2 = H + 4.5, H + 9.5
    MINI = {"progenetic-miniature": INK, "proportioned-dwarf": SINORT}
    dth = 2 * math.pi * (1 - GAP) / N
    for t in tips:
        e = Tt.loc[t.name]; th = ang[t]
        x, y = xy(R1, th)
        if e.estado == "presente": a.plot([x], [y], "D", ms=2.4, color=PRES, mew=0, zorder=4)
        else: a.plot([x], [y], "D", ms=2.6, mfc=BG, mec=AUS, mew=0.8, zorder=4)
        if e.miniatura in MINI:
            w = np.degrees(dth)
            a.add_patch(matplotlib.patches.Wedge((0, 0), R2 + 2.2, np.degrees(th) - w / 2, np.degrees(th) + w / 2, width=2.2,
                        facecolor=MINI[e.miniatura], edgecolor="none", zorder=3))
    # family arcs (families with >= 4 tips), label from tip name: Suborder_Family_Genus_...
    INC = {"Danionella": "Danionidae", "Thryssocypris": "Danionidae", "Paedocypris": "Paedocyprididae",
           "Sundadanio": "Sundadanionidae", "Fangfangia": "Sundadanionidae"}   # incertae sedis in the tip labels
    fam = [INC.get(t.name.split("_")[2], t.name.split("_")[1]) if "incertae" in t.name else t.name.split("_")[1] for t in tips]
    R3 = R2 + 4.5
    runs = []; i = 0
    while i < N:
        j = i
        while j + 1 < N and fam[j + 1] == fam[i]: j += 1
        runs.append([fam[i], i, j]); i = j + 1
    # merge runs of the same family interrupted by <= 6 tips of other families (one arc and one label per family)
    m = []
    for r in runs:
        if m and r[0] == m[-1][0]: m[-1][2] = r[2]; continue
        hit = [q for q in m[-3:] if q[0] == r[0] and r[1] - q[2] <= 7]
        if hit: q = hit[-1]; q[2] = r[2]; continue
        m.append(r)
    for f_, i, j in m:
        if j - i + 1 < 3: continue
        t0, t1 = ang[tips[j]] + dth * 0.3, ang[tips[i]] - dth * 0.3; tt = np.linspace(t0, t1, 80)
        a.plot(R3 * np.cos(tt), R3 * np.sin(tt), color=INK2, lw=0.9, solid_capstyle="butt")
        tm = (t0 + t1) / 2; x, y = xy(R3 + 2.0, tm); rot = np.degrees(tm)
        if math.cos(tm) < 0: rot += 180
        a.text(x, y, f_, rotation=rot, rotation_mode="anchor", fontsize=5.6, color=INK,
               ha="left" if math.cos(tm) >= 0 else "right", va="center")
    # loss table (center-left inset) and legend
    NOM = {1: None}
    txt = []
    for k, r_ in perd.iterrows():
        g = sorted({p.split("_")[0] for p in r_.puntas.split(", ")})
        nombre = {"Danionella": "D. dracula", "Danio": "Danio", "Paedocypris": "Paedocypris"}.get(g[0], "Cobitis + Misgurnus") if len(g) == 1 else "Cobitis + Misgurnus"
        lim = f"< {r_.edad_max_Ma:g} Ma" if r_.edad_min_Ma == 0 else f"{r_.edad_max_Ma:g}–{r_.edad_min_Ma:g} Ma"
        txt.append((k + 1, nombre, lim))
    y0 = -H - 16
    lx = -H - 22
    hdl = [Line2D([], [], marker="D", ls="", ms=4, color=PRES, mew=0, label="DCX present"),
           Line2D([], [], marker="D", ls="", ms=4, mfc=BG, mec=AUS, mew=1.0, label="DCX absent"),
           Line2D([], [], color=AUS, lw=1.4, label="lineage without DCX"),
           matplotlib.patches.Patch(color=INK, label="progenetic miniature"),
           matplotlib.patches.Patch(color=SINORT, label="proportioned dwarf")]
    a.legend(handles=hdl, loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=6.2, handlelength=1.4)
    ty = 0.955
    fig.text(0.70, ty + 0.035, "Losses (Dollo, stem interval)", fontsize=6.4, weight="bold")
    for k, (n, nom, lim) in enumerate(txt):
        fig.text(0.70, ty - 0.022 * k, f"{n}  ", fontsize=6.2, weight="bold")
        fig.text(0.725, ty - 0.022 * k, f"{nom}, {lim}", fontsize=6.2, style="italic" if nom != "Cobitis + Misgurnus" else "italic")
    L = R3 + 16; a.set_xlim(-L, L); a.set_ylim(-L, L)
    guardar(fig, "fig6_census_tree")

TODAS = ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "figS1", "figS2", "figS3", "figS4"]
for _f in (sys.argv[1:] or TODAS): globals()[_f]()
print("DONE")
