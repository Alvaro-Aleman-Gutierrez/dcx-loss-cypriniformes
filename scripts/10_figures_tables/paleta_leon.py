"""Paleta derivada de Leon et al. 2024, Mol Biol Evol 41(9):msae166 (CC BY 4.0).
HEX medidos sobre las figuras de PMC11371425 (mediana de pixeles de relleno solido). Uso:
    from paleta_leon import LEON, aplicar_estilo; aplicar_estilo()
"""
LEON = {
    "presente": "#664BBC",  # violeta oscuro (Fig. 6)      -> DCX presente   (circulo lleno)
    "ausente":  "#FA41CA",  # magenta Humboldt (Figs 1,2,5) -> DCX ausente / perdida (circulo vacio, x)
    "parcial":  "#FEDC00",  # amarillo Magellanic          -> parcial (siempre con borde #525252)
    "externo":  "#71D226",  # verde Galapagos              -> grupo externo
    "sintenia": "#66FEFF",  # cian claro (Fig. 6)          -> vecinos de sintenia (con borde #525252)
    "hueco":    "#E5C0DA",  # tinte magenta (Fig. 2d)      -> region perdida
    "ctrl1": "#EBEBEB", "ctrl2": "#C4BAB9", "ctrl3": "#000000",   # controles internos (neutros del articulo)
    "texto": "#000000", "texto2": "#525252", "rejilla": "#EBEBEB", "fondo": "#FFFFFF",
}
# Secuencial (solo si aparece una variable continua): beige -> salmon -> rojo -> rojo oscuro,
# los 4 tonos estan en el articulo (mapa SST Fig. 6a y barras GO Fig. 3a); L* monotona 91 -> 31.
LEON_SECUENCIAL = ["#FDE2CD", "#F0ACA9", "#FE0000", "#9A0000"]
def aplicar_estilo():
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "text.color": LEON["texto"], "axes.labelcolor": LEON["texto"], "axes.edgecolor": LEON["texto2"],
        "xtick.color": LEON["texto"], "ytick.color": LEON["texto"], "axes.facecolor": LEON["fondo"],
        "figure.facecolor": LEON["fondo"], "grid.color": LEON["rejilla"], "legend.frameon": False,
        "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600})
    cmap = LinearSegmentedColormap.from_list("leon_seq", LEON_SECUENCIAL)
    try: mpl.colormaps.register(cmap)
    except ValueError: pass
    return cmap
