#!/usr/bin/env python3
"""Hoja de verificacion de accesibilidad para cada figura: original | escala de grises (luminancia CIE)
| deuteranopia | protanopia | tritanopia (Machado et al. 2009, severidad 1, en RGB lineal).
Salida: figuras/verificacion/<figura>_verif.png
Uso: python3 -I verificar_figuras.py [figuras/fig1_arbol_dcx.png ...]   (sin argumentos: todas)"""
import sys, os, glob
import numpy as np
from PIL import Image, ImageDraw, ImageFont
P = "${PROJECT}/figuras"; O = f"{P}/verificacion"; os.makedirs(O, exist_ok=True)
M = {"deuteranopia": [[0.367322,0.860646,-0.227968],[0.280085,0.672501,0.047413],[-0.011820,0.042940,0.968881]],
     "protanopia":   [[0.152286,1.052583,-0.204868],[0.114503,0.786281,0.099216],[-0.003882,-0.048116,1.051998]],
     "tritanopia":   [[1.255528,-0.076749,-0.178779],[-0.078411,0.930809,0.147602],[0.004733,0.691367,0.303900]]}
def a_lin(x): return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
def a_srgb(x): return np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(np.clip(x, 0, 1), 1 / 2.4) - 0.055)
def vistas(img):
    a = a_lin(np.asarray(img.convert("RGB"), dtype=float) / 255)
    out = [("original", img.convert("RGB"))]
    Y = a @ np.array([0.2126, 0.7152, 0.0722])
    g = (a_srgb(Y) * 255).round().astype(np.uint8)
    out.append(("escala de grises", Image.fromarray(np.stack([g] * 3, -1))))
    for k, m in M.items():
        s = np.clip(a @ np.array(m).T, 0, 1)
        out.append((k, Image.fromarray((a_srgb(s) * 255).round().astype(np.uint8))))
    return out
fs = sys.argv[1:] or sorted(f for f in glob.glob(f"{P}/fig*.png"))
for f in fs:
    img = Image.open(f); img.thumbnail((1100, 1100))
    vs = vistas(img); w, h = img.size; pad = 40
    hoja = Image.new("RGB", (w * 3 + 40, (h + pad) * 2 + 20), "white"); d = ImageDraw.Draw(hoja)
    try: fnt = ImageFont.truetype("DejaVuSans.ttf", 26)
    except Exception: fnt = ImageFont.load_default()
    for i, (nom, v) in enumerate(vs):
        x, y = (i % 3) * (w + 20), (i // 3) * (h + pad + 10)
        d.text((x + 10, y + 5), nom, fill="#222222", font=fnt); hoja.paste(v, (x, y + pad))
    dst = f"{O}/{os.path.basename(f)[:-4]}_verif.png"; hoja.save(dst); print("verificacion:", dst)
