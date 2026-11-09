#!/usr/bin/env python3
"""Validador de paleta (equivalente en Python a validate_palette.js de la skill dataviz).
Para cada par: dE OKLab x100 en vision normal y bajo protanopia/deuteranopia/tritanopia
(Machado 2009, severidad 1), y diferencia de luminosidad L* (lectura en escala de grises).
Umbrales: dE CVD >= 8 objetivo, 6-8 solo con codificacion secundaria; dL* >= 15 para escala de grises.
Uso: python3 validar_paleta.py "#332288,#CC6677,#DDCC77" [nombres]"""
import sys, itertools, math
M = {"protan": [[0.152286,1.052583,-0.204868],[0.114503,0.786281,0.099216],[-0.003882,-0.048116,1.051998]],
     "deutan": [[0.367322,0.860646,-0.227968],[0.280085,0.672501,0.047413],[-0.011820,0.042940,0.968881]],
     "tritan": [[1.255528,-0.076749,-0.178779],[-0.078411,0.930809,0.147602],[0.004733,0.691367,0.303900]]}
def lin(c): return c/12.92 if c<=0.04045 else ((c+0.055)/1.055)**2.4
def rgb(h): h=h.lstrip("#"); return [lin(int(h[i:i+2],16)/255) for i in (0,2,4)]
def oklab(r):
    l=0.4122214708*r[0]+0.5363325363*r[1]+0.0514459929*r[2]
    m=0.2119034982*r[0]+0.6806995451*r[1]+0.1073969566*r[2]
    s=0.0883024619*r[0]+0.2817188376*r[1]+0.6299787005*r[2]
    l,m,s=[max(x,0)**(1/3) for x in (l,m,s)]
    return (0.2104542553*l+0.7936177850*m-0.0040720468*s, 1.9779984951*l-2.4285922050*m+0.4505937099*s,
            0.0259040371*l+0.7827717662*m-0.8086757660*s)
def sim(r,k): return [min(max(sum(M[k][i][j]*r[j] for j in range(3)),0),1) for i in range(3)]
def de(a,b): return 100*math.dist(oklab(a),oklab(b))
def Lstar(r):
    Y=0.2126*r[0]+0.7152*r[1]+0.0722*r[2]
    return 116*(Y**(1/3) if Y>0.008856 else 7.787*Y+16/116)-16
cols=sys.argv[1].split(","); nom=(sys.argv[2].split(",") if len(sys.argv)>2 else cols)
print("L* (gris):  "+"  ".join(f"{n}={Lstar(rgb(c)):.0f}" for n,c in zip(nom,cols)))
print(f"{'par':32s} {'normal':>7s} {'protan':>7s} {'deutan':>7s} {'tritan':>7s} {'dL*':>5s}  veredicto")
for (n1,c1),(n2,c2) in itertools.combinations(zip(nom,cols),2):
    a,b=rgb(c1),rgb(c2); v=[de(a,b)]+[de(sim(a,k),sim(b,k)) for k in M]; dl=abs(Lstar(a)-Lstar(b))
    mn=min(v[1:]); ver=("OK" if mn>=8 else "SOLO con forma/etiqueta" if mn>=6 else "FALLA CVD")
    ver+= "" if dl>=15 else " | gris: indistinguible sin forma"
    print(f"{n1+' / '+n2:32s} "+" ".join(f"{x:7.1f}" for x in v)+f" {dl:5.0f}  {ver}")
