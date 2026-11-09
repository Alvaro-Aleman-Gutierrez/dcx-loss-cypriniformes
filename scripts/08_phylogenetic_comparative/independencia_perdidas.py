#!/usr/bin/env python3
"""Numero MINIMO de perdidas independientes de DCX bajo parsimonia de Dollo,
evaluado sobre topologias alternativas de la literatura.
Motivo: la tesis del articulo ("N perdidas independientes") es una afirmacion sobre la TOPOLOGIA,
no sobre los genomas. La posicion de Paedocypris esta discutida, asi que hay que decir bajo que
colocaciones el numero se mantiene y bajo cuales baja.
Dollo: el gen estaba en el ancestro y solo puede perderse. Perdidas minimas = numero de clados
maximales cuyos descendientes son TODOS ausentes (un clado con algun presente no puede ser una perdida unica).
"""
from io import StringIO
from Bio import Phylo
import itertools, os
OUT = "${DATA}/independencia"; os.makedirs(OUT, exist_ok=True)
AUSENTE = {
 "Danio_rerio","Danio_aesculapii","Danio_kyathit","Danio_albolineatus","Danio_tinwini",
 "Danio_choprae","Danio_jaintianensis","Danio_nigrofasciatus",      # 8/8, lecturas crudas
 "Danionella_dracula",                                              # lecturas crudas
 "Paedocypris_micromegethes","Paedocypris_carbunculus","Paedocypris_sp_Singkep",  # lecturas (micromegethes)
 "Misgurnus_anguillicaudatus","Paramisgurnus_dabryanus",            # lecturas crudas
}
PRESENTE = {
 "Danionella_cerebrum","Rasbora_lateristriata","Pimephales_promelas","Gobio_gobio",
 "Ctenopharyngodon_idella","Labeo_rohita","Barbatula_barbatula","Triplophysa_rosa",
 "Astyanax_mexicanus",
}
DANIO = "((((Danio_rerio,Danio_aesculapii),(Danio_kyathit,(Danio_tinwini,Danio_nigrofasciatus))),Danio_albolineatus),(Danio_choprae,Danio_jaintianensis))"
PAEDO = "(Paedocypris_micromegethes,(Paedocypris_carbunculus,Paedocypris_sp_Singkep))"
DANIONELLA = "(Danionella_cerebrum,Danionella_dracula)"
COBITOIDEI = "((Misgurnus_anguillicaudatus,Paramisgurnus_dabryanus),(Barbatula_barbatula,Triplophysa_rosa))"
CYPRINIDAE = "(Ctenopharyngodon_idella,Labeo_rohita)"
LEUCISCIDAE = "(Gobio_gobio,Pimephales_promelas)"
DANIONINAE = f"({DANIONELLA},{DANIO})"            # Danionella hermana del resto de Danioninae
DANIONIDAE = f"(Rasbora_lateristriata,{DANIONINAE})"
def cyprinoidei(variante):
    if variante == "A":   # Paedocypris hermano de Danionidae
        return f"(({PAEDO},{DANIONIDAE}),({CYPRINIDAE},{LEUCISCIDAE}))"
    if variante == "B":   # Paedocypris hermano de TODO Cyprinoidei (posicion basal)
        return f"({PAEDO},({DANIONIDAE},({CYPRINIDAE},{LEUCISCIDAE})))"
    if variante == "C":   # Paedocypris dentro de Danioninae, hermano de Danionella+Danio
        return f"((Rasbora_lateristriata,({PAEDO},{DANIONINAE})),({CYPRINIDAE},{LEUCISCIDAE}))"
    if variante == "D":   # Paedocypris hermano SOLO del clado Danio
        dn = f"({DANIONELLA},({PAEDO},{DANIO}))"
        return f"((Rasbora_lateristriata,{dn}),({CYPRINIDAE},{LEUCISCIDAE}))"
    if variante == "E":   # Paedocypris hermano de Danionella
        dn = f"(({PAEDO},{DANIONELLA}),{DANIO})"
        return f"((Rasbora_lateristriata,{dn}),({CYPRINIDAE},{LEUCISCIDAE}))"
    if variante == "F":   # Paedocypris hermano de Cobitoidei (hipotesis extrema, poco apoyada)
        return f"({DANIONIDAE},({CYPRINIDAE},{LEUCISCIDAE}))"
DESC = {
 "A": "Paedocypris hermano de Danionidae (Cyprinoidei)",
 "B": "Paedocypris en posicion basal dentro de Cyprinoidei",
 "C": "Paedocypris dentro de Danioninae, hermano de (Danionella+Danio)",
 "D": "Paedocypris hermano SOLO del clado Danio",
 "E": "Paedocypris hermano de Danionella",
 "F": "Paedocypris fuera de Cyprinoidei, hermano de Cobitoidei",
}
def arbol(v):
    if v == "F":
        return f"(Astyanax_mexicanus,(({PAEDO},{COBITOIDEI}),{cyprinoidei('F')}));"
    return f"(Astyanax_mexicanus,({COBITOIDEI},{cyprinoidei(v)}));"
def perdidas_minimas(nw):
    t = Phylo.read(StringIO(nw), "newick")
    tips = [x.name for x in t.get_terminals()]
    desconocidos = [x for x in tips if x not in AUSENTE and x not in PRESENTE]
    assert not desconocidos, f"sin estado: {desconocidos}"
    faltan = (AUSENTE | PRESENTE) - set(tips)
    assert not faltan, f"taxones con estado pero AUSENTES del arbol: {sorted(faltan)}"
    def todos_ausentes(cl):
        return all(x.name in AUSENTE for x in cl.get_terminals())
    perdidas = []
    def rec(cl, padre_ausente):
        if todos_ausentes(cl):
            if not padre_ausente:
                perdidas.append(sorted(x.name for x in cl.get_terminals()))
            return
        for h in cl.clades: rec(h, False)
    rec(t.root, False)
    return perdidas, tips
lineas = []
def say(s):
    print(s, flush=True); lineas.append(s)
say("# Perdidas minimas de DCX bajo parsimonia de Dollo, por topologia")
say(f"# taxones ausentes: {len(AUSENTE)}   presentes: {len(PRESENTE)}")
res = {}
for v in "ABCDEF":
    nw = arbol(v)
    perd, tips = perdidas_minimas(nw)
    res[v] = len(perd)
    say("")
    say(f"== Topologia {v}: {DESC[v]}   -> PERDIDAS MINIMAS = {len(perd)}")
    for p in perd:
        etiqueta = ("clado Danio" if all(x.startswith("Danio_") for x in p) else
                    "Paedocypris" if all(x.startswith("Paedocypris") for x in p) else
                    "Cobitidae" if any("isgurnus" in x for x in p) else
                    "Danionella dracula" if p == ["Danionella_dracula"] else "mixto")
        say(f"   - {etiqueta}: {', '.join(p)}")
    open(f"{OUT}/topologia_{v}.nwk", "w").write(nw + "\n")
say("")
say("## RESUMEN")
for v in "ABCDEF": say(f"   {v} ({DESC[v]}): {res[v]} perdidas")
mn, mx = min(res.values()), max(res.values())
say(f"   => rango: {mn}-{mx} perdidas independientes segun la colocacion de Paedocypris")
if mn == mx:
    say(f"   => CONCLUSION: el numero ({mn}) es ROBUSTO a todas las colocaciones evaluadas.")
else:
    bajos = [v for v in res if res[v] == mn]
    say(f"   => CONCLUSION: NO es robusto. Baja a {mn} bajo: " + "; ".join(f"{v} ({DESC[v]})" for v in bajos))
    say(f"      Para afirmar {mx} hay que excluir esas topologias con datos, o declarar 'al menos {mn}'.")
open(f"{OUT}/independencia_results.txt", "w").write("\n".join(lineas) + "\n")
