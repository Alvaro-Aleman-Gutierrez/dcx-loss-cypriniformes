#!/usr/bin/env python3
"""Sintenia del bloque de DCX anclada en la ANOTACION, no en consultas genericas.
Correccion de un error propio: en syn_ortho.py las consultas "HIP1R"/"KNTC1" eran genericas y el mejor
locus del genoma resulto ser un parálogo en otro cromosoma (hip1ra). Eso me llevo a concluir que el
bloque no existia. La anotacion de Gobio lo desmiente: en NC_137100.1, en 780 kb,
   kntc1 - hip1rb - LOC307125817 - DCX - pak3loc - LOC307124806 - zfpl1 - rwdd2b - spata22 - lhfpl7
y el MISMO orden en Pimephales (NC_142772.1, 586 kb).
Aqui las consultas son las 10 proteinas REALES de esos loci de Gobio, y el paso reciproco exige que el
mejor hit en el proteoma COMPLETO de Gobio sea exactamente la proteina de partida: asi hip1ra no puede
pasar por hip1rb.
"""
import os, subprocess, collections, csv
P="${PROJECT}"; O="${DATA}/syn2"
B=f"{O}/bloque"; B2="${DATA}/p1b/genomes"
T=os.environ.get("THREADS","6")
GEN=["Danio_rerio","Danionella_dracula","Paedocypris_micromegethes","Paedocypris_sp_Pulau_Singkep",
     "Misgurnus_anguillicaudatus","Paramisgurnus_dabryanus",
     "Gobio_gobio","Pimephales_promelas","Rasbora_lateristriata","Ctenopharyngodon_idella",
     "Triplophysa_rosa","Barbatula_barbatula","Danionella_cerebrum","Labeo_rohita"]
SIN={"Danio_rerio","Danionella_dracula","Paedocypris_micromegethes","Paedocypris_sp_Pulau_Singkep",
     "Misgurnus_anguillicaudatus","Paramisgurnus_dabryanus"}
filas=list(csv.DictReader(open(f"{B}/bloque_genes.tsv"),delimiter="\t"))
nombre={f"Gobio|{r['proteina']}": r["gen"] for r in filas}
orden=[r["gen"] for r in filas]
# proteina -> gen de Gobio (GFF completo): el RBH compara GENES, no proteinas, porque las
# isoformas del mismo gen (p.ej. DCX XP_079671135.1 / XP_079671136.1) se rechazaban como "otro gen"
import re
P2G={}
for l in open("${DATA}/syn2/gff/Gobio.gff"):
    if "\tCDS\t" in l:
        a=re.search(r"protein_id=([^;\s]+)",l); g=re.search(r";gene=([^;\s]+)",l)
        if a and g: P2G["Gobio|"+a.group(1)]=g.group(1)
assert P2G.get("Gobio|XP_079671136.1")=="LOC307125267"
ETQ={"LOC307125267":"DCX","LOC307124699":"pak3loc"}
def et(g): return ETQ.get(g,g)
def run(c,**k): return subprocess.run(c,check=True,capture_output=True,text=True,**k).stdout
def db(sp):
    for c in (f"{P}/genomes/db_{sp}", f"{B2}/db_{sp}"):
        if os.path.exists(c+".ndb") or os.path.exists(c+".nal"): return c
out=open(f"{B}/bloque_results.txt","w")
def say(*a): print(*a,flush=True); out.write(" ".join(map(str,a))+"\n"); out.flush()
say("# Bloque de DCX con anclajes de la anotacion de Gobio (orden de referencia):")
say("# " + " - ".join(et(g) for g in orden))
tabla={}
for sp in GEN:
    d=db(sp)
    if not d: say(f"\n== {sp}: sin base"); continue
    tf=f"{B}/{sp}.tbn.tsv"
    if not os.path.exists(tf) or not os.path.getsize(tf):
        open(tf,"w").write(run(["tblastn","-query",f"{B}/q_bloque.faa","-db",d,"-evalue","1e-10","-seg","no",
            "-outfmt","6 qseqid sseqid pident length sstart send bitscore sseq","-num_threads",T,
            "-max_target_seqs","15","-max_hsps","5"]))
    por=collections.defaultdict(list)
    for l in open(tf):
        f_=l.rstrip("\n").split("\t")
        if len(f_)<8: continue
        c=f_[1].split("|")[1] if f_[1].count("|")>=2 else f_[1]
        por[(f_[0],c)].append(f_)
    loci=[]
    for (q,c),hs in por.items():
        hs.sort(key=lambda x:min(int(x[4]),int(x[5])))
        grp=[]
        for h in hs:
            if grp and min(int(h[4]),int(h[5]))-max(int(grp[-1][4]),int(grp[-1][5]))<100000: grp.append(h)
            else:
                if grp: loci.append((q,c,grp))
                grp=[h]
        if grp: loci.append((q,c,grp))
    faa=f"{B}/{sp}.loci.faa"; meta=[]
    with open(faa,"w") as fo:
        for i,(q,c,grp) in enumerate(loci):
            mej=max(grp,key=lambda x:float(x[6]))
            fo.write(f">L{i}\n{mej[7].replace('-','').replace('*','X')}\n")
            meta.append(dict(i=i,q=q,c=c,bs=sum(float(x[6]) for x in grp),pid=float(mej[2]),
                             pos=min(min(int(h[4]),int(h[5])) for h in grp)))
    if not meta: say(f"\n== {sp}: sin HSPs"); continue
    mr={}
    for l in run(["blastp","-query",faa,"-db",f"{B}/gobiodb","-outfmt","6 qseqid sseqid bitscore",
                  "-max_target_seqs","1","-max_hsps","1","-num_threads",T]).splitlines():
        a,b_,_=l.split("\t")
        if a not in mr: mr[a]=b_
    mejor={}
    for m in meta:
        R=mr.get(f"L{m['i']}")
        m["ok"]= R is not None and P2G.get(R,R)==P2G.get(m["q"],m["q"])   # RBH: vuelve al MISMO GEN de Gobio
        m["R"]=R or "-"
        g=nombre[m["q"]]
        if m["ok"] and (g not in mejor or m["bs"]>mejor[g]["bs"]): mejor[g]=m
    say(f"\n== {sp} [{'SIN DCX' if sp in SIN else 'con DCX'}]  verificados {len(mejor)}/{len(orden)}")
    cc=collections.Counter(m["c"] for m in mejor.values())
    princ=cc.most_common(1)[0] if cc else ("-",0)
    for g in orden:
        m=mejor.get(g)
        if m: say(f"   {et(g):16s} {m['c']:22s} {m['pos']:>12,}  bit={m['bs']:6.0f} pid={m['pid']:.0f}%")
        else:
            cand=[x for x in meta if nombre[x['q']]==g]
            if cand:
                b=max(cand,key=lambda x:x["bs"])
                say(f"   {et(g):16s} AUSENTE tras RBH (mejor candidato bit={b['bs']:.0f} -> {b['R']} = {P2G.get(b['R'],'?')})")
            else: say(f"   {et(g):16s} AUSENTE (sin HSP)")
    say(f"   >> contig principal {princ[0]} con {princ[1]}/{len(orden)} genes del bloque")
    tabla[sp]={g:(mejor[g]["c"],mejor[g]["pos"]) if g in mejor else None for g in orden}
import json; json.dump({k:{g:(v if v is None else list(v)) for g,v in d.items()} for k,d in tabla.items()},
                       open(f"{B}/bloque_tabla.json","w"), indent=1)
out.write("\nFIN bloque\n"); out.close(); print("FIN")
