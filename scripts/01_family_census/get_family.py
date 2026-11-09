import urllib.request,urllib.parse,json,time,sys,re
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u,**p):
    for i in range(4):
        try: return urllib.request.urlopen(E+u+"?"+urllib.parse.urlencode(p),timeout=60).read().decode()
        except Exception as e: time.sleep(2)
    raise
species=sys.argv[2:]; out=open(sys.argv[1],"w"); info=[]
for sp in species:
    term=f'(doublecortin[title] OR dclk*[sym] OR dcx[sym] OR dcdc*[sym]) AND "{sp}"[orgn] AND alive[prop] AND genetype protein coding[prop]'
    ids=re.findall(r"<Id>(\d+)</Id>",get("esearch.fcgi",db="gene",term=term,retmax=100))
    if not ids: print("none",sp); continue
    sm=json.loads(get("esummary.fcgi",db="gene",id=",".join(ids),retmode="json"))["result"]
    for g in ids:
        r=sm[g]; desc=r["description"]
        if "reductase" in desc or "RP1" in desc.upper(): continue
        pids=re.findall(r"<Id>(\d+)</Id>",get("elink.fcgi",dbfrom="gene",db="protein",linkname="gene_protein_refseq",id=g))
        if not pids: continue
        ps=json.loads(get("esummary.fcgi",db="protein",id=",".join(pids[:60]),retmode="json"))["result"]
        best=max((ps[p] for p in ps["uids"]),key=lambda x:x["slen"])
        fa=get("efetch.fcgi",db="protein",id=best["caption"],rettype="fasta")
        seq="".join(fa.split("\n")[1:])
        tag=f'{sp.split()[0][:3]}{sp.split()[1][:3]}_{r["name"]}_{g}'.replace(" ","")
        out.write(f">{tag}\n{seq}\n"); info.append((tag,best["caption"],best["slen"],desc))
        time.sleep(0.35)
for i in info: print(*i,sep="\t")
