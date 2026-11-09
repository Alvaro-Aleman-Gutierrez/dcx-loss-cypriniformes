import urllib.request,urllib.parse,json,re,sys,time
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u,**p): return urllib.request.urlopen(E+u+"?"+urllib.parse.urlencode(p),timeout=60).read().decode()
gid,org,win=sys.argv[1],sys.argv[2],int(sys.argv[3])
r=json.loads(get("esummary.fcgi",db="gene",id=gid,retmode="json"))["result"][gid]
gi=r["genomicinfo"][0]; chrn=r["chromosome"]; a,b=int(gi["chrstart"]),int(gi["chrstop"]); lo,hi=min(a,b)-win,max(a,b)+win
print("#",r["name"],chrn,gi["chraccver"],a,b)
ids=re.findall(r"<Id>(\d+)</Id>",get("esearch.fcgi",db="gene",retmax=300,term=f'"{org}"[orgn] AND {chrn}[chr] AND {max(0,lo)}:{hi}[chrpos] AND alive[prop]'))
sm=json.loads(get("esummary.fcgi",db="gene",id=",".join(ids),retmode="json"))["result"] if ids else {"uids":[]}
rows=[]
for u in sm["uids"]:
    x=sm[u]; g=[k for k in x.get("genomicinfo",[]) if k.get("chraccver")==gi["chraccver"]]
    if not g: continue
    rows.append((int(g[0]["chrstart"]),u,x["name"],x["description"]))
for s,u,n,d in sorted(rows): print(s,u,n,d,sep="\t")
