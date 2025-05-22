from Bio import Entrez as E,SeqIO as S
import matplotlib.pyplot as p,csv
e,k=input("e:"),input("k:")
t=input("t:")
a,b=int(input("a:")),int(input("b:"))
E.email,E.api_key=e,k
s=E.read(E.esearch(db="nucleotide",term=f"txid{t}[Organism]",usehistory="y",retmax=100))
r=list(S.parse(E.efetch(db="nucleotide",rettype="gb",retmode="text",webenv=s["WebEnv"],query_key=s["QueryKey"],retmax=100),"genbank"))
f=[x for x in r if a<=len(x.seq)<=b]
csv.writer(open(f"t{t}.csv","w",newline="")).writerows([["A","L","D"]]+[[x.id,len(x.seq),x.description]for x in f])
f.sort(key=lambda x:-len(x.seq))
p.plot([x.id for x in f],[len(x.seq) for x in f],'o-');p.xticks(rotation=90);p.tight_layout();p.savefig(f"t{t}.png")