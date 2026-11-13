#!/usr/bin/env python3
"""Sequence-level verification of NCBI GeneID <-> OGS v3.2 (GB) mapping.
Inputs : data/raw/BeeBase/amel_OGSv3.2_pep.fa.gz, data/raw/BeeBase/amel_OGSv3.2.gff3.gz,
         data/raw/NCBI/target_refseq_proteins.fa (RefSeq proteins of the target GeneIDs; header contains accession),
         data/raw/NCBI/gene_to_refseq_protein.tsv (gene_id, protein_accession)
Method : for every target GeneID, each RefSeq isoform is compared with every OGS3.2 peptide by shared 5-mer count
         (index search); the best hit is then aligned (global Needleman-Wunsch, BLOSUM-free identity scoring) to give
         % identity over the alignment. Mapping is CONFIRMED when the alias-declared GB ID is the best hit with >=95%
         identity; REASSIGNED when another GB is the best hit; UNRESOLVED otherwise.
Output : results/gene_mapping/ogs32_sequence_verification_LOCALFILES.tsv (does not overwrite the browser-derived table)"""
import gzip, os, re, sys, collections, csv
import numpy as np, pandas as pd
R=sys.argv[1] if len(sys.argv)>1 else '.'
def fasta(p):
    op=gzip.open if p.endswith('.gz') else open; s={}; k=None
    with op(p,'rt') as f:
        for l in f:
            l=l.strip()
            if l.startswith('>'): k=l[1:].split()[0]; s[k]=[]
            elif k: s[k].append(l)
    return {k:''.join(v).rstrip('*') for k,v in s.items()}
ogs={re.sub(r'.*\|','',k):v for k,v in fasta(f'{R}/data/raw/BeeBase/amel_OGSv3.2_pep.fa.gz').items()}
ncbi=fasta(f'{R}/data/raw/NCBI/target_refseq_proteins.fa')
g2p=pd.read_csv(f'{R}/data/raw/NCBI/gene_to_refseq_protein.tsv',sep='\t',dtype=str)
tg=pd.read_csv(f'{R}/data/metadata/gene_targets_ncbi.tsv',sep='\t',dtype=str)
K=5; idx=collections.defaultdict(set)
for name,s in ogs.items():
    for i in range(len(s)-K+1): idx[s[i:i+K]].add(name)
def nw_identity(a,b,gap=-2):
    n,m=len(a),len(b); S=np.zeros((n+1,m+1),np.int32); S[:,0]=np.arange(n+1)*gap; S[0,:]=np.arange(m+1)*gap
    A=np.frombuffer(a.encode(),np.uint8); B=np.frombuffer(b.encode(),np.uint8)
    for i in range(1,n+1):
        match=np.where(B==A[i-1],2,-1); diag=S[i-1,:-1]+match; up=S[i-1,1:]+gap
        row=np.maximum(diag,up); cur=np.empty(m+1,np.int32); cur[0]=S[i,0]
        for j in range(1,m+1): cur[j]=max(row[j-1],cur[j-1]+gap)
        S[i]=cur
    i,j,idn,L=n,m,0,0
    while i>0 and j>0:
        if S[i,j]==S[i-1,j-1]+(2 if a[i-1]==b[j-1] else -1): idn+=a[i-1]==b[j-1]; i-=1; j-=1
        elif S[i,j]==S[i-1,j]+gap: i-=1
        else: j-=1
        L+=1
    L+=i+j; return idn/L if L else 0
rows=[]
for _,t in tg.iterrows():
    accs=g2p[g2p.gene_id==t.ncbi_gene_id].protein_accession.tolist()
    best=None
    for acc in accs:
        q=next((v for k,v in ncbi.items() if k.startswith(acc) or acc in k),None)
        if not q: continue
        c=collections.Counter()
        for i in range(len(q)-K+1):
            for nm in idx.get(q[i:i+K],()): c[nm]+=1
        for nm,sc in c.most_common(2):
            if best is None or sc>best[2]: best=(acc,nm,sc,q)
    declared=str(t.ogs32_id).rstrip('?')
    if best is None: rows.append(dict(symbol=t.symbol,gene_id=t.ncbi_gene_id,declared_ogs32=t.ogs32_id,best_ogs32_peptide='',status='UNRESOLVED (no RefSeq protein / no k-mer hit)')); continue
    acc,nm,sc,q=best; gb=nm.split('-')[0]; s=ogs[nm]
    ident=nw_identity(q,s) if len(q)*len(s)<4e6 else float('nan')
    decl_ident=float('nan')
    if declared.startswith('GB') and declared!=gb:
        dp=[k for k in ogs if k.startswith(declared+'-')]
        if dp: decl_ident=max(nw_identity(q,ogs[k]) for k in dp if len(q)*len(ogs[k])<4e6)
    if gb==declared and ident>=0.95: st='CONFIRMED'
    elif gb==declared: st=f'CONFIRMED-PARTIAL (best hit, identity {ident:.2f})'
    elif declared.startswith('GB'): st=f'CONFLICT: best hit {gb}; declared {declared} identity {decl_ident:.2f}'
    else: st=f'ASSIGNED BY SEQUENCE: {gb}'
    rows.append(dict(symbol=t.symbol,gene_id=t.ncbi_gene_id,declared_ogs32=t.ogs32_id,refseq_protein=acc,best_ogs32_peptide=nm,shared_5mers=sc,
                     refseq_len=len(q),ogs_len=len(s),global_identity=round(ident,3),status=st))
os.makedirs(f'{R}/results/gene_mapping',exist_ok=True)
pd.DataFrame(rows).to_csv(f'{R}/results/gene_mapping/ogs32_sequence_verification_LOCALFILES.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).status.str.split(':').str[0].str.split(' \\(').str[0].value_counts())
