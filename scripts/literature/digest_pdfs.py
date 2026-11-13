import re,glob,os,json,collections
D='data/processed/literature_text'; out=[]
kw={'Vg':r'vitellogenin|\bVg\b','MRJP':r'royal jelly protein|\bMRJP|\bmrjp','antiox':r'catalase|superoxide|\bSOD\b|glutathione|peroxiredoxin|thioredoxin|antioxidant','immune':r'defensin|abaecin|hymenoptaecin|apidaecin|antimicrobial pept|immun','limit':r'limitation|further stud|future (?:stud|research|work)|remain(?:s)? unclear|not (?:yet )?known|warrant'}
def sents(t,pat,n):
    S=re.split(r'(?<=[.!?])\s+',t); r=[s.strip()[:300] for s in S if re.search(pat,s,re.I) and 40<len(s)<600]; return r[:n]
for f in sorted(glob.glob(D+'/*.txt')):
    t=open(f,errors='ignore').read(); tf=re.sub(r'\s+',' ',t)
    doi=re.findall(r'10\.\d{4,9}/[^\s,;)\]]+',tf); doi=collections.Counter([d.rstrip('.') for d in doi]).most_common(1)
    m=re.search(r'(?i)\babstract\b(.{0,1600})',tf); ab=(m.group(1) if m else tf[:1600])
    rec=dict(file=os.path.basename(f),doi=doi[0][0] if doi else '',head=tf[:350],abstract=ab[:1300],counts={k:len(re.findall(p,tf,re.I)) for k,p in kw.items()},
             Vg=sents(tf,kw['Vg'],3),MRJP=sents(tf,kw['MRJP'],3),antiox=sents(tf,kw['antiox'],2),immune=sents(tf,kw['immune'],2),limit=sents(tf,kw['limit'],3))
    out.append(rec)
json.dump(out,open('data/processed/literature_digest.json','w'),indent=1)
print(len(out))
