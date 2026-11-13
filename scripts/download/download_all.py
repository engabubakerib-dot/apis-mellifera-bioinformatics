#!/usr/bin/env python3
"""Downloads every file in data/download_list.tsv into data/raw/<accession>/ (stdlib only; replaces the PowerShell script,
so no execution-policy change is needed). Rules: never overwrites an existing file; files whose SHA-256 matches are skipped;
partial downloads go to <name>.part (an existing .part from an earlier attempt is kept as <name>.part.prev, never deleted);
a checksum mismatch is kept as <name>.new and reported. Writes only inside the project. Log: logs/download_python.log"""
import csv,hashlib,os,sys,time,urllib.request
R=os.path.abspath(sys.argv[1] if len(sys.argv)>1 else '.'); LOG=os.path.join(R,'logs','download_python.log'); os.makedirs(os.path.dirname(LOG),exist_ok=True)
def log(m):
    s=time.strftime('%Y-%m-%dT%H:%M:%S')+'  '+m; print(s,flush=True); open(LOG,'a',encoding='utf-8').write(s+'\n')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
ok=fail=0
for r in csv.DictReader(open(os.path.join(R,'data','download_list.tsv'),encoding='utf-8'),delimiter='\t'):
    d=os.path.join(R,'data','raw',r['accession']); os.makedirs(d,exist_ok=True); out=os.path.join(d,r['file_name'])
    if os.path.exists(out):
        h=sha(out)
        if r['sha256'] and h==r['sha256']: log(f"OK present+verified {r['file_name']}"); ok+=1; continue
        if not r['sha256']: log(f"present (no reference hash) {r['file_name']} sha256={h}"); ok+=1; continue
        log(f"present but checksum differs - NOT overwritten: {r['file_name']} {h}"); fail+=1; continue
    part=out+'.part'
    if os.path.exists(part):
        prev=part+'.prev'; i=1
        while os.path.exists(prev): prev=f'{part}.prev{i}'; i+=1
        os.rename(part,prev); log(f'earlier partial kept as {os.path.basename(prev)}')
    log('GET '+r['url'])
    try:
        with urllib.request.urlopen(urllib.request.Request(r['url'],headers={'User-Agent':'Apis-project/1.0'}),timeout=120) as resp, open(part,'wb') as f:
            while True:
                b=resp.read(1<<20)
                if not b: break
                f.write(b)
    except Exception as e: log(f"  FAILED {r['file_name']}: {e}"); fail+=1; continue
    h=sha(part); size=os.path.getsize(part)
    if r['sha256'] and h!=r['sha256']:
        os.rename(part,out+'.new'); log(f"  CHECKSUM MISMATCH {r['file_name']} got {h}; kept as .new"); fail+=1; continue
    os.rename(part,out); log(f"  verified {r['file_name']} bytes={size} sha256={h}"); ok+=1
log(f'DONE ok={ok} failed={fail}')
