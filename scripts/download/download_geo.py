#!/usr/bin/env python3
"""Download the GEO files listed in data/DOWNLOAD_MANIFEST.tsv into data/raw/<accession>/ and verify SHA-256.
Resumable (skips files whose checksum already matches). Never overwrites a file whose checksum differs:
it saves the new copy as <name>.new and reports the mismatch.
Usage: python3 scripts/download/download_geo.py [--project-root .]"""
import csv, hashlib, os, sys, urllib.request, argparse
ap=argparse.ArgumentParser(); ap.add_argument('--project-root',default='.'); a=ap.parse_args()
root=a.project_root; man=os.path.join(root,'data','DOWNLOAD_MANIFEST.tsv')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
for r in csv.DictReader(open(man),delimiter='\t'):
    url=r['url']
    if not url.startswith('https://ftp.ncbi.nlm.nih.gov/'): continue
    d=os.path.join(root,'data','raw',r['accession']); os.makedirs(d,exist_ok=True)
    out=os.path.join(d,r['file_name'])
    if os.path.exists(out) and sha(out)==r['sha256']:
        print('OK (present)',out); continue
    tmp=out+('.new' if os.path.exists(out) else '.part')
    print('GET',url); urllib.request.urlretrieve(url,tmp)
    s=sha(tmp)
    if s!=r['sha256']: print('  WARNING checksum differs from manifest:',s,'kept as',tmp); continue
    if tmp.endswith('.part'): os.replace(tmp,out)
    print('  verified',out)
