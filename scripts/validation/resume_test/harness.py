"""Portable (Windows + Linux) test harness for scripts/rawreads/fetch_and_quant.py.
Everything is created under a caller-supplied base folder; the project's own data are never touched.
kallisto is replaced by a Python stand-in (fake_kallisto.py) for the mock tests; FASTQs are local synthetic files."""
import os,sys,gzip,json,hashlib,importlib.util,random,subprocess,time
HERE=os.path.dirname(os.path.abspath(__file__))
SCRIPT=os.environ.get('FQ_SCRIPT') or os.path.normpath(os.path.join(HERE,'..','..','rawreads','fetch_and_quant.py'))
FAKE_ID='FAKE_KALLISTO'
FAKE_SRC=r'''import sys,json,gzip,os,time
a=sys.argv
if a[1]=='version': print('kallisto, version FAKE'); sys.exit()
od=a[a.index('-o')+1]; fq=a[-1]; mode=os.environ.get('FAKEK','ok')
with gzip.open(fq) as f: n=sum(1 for _ in f)//4
if mode=='hang':
    with open(os.path.join(od,'abundance.tsv'),'w') as f: f.write('partial')
    time.sleep(120)
if mode=='fail': sys.stderr.write('simulated failure'); sys.exit(1)
if mode=='partial':
    with open(os.path.join(od,'abundance.tsv'),'w') as f: f.write('target_id\tlen\n')
    sys.exit(137)
# mimic kallisto 0.51.1 on Windows: "call" holds raw backslashes -> the file is NOT valid JSON
call='C:\\Users\\fake\\kallisto.exe quant -o '+od.replace('/','\\')+' '+fq.replace('/','\\')
with open(os.path.join(od,'run_info.json'),'w') as f: f.write('{\n\t"n_targets": 1,\n\t"n_processed": %d,\n\t"n_pseudoaligned": %d,\n\t"p_pseudoaligned": 100.0,\n\t"kallisto_version": "FAKE",\n\t"call": "%s"\n}\n'%(n,n,call))
with open(os.path.join(od,'abundance.tsv'),'w') as f: f.write('target_id\tlength\teff_length\test_counts\ttpm\nGB1-RA\t100\t50\t'+str(n)+'\t1e6\n')
'''
class Ctx:
    def __init__(s,base):
        s.base=os.path.abspath(base); s.R=os.path.join(s.base,'R'); s.src=os.path.join(s.base,'src')
def setup(base):
    c=Ctx(base); random.seed(1)
    for d in ['data/metadata','data/raw/BeeBase','data/reference','logs','data/raw/reads_subsample']: os.makedirs(os.path.join(c.R,d),exist_ok=True)
    os.makedirs(c.src,exist_ok=True)
    for acc,n in [('RUN1',3000),('RUN2',3000),('RUN3',3000),('SHORT',1500)]:
        with gzip.open(os.path.join(c.src,f'{acc}.fq.gz'),'wt') as f:
            for i in range(n): f.write(f'@{acc}.{i}\n'+''.join(random.choice('ACGT') for _ in range(60))+'\n+\n'+'I'*60+'\n')
    hdr='run_accession\tsample_alias\tsample_accession\tfastq_ftp\n'
    with open(os.path.join(c.R,'data/metadata/ena_PRJNA713798_runs.tsv'),'w',encoding='utf-8') as f: f.write(hdr+''.join(f'{a}\tGSM{a}\tS\t{src_url(c,a)}\n' for a in ['RUN1','RUN2','RUN3']))
    with open(os.path.join(c.R,'data/metadata/ena_PRJNA740269_runs.tsv'),'w',encoding='utf-8') as f: f.write(hdr+f'SHORT\tGSMS\tS\t{src_url(c,"SHORT")}\n')
    with open(os.path.join(c.R,'data/metadata/ena_PRJNA1017469_runs.tsv'),'w',encoding='utf-8') as f: f.write(hdr)
    for p in ['data/metadata/ena_PRJNA1017469_sample_attributes.tsv','data/raw/BeeBase/amel_OGSv3.2_cds.fa.gz','data/reference/amel_OGSv3.2_cds.kidx']:
        with open(os.path.join(c.R,p),'w') as f: f.write('x\n')
    with open(os.path.join(c.base,'fake_kallisto.py'),'w') as f: f.write(FAKE_SRC)
    return c
def src_url(c,acc,src=None):
    # the pipeline prepends https:// to fastq_ftp; the test get() maps it back to a local file
    return 'LOCAL'+os.path.join(src or c.src,f'{acc}.fq.gz')
def add_run(c,table,acc,src_acc):
    with open(os.path.join(c.R,'data/metadata',table),'a',encoding='utf-8') as f: f.write(f'{acc}\tGSM{acc}\tS\t{src_url(c,src_acc)}\n')
class Breaking:
    """File-like stream that raises a connection reset after `cut` bytes (cut>=0) or silently ends early (cut<0)."""
    def __init__(s,data,cut): s.d=data; s.p=0; s.cut=cut
    def read(s,n=-1):
        if s.cut<0 and s.p>=-s.cut: return b''
        if s.cut>=0 and s.p>=s.cut: raise ConnectionResetError('simulated connection reset')
        n=65536 if n is None or n<0 else n; n=min(n,abs(s.cut)-s.p,8192); out=s.d[s.p:s.p+n]; s.p+=len(out); return out
    def readable(s): return True
    def close(s): pass
    def __enter__(s): return s
    def __exit__(s,*a): pass
class SubprocShim:
    """Routes calls to FAKE_ID to `python fake_kallisto.py`; everything else goes to the real subprocess module."""
    def __init__(s,c): s.c=c
    def __getattr__(s,k): return getattr(subprocess,k)
    def run(s,cmd,*a,**k):
        if cmd and cmd[0]==FAKE_ID: cmd=[sys.executable,os.path.join(s.c.base,'fake_kallisto.py')]+list(cmd[1:])
        return subprocess.run(cmd,*a,**k)
def load_module(root):
    os.environ['FQ_TEST_HARNESS']='1'
    spec=importlib.util.spec_from_file_location('fq_under_test',os.environ.get('FQ_SCRIPT') or SCRIPT); m=importlib.util.module_from_spec(spec)
    old=sys.argv; sys.argv=['x',root]
    try: spec.loader.exec_module(m)
    finally: sys.argv=old
    return m
def run(c,N=2000,breaks=None,fakek='ok'):
    os.environ['FAKEK']=fakek
    m=load_module(c.R); m.N=N; m.refseq_and_mapping=lambda:None; m.kallisto=lambda:FAKE_ID; m.subprocess=SubprocShim(c); breaks=breaks or {}
    def get(url,*a,**k):
        if '://LOCAL' not in url: raise RuntimeError('no network in mock test: '+url)
        path=url.split('://LOCAL',1)[1]; acc=os.path.basename(path).split('.')[0]
        with open(path,'rb') as f: data=f.read()
        if acc in breaks: return Breaking(data,breaks[acc])
        return open(path,'rb')
    m.get=get; m.main()
def manifest_rows(c):
    p=os.path.join(c.R,'results/rawreads/run_manifest.tsv')
    with open(p,encoding='utf-8') as f: rows=[l.rstrip('\n').split('\t') for l in f]
    return [dict(zip(rows[0],r)) for r in rows[1:]]
def state(c):
    rows=manifest_rows(c); last={}
    for d in rows: last[(d['study'],d['subset'],d['run'])]=d['state']
    return last,len(rows)
def sha(p):
    with open(p,'rb') as f: return hashlib.sha256(f.read()).hexdigest()
def snap(c):
    out={}
    for dp,_,fs in os.walk(os.path.join(c.R,'results/rawreads/kallisto')):
        for f in fs: p=os.path.join(dp,f); out[os.path.relpath(p,c.R).replace('\\','/')]=(os.path.getmtime(p),sha(p))
    d=os.path.join(c.R,'data/raw/reads_subsample')
    for f in os.listdir(d): p=os.path.join(d,f); out['fq/'+f]=(os.path.getmtime(p),sha(p))
    return out
def kill_tree(p):
    if os.name=='nt': subprocess.run(['taskkill','/F','/T','/PID',str(p.pid)],capture_output=True)
    else:
        import signal; os.killpg(p.pid,signal.SIGKILL)
    p.wait()
