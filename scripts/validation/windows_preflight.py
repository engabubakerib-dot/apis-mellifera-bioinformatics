r"""Windows preflight for the raw-read pipeline. It does NOT start the pipeline and does NOT modify project data.

Usage (from the project folder, with the launcher's interpreter):
  python scripts\validation\windows_preflight.py <project_dir>                 offline tests only
  python scripts\validation\windows_preflight.py <project_dir> --network       offline + ENA network test
Options: --scratch <dir>  (default: <project>\preflight_scratch\<timestamp>)

Every write goes into the scratch folder. Project files are only READ:
  scripts (hashes), data\metadata\ena_*_runs.tsv, data\reference\amel_OGSv3.2_cds.kidx (copied),
  data\raw\BeeBase\amel_OGSv3.2_cds.fa.gz (read to make synthetic reads), tools\kallisto\**\kallisto.exe (executed on scratch files).
  O12 also reads data\metadata\run_exclusions.tsv, results\rawreads\run_manifest.tsv and the pilot outputs (sha256 only).
The script under test is the DEPLOYED scripts\rawreads\fetch_and_quant.py (the file RUN_ALL_WINDOWS.bat runs).
If every test passes AND the network test ran, <project>\preflight_scratch\PREFLIGHT_PASSED.json is written with that file's SHA-256;
fetch_and_quant.py refuses to start the real pipeline without a marker matching its own SHA-256.
Full output: <scratch>\preflight_log.txt ; per-test table: <scratch>\preflight_results.tsv"""
import os, sys, time, json, hashlib, shutil, platform, subprocess, traceback, importlib, gzip, threading, locale, re
args=sys.argv[1:]
PROJ=os.path.abspath(args[0] if args and not args[0].startswith('--') else '.')
NETWORK='--network' in args
TS=time.strftime('%Y%m%d_%H%M%S')
SCR=os.path.abspath(args[args.index('--scratch')+1]) if '--scratch' in args else os.path.join(PROJ,'preflight_scratch',TS)
os.makedirs(SCR,exist_ok=False)
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(HERE,'resume_test'))
LOG=os.path.join(SCR,'preflight_log.txt'); RES=os.path.join(SCR,'preflight_results.tsv'); RESULTS=[]
SCRIPTS={'DEPLOYED':os.path.join(PROJ,'scripts','rawreads','fetch_and_quant.py')}
WIN=os.name=='nt'
class Tee:
    def __init__(s,*f): s.f=f
    def write(s,x):
        for f in s.f: f.write(x); f.flush()
    def flush(s):
        for f in s.f: f.flush()
_logf=open(LOG,'a',encoding='utf-8'); sys.stdout=Tee(sys.__stdout__,_logf); sys.stderr=Tee(sys.__stderr__,_logf)
def out(m=''): print(m,flush=True)
def rec(test,status,detail=''):
    RESULTS.append((test,status,detail)); out(f'{status:4} {test}'+(f'   [{detail}]' if detail else ''))
def T(test,cond,detail=''): rec(test,'PASS' if cond else 'FAIL',detail); return bool(cond)
def NA(test,why): rec(test,'N/A',why)
def section(title,fn):
    out(''); out('='*100); out(title); out('='*100)
    try: fn()
    except Exception as e: rec(title+' :: unexpected exception','FAIL',repr(e)); out(traceback.format_exc())
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def load(which,root):
    os.environ['FQ_SCRIPT']=SCRIPTS[which]; import harness as H; return H.load_module(root)

out(f'WINDOWS PREFLIGHT  {TS}   network tests: {"ON" if NETWORK else "OFF (offline only)"}')
out(f'project : {PROJ}'); out(f'scratch : {SCR}  (all writes go here)'); out(f'log     : {LOG}')

# ===================================================================================== OFFLINE TESTS
def o1():
    exp=os.path.join(os.environ.get('USERPROFILE',os.path.expanduser('~')),'anaconda3','python.exe')
    out(f'sys.executable={sys.executable}  version={sys.version}  platform={platform.platform()}')
    T('O1.1 interpreter = launcher interpreter %USERPROFILE%\\anaconda3\\python.exe',os.path.normcase(sys.executable)==os.path.normcase(exp),sys.executable) if WIN else NA('O1.1 interpreter = launcher interpreter','not Windows')
    T('O1.2 Python version >= 3.8',sys.version_info>=(3,8),sys.version.split()[0])
    bad=[]
    for m in ['os','sys','gzip','json','hashlib','shutil','subprocess','urllib.request','zipfile','tarfile','io','time','platform','xml.etree.ElementTree','ssl','zlib']:
        try: importlib.import_module(m)
        except Exception as e: bad.append(f'{m}: {e}')
    T('O1.3 all modules used by fetch_and_quant.py import (stdlib only, no extra packages)',not bad,'; '.join(bad))
    import ssl; T('O1.4 SSL module available (needed for https to ENA)',True,ssl.OPENSSL_VERSION)
    out(f'filesystem encoding={sys.getfilesystemencoding()}  preferred locale encoding={locale.getpreferredencoding(False)}  stdout={getattr(sys.__stdout__,"encoding",None)}')
section('O1  Python executable, version, dependencies',o1)

def o2():
    exp={}
    with open(os.path.join(HERE,'preflight_expected_sha256.txt'),encoding='utf-8') as f:
        for l in f:
            if l.strip() and not l.startswith('#'): h,p=l.split(None,1); exp[p.strip()]=h
    for rel,h in exp.items():
        p=os.path.join(PROJ,*rel.split('/')); got=sha(p) if os.path.exists(p) else 'MISSING'
        T(f'O2.1 version on disk: {rel}',got==h,f'sha256 {got[:16]} expected {h[:16]}')
    for rel in ['data/metadata/ena_PRJNA713798_runs.tsv','data/reference/amel_OGSv3.2_cds.kidx','data/raw/BeeBase/amel_OGSv3.2_cds.fa.gz','data/metadata/GSE178742_sample_metadata_verified.tsv']:
        T(f'O2.2 required input present: {rel}',os.path.exists(os.path.join(PROJ,*rel.split('/'))))
    longest=max([os.path.join(PROJ,'results','rawreads','kallisto','PRJNA1017469','A','SRR26123456.inprogress.failed1','kallisto_stderr.txt'),
                 os.path.join(PROJ,'data','raw','reads_subsample','SRR14817123_R1_subsetA.fastq.gz.part.prev1')],key=len)
    lp=None
    if WIN:
        try:
            import winreg; k=winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,r'SYSTEM\CurrentControlSet\Control\FileSystem'); lp=winreg.QueryValueEx(k,'LongPathsEnabled')[0]
        except Exception: lp=None
    T('O2.3 longest path the pipeline will create is < 260 characters (or long paths enabled)',len(longest)<250 or lp==1,f'{len(longest)} chars; LongPathsEnabled={lp}')
    for f in sorted(os.listdir(os.path.join(PROJ,'data','metadata'))):
        if f.endswith('.tsv'):
            p=os.path.join(PROJ,'data','metadata',f)
            try:
                with open(p,encoding='utf-8',errors='strict') as fh: fh.read(); ok=True; d=''
            except Exception as e: ok=False; d=repr(e)
            T(f'O2.4 strict UTF-8 read of data/metadata/{f}',ok,d)
section('O2  Project paths, versions on disk, file encodings (read-only)',o2)

KEXE=[]
def o3():
    ks=[os.path.join(dp,f) for dp,_,fs in os.walk(os.path.join(PROJ,'tools','kallisto')) for f in fs if f in ('kallisto.exe','kallisto')]
    T('O3.1 kallisto executable discovered by the pipeline rule (tools\\kallisto\\**\\kallisto.exe)',bool(ks),ks[0] if ks else 'none -> the pipeline would try to download it')
    if not ks: return
    KEXE.append(ks[0])
    if WIN: T('O3.2 discovered file is kallisto.exe (Windows binary)',ks[0].lower().endswith('kallisto.exe'),os.path.basename(ks[0]))
    v=subprocess.run([ks[0],'version'],capture_output=True,text=True,errors='replace')
    T('O3.3 kallisto runs and reports version 0.51.1',v.returncode==0 and '0.51.1' in v.stdout,(v.stdout+v.stderr).strip()[:120])
    flags={0x400000:'RECALL_ON_DATA_ACCESS(cloud-only)',0x40000:'RECALL_ON_OPEN',0x1000:'OFFLINE',0x80000:'PINNED',0x100000:'UNPINNED'}
    for p in [ks[0],os.path.join(PROJ,'data','reference','amel_OGSv3.2_cds.kidx')]:
        a=getattr(os.stat(p),'st_file_attributes',0); names=[n for b,n in flags.items() if a & b]
        T(f'O3.4 OneDrive: {os.path.relpath(p,PROJ)} is stored locally (not cloud-only)',not (a & (0x400000|0x40000|0x1000)),', '.join(names) or 'no placeholder flags')
    free=shutil.disk_usage(PROJ).free/1e9
    T('O3.5 free disk space > 10 GB (209 subsets x ~28 MB FASTQ + outputs ~ 7 GB)',free>10,f'{free:.1f} GB free')
section('O3  kallisto discovery and version, OneDrive file state, disk space',o3)

def o4():
    d=os.path.join(SCR,'o4_rename'); os.makedirs(d)
    for which in SCRIPTS:
        m=load(which,os.path.join(SCR,f'o4_root_{which}')); w=os.path.join(d,which); os.makedirs(w)
        # closed file -> folder rename works immediately
        f1=os.path.join(w,'closed'); os.makedirs(f1)
        with open(os.path.join(f1,'a.txt'),'w',encoding='utf-8') as f: f.write('x')
        try: os.rename(f1,f1+'.renamed'); ok=True; e=''
        except Exception as ex: ok=False; e=repr(ex)
        T(f'O4.1 [{which}] folder with a CLOSED file renames immediately',ok,e)
        # open handle inside folder -> Windows blocks rename (shows why explicit closing matters)
        f2=os.path.join(w,'open'); os.makedirs(f2); h=open(os.path.join(f2,'a.txt'),'w')
        try: os.rename(f2,f2+'.renamed'); blocked=False
        except PermissionError: blocked=True
        h.close()
        if WIN: T(f'O4.2 [{which}] Windows blocks renaming a folder while a file in it is open (control)',blocked,'PermissionError raised' if blocked else 'rename succeeded')
        else: NA(f'O4.2 [{which}] open-handle rename block','Linux does not block renames of open files')
        # keep_aside never overwrites, numbers duplicates
        p=os.path.join(w,'x.fastq.gz');
        for i in range(3):
            with open(p,'w') as f: f.write(str(i))
            m.keep_aside(p,'.prev')
        names=sorted(os.listdir(w)); exp=['x.fastq.gz.prev','x.fastq.gz.prev1','x.fastq.gz.prev2']
        vals=[open(os.path.join(w,n)).read() for n in exp if os.path.exists(os.path.join(w,n))]
        T(f'O4.3 [{which}] keep_aside renames to .prev/.prev1/.prev2 without overwriting',all(n in names for n in exp) and vals==['0','1','2'],str([n for n in names if n.startswith('x.')]))
section('O4  File closing and rename behaviour (scratch only)',o4)

def o5():
    for which in SCRIPTS:
        m=load(which,os.path.join(SCR,f'o5_root_{which}'))
        if not hasattr(m,'safe_rename'):
            NA(f'O5 [{which}] lock-retry tests','this version has no safe_rename (renames are not retried on a lock)')
            if WIN:
                d=os.path.join(SCR,f'o5_{which}'); os.makedirs(d); src=os.path.join(d,'held'); open(src,'w').close(); h=open(src,'r')
                try: m.keep_aside(src,'.prev'); r='renamed'
                except PermissionError as e: r='PermissionError (no retry)'
                h.close(); T(f'O5.0 [{which}] behaviour when a file is briefly locked: fails immediately',r=='PermissionError (no retry)',r)
            continue
        d=os.path.join(SCR,f'o5_{which}'); os.makedirs(d)
        # transient lock released after 3 s -> retried and succeeds
        src=os.path.join(d,'transient'); open(src,'w').write('data'); h=open(src,'r')
        threading.Timer(3.0,h.close).start(); t0=time.time()
        try: m.safe_rename(src,src+'.moved'); ok=os.path.exists(src+'.moved') and not os.path.exists(src); e=''
        except Exception as ex: ok=False; e=repr(ex)
        dt=time.time()-t0
        if WIN: T(f'O5.1 [{which}] transient lock (3 s) -> retried, then rename succeeds',ok and dt>=2,f'{dt:.1f} s {e}')
        else: T(f'O5.1 [{which}] rename succeeds (Linux: no lock)',ok,e)
        # persistent lock -> gives up with PermissionError; source intact, nothing created
        src=os.path.join(d,'held'); open(src,'w').write('keep'); h=open(src,'r')
        try: m.safe_rename(src,src+'.moved',tries=3); r='renamed'
        except PermissionError: r='PermissionError'
        h.close()
        if WIN: T(f'O5.2 [{which}] persistent lock -> PermissionError after retries; source intact; no destination',r=='PermissionError' and open(src).read()=='keep' and not os.path.exists(src+'.moved'),r)
        else: NA(f'O5.2 [{which}] persistent lock','Linux does not lock open files')
        # never overwrite an existing destination
        a=os.path.join(d,'a'); b=os.path.join(d,'b'); open(a,'w').write('A'); open(b,'w').write('B')
        try: m.safe_rename(a,b); r='overwrote'
        except FileExistsError: r='FileExistsError'
        T(f'O5.3 [{which}] destination exists -> refused (FileExistsError); both files unchanged',r=='FileExistsError' and open(a).read()=='A' and open(b).read()=='B',r)
section('O5  Windows file-lock retry behaviour (scratch only)',o5)

def o6():
    for which in SCRIPTS:
        root=os.path.join(SCR,f'o6_root_{which}'); m=load(which,root)
        m.manifest_add('GSE168740','A','SRRTEST1','DOWNLOAD_STARTED',reads_target=500000,note='url\twith tab\nand newline; non-ASCII Kadle\u010dkov\u00e1')
        m.manifest_add('GSE168740','A','SRRTEST1','FASTQ_VALIDATED',reads_target=500000,reads_written=500000,fastq_eof=False,fastq_sha256='ab'*32,fastq_bytes=123,reads_recounted=500000)
        p=m.manifest_path()
        with open(p,encoding='utf-8',errors='strict') as f: L=f.read().split('\n')
        rows=[l.split('\t') for l in L if l]
        T(f'O6.1 [{which}] manifest written with header + 2 rows, all with {len(m.MANIFEST_COLS)} columns (tabs/newlines in notes neutralised)',len(rows)==3 and all(len(r)==len(m.MANIFEST_COLS) for r in rows),f'{len(rows)} lines')
        T(f'O6.2 [{which}] manifest is UTF-8 and non-ASCII text round-trips',any('Kadle\u010dkov\u00e1' in l for l in L))
        if hasattr(m,'last_validated_fastq'):
            v=m.last_validated_fastq('GSE168740','A','SRRTEST1')
            T(f'O6.3 [{which}] manifest read-back: last FASTQ_VALIDATED record found',v and v['fastq_sha256']=='ab'*32 and v['reads_written']=='500000')
            m.manifest_add('GSE168740','A','SRRTEST1','DOWNLOAD_STARTED')
            T(f'O6.4 [{which}] a later DOWNLOAD_STARTED invalidates the earlier FASTQ_VALIDATED record',m.last_validated_fastq('GSE168740','A','SRRTEST1') is None)
section('O6  Manifest writing and reading (scratch only)',o6)

def o7():
    d=os.path.join(SCR,'o7'); os.makedirs(d)
    for which in SCRIPTS:
        m=load(which,os.path.join(SCR,f'o7_root_{which}'))
        big=os.path.join(d,f'big_{which}.bin')
        with open(big,'wb') as f: f.write(os.urandom(5*1024*1024+17))
        h=hashlib.sha256(open(big,'rb').read()).hexdigest()
        T(f'O7.1 [{which}] sha256_file equals hashlib on a 5 MB file',m.sha256_file(big)==h)
        good=os.path.join(d,f'good_{which}.fastq.gz')
        with gzip.open(good,'wb') as f: f.write(b''.join(b'@r%d\nACGT\n+\nIIII\n'%i for i in range(1000)))
        T(f'O7.2 [{which}] strict read count of a valid FASTQ = 1000',m.count_reads_strict(good)==1000)
        bad=os.path.join(d,f'bad_{which}.fastq.gz')
        with gzip.open(bad,'wb') as f: f.write(b'@r1\nACGT\n+\nIIII\n@r2\nACGT\n')
        try: m.count_reads_strict(bad); r='accepted'
        except ValueError: r='rejected'
        T(f'O7.3 [{which}] FASTQ with a truncated record is rejected',r=='rejected',r)
        trunc=os.path.join(d,f'trunc_{which}.fastq.gz'); data=open(good,'rb').read()
        with open(trunc,'wb') as f: f.write(data[:len(data)//2])
        try: m.count_reads_strict(trunc); r='accepted'
        except Exception as e: r='rejected: '+type(e).__name__
        T(f'O7.4 [{which}] gzip cut in half is rejected',r.startswith('rejected'),r)
section('O7  Checksum and read-count validation (scratch only)',o7)

def o8():
    import tests as TS
    for which in SCRIPTS:
        os.environ['FQ_SCRIPT']=SCRIPTS[which]
        for name,ok in TS.run_all(os.path.join(SCR,f'o8_mock_{which}'),report=lambda m:None): T(f'O8 [{which}] {name}',ok)
        for name,ok in TS.run_kill(os.path.join(SCR,f'o8_kill_{which}'),report=lambda m:None): T(f'O8 [{which}] {name} ({"taskkill /F /T" if WIN else "SIGKILL"})',ok)
section('O8  Recovery from incomplete files/folders; validated outputs never overwritten (mock kallisto, local FASTQs)',o8)

ORIG_REF=('reports/archive/pre_exclusion_20261112/scripts/rawreads/fetch_and_quant.py','c3e26f2834eb9aa35bcf1255b639dff0f24999b53b8897c29ac93c2a11a89b5c')
def o12():
    import exclusion_tests as ET
    orig=os.path.join(PROJ,*ORIG_REF[0].split('/'))
    got=sha(orig) if os.path.exists(orig) else 'MISSING'
    if not T('O12.0 archived pre-change script present and unaltered (reference for byte-identity)',got==ORIG_REF[1],f'{ORIG_REF[0]} sha256 {got[:16]}'): return
    os.environ['FQ_SCRIPT']=SCRIPTS['DEPLOYED']
    for name,ok,detail in ET.run_all(os.path.join(SCR,'o12_exclusions'),PROJ,orig,report=lambda m:None): T(f'O12 {name}',ok,detail)
section('O12 Run exclusions (approved 2026-11-12): real ENA tables, other studies byte-identical, audit rows, pilot reuse (project read-only)',o12)

ACC_ENV=['FQ_WORKERS','FQ_MAX_RESUMES','FQ_READ_TIMEOUT','FQ_JOB_MAX_SECONDS','FQ_RETRY_PASSES','FQ_SELECTION','FQ_SLOW_403','FQ_BACKOFF']
def o13():
    import accel_tests as AT
    os.environ['FQ_SCRIPT']=SCRIPTS['DEPLOYED']; AT.SCRIPT=SCRIPTS['DEPLOYED']
    try:
        for name,ok,detail in AT.run_all(os.path.join(SCR,'o13_accel'),PROJ,report=lambda m:None): T(f'O13 {name}',ok,detail)
    finally:
        for k in ACC_ENV: os.environ.pop(k,None)
section('O13 Acceleration options: byte-range resume, 1-3 parallel downloads, fail-fast, retry pass, run selection (mock ENA; project read-only)',o13)

def o9():
    if not KEXE: NA('O9 real kallisto.exe offline pipeline test','kallisto not found (see O3)'); return
    import harness as H
    cds=os.path.join(PROJ,'data','raw','BeeBase','amel_OGSv3.2_cds.fa.gz')
    for which in SCRIPTS:
        base=os.path.join(SCR,f'o9_real_kallisto_{which}'); R=os.path.join(base,'R'); src=os.path.join(base,'src')
        for x in ['data/metadata','data/reference','data/raw/BeeBase','logs']: os.makedirs(os.path.join(R,*x.split('/')),exist_ok=True)
        os.makedirs(src)
        seqs=[]; cur=[]
        with gzip.open(cds,'rt') as f:
            for l in f:
                if l.startswith('>'):
                    if cur and len(''.join(cur))>=300: seqs.append(''.join(cur))
                    cur=[]
                    if len(seqs)>=300: break
                else: cur.append(l.strip())
        import random; random.seed(7)
        with gzip.open(os.path.join(src,'SRRSYN1.fq.gz'),'wt') as f:
            for i in range(3000):
                s=random.choice(seqs); st=random.randint(0,len(s)-100); f.write(f'@SRRSYN1.{i+1}\n{s[st:st+100]}\n+\n'+'I'*100+'\n')
        hdr='run_accession\tsample_alias\tsample_accession\tfastq_ftp\n'
        with open(os.path.join(R,'data','metadata','ena_PRJNA713798_runs.tsv'),'w',encoding='utf-8') as f: f.write(hdr+'SRRSYN1\tGSMSYN\tS\tLOCAL'+os.path.join(src,'SRRSYN1.fq.gz')+'\n')
        for p in ['ena_PRJNA740269_runs.tsv','ena_PRJNA1017469_runs.tsv']:
            with open(os.path.join(R,'data','metadata',p),'w',encoding='utf-8') as f: f.write(hdr)
        with open(os.path.join(R,'data','metadata','ena_PRJNA1017469_sample_attributes.tsv'),'w') as f: f.write('x\n')
        idx=os.path.join(PROJ,'data','reference','amel_OGSv3.2_cds.kidx'); h0=sha(idx)
        shutil.copy2(idx,os.path.join(R,'data','reference','amel_OGSv3.2_cds.kidx'))
        def go():
            m=load(which,R); m.N=1000; m.refseq_and_mapping=lambda:None; m.kallisto=lambda:KEXE[0]
            def get(url,*a,**k):
                if '://LOCAL' not in url: raise RuntimeError('offline test: no network '+url)
                return open(url.split('://LOCAL',1)[1],'rb')
            m.get=get; m.main()
        c=H.Ctx(base); go(); st,_=H.state(c)
        T(f'O9.1 [{which}] real kallisto.exe: subsets A and B COMPLETE in the OneDrive scratch folder',st.get(('GSE168740','A','SRRSYN1'))=='COMPLETE' and st.get(('GSE168740','B','SRRSYN1'))=='COMPLETE',str(st))
        cj=os.path.join(R,'results','rawreads','kallisto','GSE168740','A','SRRSYN1','_COMPLETE.json')
        if os.path.exists(cj):
            with open(cj,encoding='utf-8') as f: rj=json.load(f)
            T(f'O9.2 [{which}] _COMPLETE.json: kallisto rc 0, n_processed = reads written = recount = 1000',rj['kallisto_rc']==0 and rj['n_processed']==1000==rj['reads_written']==rj['reads_recounted'])
            T(f'O9.3 [{which}] recorded SHA-256 = SHA-256 of FASTQ re-read from disk',sha(os.path.join(R,rj['fastq']))==rj['fastq_sha256'])
            T(f'O9.4 [{which}] reads sampled from the CDS pseudoalign at a high rate (index/binary sanity)',(rj.get('p_pseudoaligned') or 0)>80,f"{rj.get('p_pseudoaligned')}%")
        before=H.snap(c); _,a=H.state(c); time.sleep(1.5); go(); after=H.snap(c); _,b=H.state(c)
        T(f'O9.5 [{which}] restart: validated runs skipped - no new manifest rows, files and mtimes identical',before==after and a==b,f'manifest rows {a}->{b}')
        T(f'O9.6 [{which}] project index unchanged (read-only use)',sha(idx)==h0)
section('O9  Real kallisto.exe through the pipeline, offline (synthetic reads from the local CDS)',o9)

def o10():
    with open(os.path.join(PROJ,'scripts','R','de_from_kallisto.R'),encoding='utf-8') as f: rs=f.read()
    T('O10.1 de_from_kallisto.R selects only accession-named folders that contain _COMPLETE.json',"grepl('^[SED]RR[0-9]+$',basename(d))&file.exists(file.path(d,'_COMPLETE.json'))" in rs)
    d=os.path.join(SCR,'o10','results','rawreads','kallisto','GSE168740','A')
    for x,c in [('SRR0000001',True),('SRR0000002.failed',True),('SRR0000003.inprogress',False),('SRR0000004',False)]:
        os.makedirs(os.path.join(d,x))
        with open(os.path.join(d,x,'abundance.tsv'),'w') as f: f.write('target_id\n')
        if c:
            with open(os.path.join(d,x,'_COMPLETE.json'),'w') as f: f.write('{}')
    sel=sorted(x for x in os.listdir(d) if re.match(r'^[SED]RR[0-9]+$',x) and os.path.exists(os.path.join(d,x,'_COMPLETE.json')))
    T('O10.2 selection rule (emulated in Python) keeps only SRR0000001',sel==['SRR0000001'],str(sel))
    rs_exe=[p for p in [os.path.join(PROJ,'envs','r','Lib','R','bin','Rscript.exe')]+sorted(__import__('glob').glob(r'C:\Program Files\R\R-*\bin\Rscript.exe'))[-1:] if os.path.exists(p)]
    if rs_exe:
        code="d=list.dirs('%s',recursive=FALSE);d=d[grepl('^[SED]RR[0-9]+$',basename(d))&file.exists(file.path(d,'_COMPLETE.json'))];cat(basename(d))"%d.replace('\\','/')
        r=subprocess.run([rs_exe[0],'--vanilla','-e',code],capture_output=True,text=True,errors='replace')
        T('O10.3 same rule executed by real Rscript keeps only SRR0000001',r.stdout.strip()=='SRR0000001',f'{rs_exe[0]}: {r.stdout.strip()} {r.stderr.strip()[-150:]}')
    else: NA('O10.3 rule executed by real Rscript','no Rscript.exe found')
section('O10 R input filter (only folders with _COMPLETE.json)',o10)


def o11():
    """tx2gene (as built by scripts/R/tx2gene.R) must cover every kallisto target_id; ignoreTxVersion must be FALSE."""
    gff=os.path.join(PROJ,'data','raw','BeeBase','amel_OGSv3.2.gff3.gz'); cds=os.path.join(PROJ,'data','raw','BeeBase','amel_OGSv3.2_cds.fa.gz')
    tx={}
    with gzip.open(gff,'rt',encoding='utf-8') as f:
        for l in f:
            if l.startswith('#'): continue
            c=l.rstrip('\n').split('\t')
            if len(c)>8 and c[2]=='mRNA':
                a=dict(x.split('=',1) for x in c[8].split(';') if '=' in x); tx['gnl|Amel_4.5|'+a['ID']]=a['Parent']
    with gzip.open(cds,'rt',encoding='utf-8') as f: hdr=[l[1:].split()[0] for l in f if l.startswith('>')]
    T('O11.1 tx2gene (GFF3 mRNA IDs + prefix gnl|Amel_4.5|) covers every CDS header = every kallisto target_id',set(hdr)<=set(tx) and len(hdr)==len(set(hdr)),f'{len(hdr)} CDS targets, {len(tx)} tx2gene rows, {len(set(tx.values()))} genes, unmapped {len(set(hdr)-set(tx))}')
    ab=sorted(__import__('glob').glob(os.path.join(SCR,'o9_real_kallisto_*','R','results','rawreads','kallisto','*','*','*','abundance.tsv')))
    if ab:
        with open(ab[0],encoding='utf-8') as f: ids=[l.split('\t')[0] for l in f][1:]
        T('O11.2 target_ids in a real kallisto.exe abundance.tsv all map to tx2gene',set(ids)<=set(tx),f'{len(ids)} ids; unmapped {len(set(ids)-set(tx))}')
    collapsed=len({re.sub(r'\..*','',h) for h in hdr})
    T('O11.3 control: tximport ignoreTxVersion=TRUE would collapse the IDs (so it must be FALSE)',collapsed<len(hdr),f'{len(hdr)} ids -> {collapsed} after stripping at the first "."')
    with open(os.path.join(PROJ,'scripts','R','de_from_kallisto.R'),encoding='utf-8') as f: rs=f.read()
    T('O11.4 de_from_kallisto.R uses ignoreTxVersion=FALSE only and the shared tx2gene.R with check_tx2gene()',('ignoreTxVersion=TRUE' not in rs) and rs.count('ignoreTxVersion=FALSE')>=2 and 'check_tx2gene(' in rs and "source('scripts/R/tx2gene.R')" in rs)
section('O11 Transcript-to-gene mapping compatibility for tximport (read-only)',o11)

# ===================================================================================== NETWORK TEST (optional)
def n1():
    import harness as H
    base=os.path.join(SCR,'n1_ena'); R=os.path.join(base,'R')
    for x in ['data/metadata','data/reference','data/raw/BeeBase','logs']: os.makedirs(os.path.join(R,*x.split('/')),exist_ok=True)
    with open(os.path.join(PROJ,'data','metadata','ena_PRJNA713798_runs.tsv'),encoding='utf-8') as f: L=f.read().strip().split('\n')
    row=[l for l in L[1:] if l.startswith('SRR13938890\t')]; T('N1.0 SRR13938890 row found in project ENA table',bool(row))
    with open(os.path.join(R,'data','metadata','ena_PRJNA713798_runs.tsv'),'w',encoding='utf-8') as f: f.write(L[0]+'\n'+row[0]+'\n')
    for p in ['ena_PRJNA740269_runs.tsv','ena_PRJNA1017469_runs.tsv']:
        with open(os.path.join(R,'data','metadata',p),'w',encoding='utf-8') as f: f.write(L[0]+'\n')
    with open(os.path.join(R,'data','metadata','ena_PRJNA1017469_sample_attributes.tsv'),'w') as f: f.write('x\n')
    shutil.copy2(os.path.join(PROJ,'data','reference','amel_OGSv3.2_cds.kidx'),os.path.join(R,'data','reference','amel_OGSv3.2_cds.kidx'))
    def go(cut=None):
        m=load('DEPLOYED',R); m.N=2000; m.refseq_and_mapping=lambda:None; m.kallisto=lambda:KEXE[0]
        if cut is not None:
            real_get=m.get; used={'n':0}
            def get(url,*a,**k):
                resp=real_get(url,*a,**k)
                if 'fastq' in url and used['n']==0:
                    used['n']+=1; data=resp.read(cut*2); resp.close(); return H.Breaking(data,cut)
                return resp
            m.get=get
        m.main()
    c=H.Ctx(base); t0=time.time(); go(cut=60000); st,_=H.state(c)
    T('N1.1 real ENA stream cut after 60 kB -> subset A DOWNLOAD_FAILED (not complete)',st.get(('GSE168740','A','SRR13938890'))=='DOWNLOAD_FAILED',str(st))
    T('N1.2 subset B (reads 2001-4000) streamed from ENA and COMPLETE',st.get(('GSE168740','B','SRR13938890'))=='COMPLETE')
    go(); st,_=H.state(c)
    T('N1.3 restart: subset A COMPLETE; interrupted partial kept as .part.prev',st.get(('GSE168740','A','SRR13938890'))=='COMPLETE' and os.path.exists(os.path.join(R,'data','raw','reads_subsample','SRR13938890_R1_subsetA.fastq.gz.part.prev')))
    cj=os.path.join(R,'results','rawreads','kallisto','GSE168740','A','SRR13938890','_COMPLETE.json')
    if os.path.exists(cj):
        with open(cj,encoding='utf-8') as f: rj=json.load(f)
        with gzip.open(os.path.join(R,rj['fastq']),'rt') as f: first=f.readline().strip()
        T('N1.4 FASTQ begins at read 1 of SRR13938890; n_processed = 2000',re.match(r'^@SRR13938890\.1( |/|$)',first) and rj['n_processed']==2000,first[:50])
        T('N1.5 real brain reads pseudoalign to the CDS index (> 20%)',(rj.get('p_pseudoaligned') or 0)>20,f"{rj.get('p_pseudoaligned')}%")
    out(f'   network test took {time.time()-t0:.0f} s')
def n2():
    """Real ENA: the stream is cut after X bytes and continued with HTTP Range (FQ_MAX_RESUMES=3) in the SAME attempt;
    the reads must be byte-identical to the uninterrupted N1 download of the same run."""
    import harness as H
    ref=os.path.join(SCR,'n1_ena','R','data','raw','reads_subsample','SRR13938890_R1_subsetA.fastq.gz')
    if not T('N2.0 reference FASTQ from N1 present',os.path.exists(ref)): return
    with gzip.open(ref,'rb') as f: want=f.read()
    for cut in (30000,90000):                                    # 30 kB: overlap re-request from byte 0; 90 kB: true 206 from byte 24,464
        base=os.path.join(SCR,f'n2_ena_{cut}'); R=os.path.join(base,'R')
        shutil.copytree(os.path.join(SCR,'n1_ena','R','data','metadata'),os.path.join(R,'data','metadata'))
        for x in ['data/reference','data/raw/BeeBase','logs']: os.makedirs(os.path.join(R,*x.split('/')),exist_ok=True)
        shutil.copy2(os.path.join(PROJ,'data','reference','amel_OGSv3.2_cds.kidx'),os.path.join(R,'data','reference','amel_OGSv3.2_cds.kidx'))
        for k in ACC_ENV: os.environ.pop(k,None)
        os.environ.update({'FQ_MAX_RESUMES':'3','FQ_READ_TIMEOUT':'90'})
        try:
            m=load('DEPLOYED',R); m.N=2000; m.refseq_and_mapping=lambda:None; m.kallisto=lambda:KEXE[0]
            real_get=m.get; used={'n':0}
            def get(url,*a,**k):
                resp=real_get(url,*a,**k)
                if 'fastq' in url and used['n']==0:
                    used['n']+=1; data=resp.read(cut*2); resp.close(); return H.Breaking(data,cut)
                return resp
            m.get=get; m.main()
        finally:
            for k in ACC_ENV: os.environ.pop(k,None)
        c=H.Ctx(base); st,_=H.state(c); rows=H.manifest_rows(c)
        fv=[d for d in rows if d['run']=='SRR13938890' and d['subset']=='A' and d['state']=='FASTQ_VALIDATED']
        T(f'N2.1 [{cut//1000} kB cut] real ENA stream cut, continued via HTTP Range in the same attempt -> COMPLETE (resumes=1)',
          st.get(('GSE168740','A','SRR13938890'))=='COMPLETE' and bool(fv) and fv[0]['note']=='resumes=1' and sum(1 for d in rows if d['run']=='SRR13938890' and d['subset']=='A' and d['state']=='DOWNLOAD_STARTED')==1,
          (fv[0]['note'] if fv else str(st)))
        p=os.path.join(R,'data','raw','reads_subsample','SRR13938890_R1_subsetA.fastq.gz')
        got=gzip.open(p,'rb').read() if os.path.exists(p) else b''
        T(f'N2.2 [{cut//1000} kB cut] reads byte-identical to the uninterrupted N1 download (2,000 reads)',got==want and got.count(b'\n')==8000,f'{len(got)} vs {len(want)} bytes')
if NETWORK:
    section('N1  NETWORK: real ENA stream (SRR13938890, first 4,000 reads), interrupted then resumed - DEPLOYED script',n1)
    section('N2  NETWORK: real ENA byte-range resume (FQ_MAX_RESUMES=3) after a cut at 30 kB and at 90 kB - DEPLOYED script',n2)
else: out(''); out('N1/N2  NETWORK tests skipped (run again with --network to include them)')

# ===================================================================================== SUMMARY
with open(RES,'w',encoding='utf-8') as f:
    f.write('test\tresult\tdetail\n'); [f.write(f'{t}\t{s}\t{d}\n'.replace('\r',' ')) for t,s,d in RESULTS]
npass=sum(r[1]=='PASS' for r in RESULTS); nfail=sum(r[1]=='FAIL' for r in RESULTS); nna=sum(r[1]=='N/A' for r in RESULTS)
net_ok=NETWORK and any(r[0].startswith('N1') for r in RESULTS) and any(r[0].startswith('N2') for r in RESULTS) and all(r[1]=='PASS' for r in RESULTS if r[0].startswith(('N1','N2')))
out(''); out('='*100); out(f'SUMMARY: {npass} PASS, {nfail} FAIL, {nna} N/A   (script under test sha256 {sha(SCRIPTS["DEPLOYED"])[:16]})')
for t,s_,d in RESULTS:
    if s_=='FAIL': out(f'  FAIL: {t}  [{d}]')
for t,s_,d in RESULTS:
    if s_=='N/A': out(f'  N/A : {t}  [{d}]')
mk=os.path.join(PROJ,'preflight_scratch','PREFLIGHT_PASSED.json')
if nfail==0 and net_ok:
    rec={'timestamp':TS,'script':os.path.relpath(SCRIPTS['DEPLOYED'],PROJ),'script_sha256':sha(SCRIPTS['DEPLOYED']),'python':sys.executable,
         'python_version':sys.version.split()[0],'platform':platform.platform(),'network_test_passed':True,'pass':npass,'fail':0,'na':nna,'log':LOG}
    if os.path.exists(mk): os.rename(mk,mk+'.superseded_'+TS)
    with open(mk+'.tmp','w',encoding='utf-8') as f: json.dump(rec,f,indent=1)
    os.replace(mk+'.tmp',mk); out('RESULT: PASS (offline + network). Marker written: '+mk)
else:
    if os.path.exists(mk):
        try:
            with open(mk,encoding='utf-8') as f: old=json.load(f)
        except Exception: old={}
        if old.get('script_sha256')==sha(SCRIPTS['DEPLOYED']) and nfail: os.rename(mk,mk+'.invalidated_'+TS); out('Existing PASS marker for this script invalidated (renamed) because tests now fail.')
    out('RESULT: '+('FAIL' if nfail else 'INCOMPLETE - network test not run')+'  -> the pipeline will refuse to start (no valid PREFLIGHT_PASSED.json)')
out(f'Full log: {LOG}'); out(f'Results table: {RES}')
