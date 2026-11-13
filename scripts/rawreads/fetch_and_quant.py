#!/usr/bin/env python3
"""Raw-read re-quantification (targeted subsampling) - stdlib only, runs on Windows (Anaconda) or Linux.
For each run of PRJNA713798 (GSE168740), PRJNA740269 (GSE178742) and PRJNA1017469:
  1. ENA run + sample metadata are saved to data/metadata/ena_<project>_runs.tsv (BioSample attributes for PRJNA1017469).
  2. The first N reads of the R1 FASTQ are streamed from ENA (no full download) into a temporary fastq.gz.
  3. kallisto (official pachterlab GitHub release) quantifies them single-end against the BeeBase OGS v3.2 CDS
     (index built once). Subsample FASTQs are KEPT (data/raw/reads_subsample/), never deleted by this script.
  4. Validation per run (see process_run): download completes without error into <fastq>.part; strict re-count of
     records on disk equals reads written; SHA-256 of the FASTQ file recorded; kallisto exit code 0, n_processed equals
     reads written, abundance.tsv/run_info.json non-empty. Only then is _COMPLETE.json written (atomically) and the
     output folder moved to its final name. Per-run states: results/rawreads/run_manifest.tsv (append-only).
  5. Positional-bias check: for GSE168740 runs a second, non-overlapping subset (reads N+1..2N, label 'B') is
     quantified separately (subset 'A' = reads 1..N).
  6. Run exclusions (2026-11-12, approved): runs listed in data/metadata/run_exclusions.tsv (PRJNA1017469 duplicate ENA
     listings of the same 47 libraries) are not downloaded or quantified; each is recorded once in the manifest as
     EXCLUDED_DUPLICATE with the kept run, reason and evidence. The script stops if that list is inconsistent.
  7. Acceleration options (2026-11-12; ALL OFF unless set - with none set the behaviour is identical to the previous version):
     FQ_SELECTION=<tsv>       process only the authorized rows of an approved run-selection file, in its priority order
     FQ_WORKERS=1..3          number of simultaneous downloads (kallisto still runs one job at a time, same settings)
     FQ_MAX_RESUMES=n         after a dropped connection, continue the SAME file from the exact byte reached (HTTP Range);
                              the server must answer 206 with Content-Range starting at that byte and an unchanged length /
                              ETag / Last-Modified, otherwise the attempt fails (never spliced). gzip CRC checks still apply.
     FQ_READ_TIMEOUT=s        stall timeout per network read (default 300)
     FQ_JOB_MAX_SECONDS=s     give up on one download after s seconds and move on (recorded DOWNLOAD_FAILED)
     FQ_RETRY_PASSES=n        at the end of the queue, retry jobs whose download failed in this run, n times
     FQ_SLOW_403=1            longer back-off after HTTP 403; with >1 worker, 3 x 403 within 10 min -> 1 download at a time
PRELIMINARY: subsample of Read 1 only; not equivalent to full-depth paired-end quantification.
Resumable: a run is skipped only if results/rawreads/kallisto/<study>/<subset>/<run>/_COMPLETE.json exists and re-validates (see process_run). Log: logs/rawreads.log
Limitation (documented): R1 only, first N reads (default 500,000; NREADS env var) - a subsample, not full-depth quantification."""
import os, sys, gzip, json, hashlib, shutil, subprocess, urllib.request, urllib.error, zipfile, tarfile, io, time, platform, xml.etree.ElementTree as ET
import threading, contextlib, http.client, re as _re, concurrent.futures
ROOT=os.path.abspath(sys.argv[1] if len(sys.argv)>1 else '.'); N=int(os.environ.get('NREADS','500000'))  # reduced from 2,000,000: measured ENA throughput on this connection ~17,600 reads/min
LOG=os.path.join(ROOT,'logs','rawreads.log'); os.makedirs(os.path.dirname(LOG),exist_ok=True)
def _envi(k,d,lo=None,hi=None):
    v=int(os.environ.get(k,str(d)))
    if (lo is not None and v<lo) or (hi is not None and v>hi): raise SystemExit(f'{k}={v} outside the allowed range {lo}..{hi}')
    return v
CFG=dict(workers=_envi('FQ_WORKERS',1,1,3),max_resumes=_envi('FQ_MAX_RESUMES',0,0,10),read_timeout=_envi('FQ_READ_TIMEOUT',300,10,600),
         job_max_s=_envi('FQ_JOB_MAX_SECONDS',0,0,7200),retry_passes=_envi('FQ_RETRY_PASSES',0,0,3),selection=os.environ.get('FQ_SELECTION',''),
         slow403=os.environ.get('FQ_SLOW_403','0')=='1',backoff=float(os.environ.get('FQ_BACKOFF','5')))
_LOCK=threading.RLock(); QUANT_LOCK=threading.Lock()
def log(m):
    s=time.strftime('%Y-%m-%dT%H:%M:%S')+'  '+m
    with _LOCK:
        print(s,flush=True); open(LOG,'a',encoding='utf-8').write(s+'\n')
_SLOT=threading.Condition(); _ACTIVE=[0]; _ALLOWED=[CFG['workers']]; _T403=[]
def _note_403():
    with _SLOT:
        now=time.time(); _T403.append(now); recent=[t for t in _T403 if now-t<600]
        if _ALLOWED[0]>1 and len(recent)>=3:
            _ALLOWED[0]=1; log('  THROTTLE: 3 x HTTP 403 within 10 min -> 1 download at a time for the rest of this run')
@contextlib.contextmanager
def download_slot():
    with _SLOT:
        while _ACTIVE[0]>=_ALLOWED[0]: _SLOT.wait()
        _ACTIVE[0]+=1
    try: yield
    finally:
        with _SLOT: _ACTIVE[0]-=1; _SLOT.notify_all()
def get(url,timeout=120,tries=4):
    for i in range(tries):
        try: return urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Apis-project/1.0'}),timeout=timeout)
        except Exception as e:
            log(f'  retry {i+1} {url}: {e}')
            if CFG['slow403'] and isinstance(e,urllib.error.HTTPError) and e.code==403: _note_403(); time.sleep(6*CFG['backoff']*(i+1))
            else: time.sleep(CFG['backoff']*(i+1))
    raise RuntimeError('failed '+url)
def get_ranged(url,start,timeout):
    """One request for bytes start.. of url (no retry loop here; ResumableStream decides)."""
    return urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Apis-project/1.0','Range':f'bytes={start}-'}),timeout=timeout)
class ResumeRefused(Exception): pass
class ResumableStream:
    """File-like view of one remote file that survives dropped connections by continuing at the exact byte reached.
    A continuation re-requests the last OVERLAP bytes already received and is accepted only if the server answers 206 with
    Content-Range starting exactly there, those bytes are byte-for-byte identical to what was received, and the file's total
    length / ETag / Last-Modified are unchanged. Anything else raises ResumeRefused (the attempt fails; nothing is spliced).
    The gzip layer above still verifies every complete member's CRC-32 and length."""
    OVERLAP=65536
    def __init__(s,url,timeout,max_resumes,deadline,tag):
        s.url,s.timeout,s.max,s.deadline,s.tag=url,timeout,max_resumes,deadline,tag; s.pos=0; s.resumes=0; s.tail=b''
        s.resp=get(url,timeout); h=getattr(s.resp,'headers',None)
        g=(lambda k: h.get(k) if h is not None else None)
        s.total=int(g('Content-Length')) if g('Content-Length') else None; s.etag=g('ETag'); s.lm=g('Last-Modified')
    def _deadline(s):
        if s.deadline and time.time()>s.deadline: raise TimeoutError(f"job time limit reached ({CFG['job_max_s']} s) at byte {s.pos}")
    def read(s,n=-1):
        while True:
            s._deadline()
            try:
                if s.resp is None: s._reopen()             # a drop while reconnecting counts as one more resume
                b=s.resp.read(n)
                if b or n==0 or s.total is None or s.pos>=s.total:
                    s.pos+=len(b); s.tail=(s.tail+b)[-s.OVERLAP:]; return b
                err=ConnectionError(f'connection closed at byte {s.pos} of {s.total}')
            except (OSError,http.client.HTTPException) as e: err=e
            if s.resumes>=s.max: raise err
            s.resumes+=1
            try: s.resp.close()
            except Exception: pass
            log(f"  {s.tag}: resume {s.resumes}/{s.max} at byte {s.pos} after: {str(err)[:120]}"); time.sleep(CFG['backoff']*s.resumes)
            s.resp=None
    def _reopen(s):
        ov=len(s.tail); start=s.pos-ov
        try: r=get_ranged(s.url,start,s.timeout)
        except urllib.error.HTTPError as e:
            if e.code==403 and CFG['slow403']: _note_403()
            raise
        st=getattr(r,'status',None) or r.getcode(); cr=r.headers.get('Content-Range','') or ''
        m=_re.match(r'bytes (\d+)-(\d+)/(\d+|\*)$',cr.strip())
        why=None
        full0=(start==0 and st==200 and s.total is not None and str(s.total)==(r.headers.get('Content-Length') or ''))   # a 200 for bytes=0- is a valid answer
        if not full0 and (st!=206 or not m or int(m.group(1))!=start): why=f'server did not continue at byte {start} (HTTP {st}, Content-Range {cr!r})'
        elif not full0 and s.total is not None and m.group(3)!='*' and int(m.group(3))!=s.total: why=f'file length changed ({m.group(3)} vs {s.total})'
        elif s.etag and r.headers.get('ETag') and r.headers.get('ETag')!=s.etag: why='ETag changed'
        elif s.lm and r.headers.get('Last-Modified') and r.headers.get('Last-Modified')!=s.lm: why='Last-Modified changed'
        if not why and ov:
            got=b''
            while len(got)<ov:
                b=r.read(ov-len(got))
                if not b: raise ConnectionError(f'connection closed during the overlap check ({len(got)} of {ov} bytes)')
                got+=b
            if got!=s.tail: why=f'the {ov} overlapping bytes before byte {s.pos} differ from those already received'
        if why:
            try: r.close()
            except Exception: pass
            raise ResumeRefused('resume refused: '+why)
        s.resp=r
    def readable(s): return True
    def close(s):
        try: s.resp.close()
        except Exception: pass
    def __enter__(s): return s
    def __exit__(s,*a): s.close()
PROJECTS=[('PRJNA713798','GSE168740'),('PRJNA740269','GSE178742'),('PRJNA1017469','PRJNA1017469')]
FIELDS='run_accession,sample_accession,sample_alias,sample_title,experiment_alias,experiment_title,library_layout,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes,instrument_model'
def ena_runs(p):
    out=os.path.join(ROOT,'data','metadata',f'ena_{p}_runs.tsv')
    if not os.path.exists(out):
        t=get(f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={p}&result=read_run&fields={FIELDS}&format=tsv&limit=1000').read().decode()
        open(out,'w',encoding='utf-8').write(t); log(f'saved {out}')
    L=open(out,encoding='utf-8').read().strip().split('\n'); h=L[0].split('\t'); return [dict(zip(h,l.split('\t'))) for l in L[1:]]
def biosample_attrs(rows,p):
    out=os.path.join(ROOT,'data','metadata',f'ena_{p}_sample_attributes.tsv')
    if os.path.exists(out): return
    recs=[];keys=[]
    for r in rows:
        try:
            x=ET.fromstring(get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{r['sample_accession']}").read())
            a={e.findtext('TAG'):e.findtext('VALUE') for e in x.iter('SAMPLE_ATTRIBUTE')}
        except Exception as e: a={'error':str(e)}
        a={'run_accession':r['run_accession'],'sample_accession':r['sample_accession'],**a}; recs.append(a); keys+= [k for k in a if k not in keys]
    with open(out,'w',encoding='utf-8') as f:
        f.write('\t'.join(keys)+'\n'); [f.write('\t'.join(str(a.get(k,'')).replace('\t',' ') for k in keys)+'\n') for a in recs]
    log(f'saved {out}')
def kallisto():
    d=os.path.join(ROOT,'tools','kallisto'); os.makedirs(d,exist_ok=True)
    exe=[os.path.join(dp,f) for dp,_,fs in os.walk(d) for f in fs if f in ('kallisto.exe','kallisto')]
    if exe: return exe[0]
    win=platform.system()=='Windows'
    rel=json.load(get('https://api.github.com/repos/pachterlab/kallisto/releases?per_page=40'))
    for r in rel:
        for a in r['assets']:
            n=a['name'].lower()
            if (win and 'windows' in n) or (not win and 'linux' in n and 'arm' not in n):
                log(f"kallisto release {r['tag_name']} asset {a['name']}"); b=get(a['browser_download_url'],600).read()
                open(os.path.join(d,a['name']),'wb').write(b); open(os.path.join(ROOT,'logs','tool_versions.txt'),'a').write(f"kallisto {r['tag_name']} {a['name']} sha256={hashlib.sha256(b).hexdigest()}\n")
                base=os.path.realpath(d)
                def safe(name):
                    t=os.path.realpath(os.path.join(d,name))
                    if not (t==base or t.startswith(base+os.sep)): raise RuntimeError('unsafe path in archive: '+name)
                if n.endswith('.zip'):
                    z=zipfile.ZipFile(io.BytesIO(b)); [safe(x) for x in z.namelist()]; z.extractall(d)
                else:
                    t=tarfile.open(fileobj=io.BytesIO(b)); ms=[m for m in t.getmembers() if not (m.issym() or m.islnk() or m.isdev())]; [safe(m.name) for m in ms]; t.extractall(d,members=ms)
                exe=[os.path.join(dp,f) for dp,_,fs in os.walk(d) for f in fs if f in ('kallisto.exe','kallisto')]
                if not win: os.chmod(exe[0],0o755)
                return exe[0]
    raise RuntimeError('no kallisto asset found')
def refseq_and_mapping():
    """Downloads the RefSeq proteins used for sequence-level mapping verification and runs verify_ogs32_mapping.py on local files."""
    d=os.path.join(ROOT,'data','raw','NCBI'); os.makedirs(d,exist_ok=True); fa=os.path.join(d,'target_refseq_proteins.fa'); g2p=os.path.join(d,'gene_to_refseq_protein.tsv')
    src=os.path.join(ROOT,'data','metadata','gene_to_refseq_protein.tsv')
    if not os.path.exists(g2p): shutil.copy(src,g2p)
    if not os.path.exists(fa):
        acc=[l.split('\t')[1].strip() for l in open(src).read().strip().split('\n')[1:]]
        b=get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id='+','.join(acc),300).read()
        open(fa,'wb').write(b); log(f'saved {fa} sha256={hashlib.sha256(b).hexdigest()} n_seq={b.count(b">")}')
    r=subprocess.run([sys.executable,os.path.join(ROOT,'scripts','gene_mapping','verify_ogs32_mapping.py'),ROOT],capture_output=True,text=True)
    log('mapping verification (local files): '+(r.stdout+r.stderr)[-600:].replace('\n',' | '))

# ---------------- checkpoint / manifest ----------------
# Persistent per-run manifest: results/rawreads/run_manifest.tsv (append-only event log; the LAST row per
# (study, subset, run) is that run's current state). A run is COMPLETE only when its folder contains _COMPLETE.json
# (written atomically LAST) and re-validation of that record succeeds. Nothing is ever deleted: stale/partial files and
# folders are renamed with a .prev/.failed suffix.
MANIFEST_COLS=['time','study','subset','run','state','reads_target','reads_written','fastq_eof','fastq_sha256','fastq_bytes','reads_recounted','kallisto_rc','n_processed','abundance_rows','note']
def manifest_path(): return os.path.join(ROOT,'results','rawreads','run_manifest.tsv')
def manifest_add(study,sub,run,state,**kw):
    with _LOCK: _manifest_add(study,sub,run,state,**kw)
def _manifest_add(study,sub,run,state,**kw):
    p=manifest_path(); os.makedirs(os.path.dirname(p),exist_ok=True); new=not os.path.exists(p)
    row={'time':time.strftime('%Y-%m-%dT%H:%M:%S'),'study':study,'subset':sub,'run':run,'state':state,**kw}
    with open(p,'a',encoding='utf-8') as f:
        if new: f.write('\t'.join(MANIFEST_COLS)+'\n')
        f.write('\t'.join(str(row.get(c,'')).replace('\t',' ').replace('\n',' ') for c in MANIFEST_COLS)+'\n'); f.flush(); os.fsync(f.fileno())
def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def count_reads_strict(p):
    lines=0
    with gzip.open(p,'rb') as f:
        for _ in f: lines+=1
    if lines%4: raise ValueError(f'{lines} lines, not a multiple of 4')
    return lines//4
def safe_rename(src,dst,tries=10):
    """os.rename that never overwrites and retries briefly when Windows/OneDrive/antivirus holds a transient lock."""
    for i in range(tries):
        if os.path.exists(dst): raise FileExistsError(dst)
        try: os.rename(src,dst); return
        except PermissionError as e:
            if i==tries-1: raise
            log(f'  rename locked ({e}); retry {i+1}'); time.sleep(2*(i+1))
def keep_aside(path,suffix):
    """Rename (never delete) an existing file/folder to <path><suffix>[n]."""
    if not os.path.exists(path): return None
    t=path+suffix; i=1
    while os.path.exists(t): t=f'{path}{suffix}{i}'; i+=1
    safe_rename(path,t); return t
def is_complete(od,fq):
    """True only if _COMPLETE.json exists AND the recorded FASTQ sha256/size and the outputs still validate."""
    c=os.path.join(od,'_COMPLETE.json')
    if not os.path.exists(c): return False,'no _COMPLETE.json'
    try:
        with open(c,encoding='utf-8') as f: rec=json.load(f)
    except Exception as e: return False,'unreadable _COMPLETE.json: '+str(e)
    for f in ('abundance.tsv','run_info.json'):
        fp=os.path.join(od,f)
        if not os.path.exists(fp) or os.path.getsize(fp)==0: return False,f'missing/empty {f}'
    if not os.path.exists(fq): return False,'FASTQ missing'
    if os.path.getsize(fq)!=rec.get('fastq_bytes'): return False,'FASTQ size differs from record'
    if sha256_file(fq)!=rec.get('fastq_sha256'): return False,'FASTQ sha256 differs from record'
    return True,'ok'
def process_run(K,idx,tmpd,summ,label,sub,r):
    run=r['run_accession']; od=os.path.join(ROOT,'results','rawreads','kallisto',label,sub,run)
    fq=os.path.join(tmpd,f'{run}_R1_subset{sub}.fastq.gz'); part=fq+'.part'
    ok,why=is_complete(od,fq)
    if ok: return                                  # validated -> skip, touch nothing
    if os.path.exists(os.path.join(od,'_COMPLETE.json')):
        manifest_add(label,sub,run,'REVALIDATION_FAILED',note=why); log(f'  {label} {sub} {run} REVALIDATION_FAILED ({why}) - left untouched for review'); return
    url='https://'+r['fastq_ftp'].split(';')[0]; skip=0 if sub=='A' else N
    # a FASTQ that already passed validation (manifest FASTQ_VALIDATED with matching size + sha256) is reused, not re-downloaded
    v=last_validated_fastq(label,sub,run)
    if v and os.path.exists(fq) and os.path.getsize(fq)==int(v['fastq_bytes']) and int(v['reads_target'])==N and sha256_file(fq)==v['fastq_sha256']:
        for pth,suf in ((part,'.prev'),(od,'.failed')):
            k=keep_aside(pth,suf)
            if k: log(f'  {run}: kept earlier incomplete {os.path.basename(pth)} as {os.path.basename(k)}')
        log(f'  {label} {sub} {run} reusing validated FASTQ (sha256 matches manifest)')
        return quantify(K,idx,summ,label,sub,r,run,od,fq,url,int(v['reads_written']),v['fastq_eof']=='True',int(v['reads_recounted']),v['fastq_sha256'],int(v['fastq_bytes']),'')
    # anything else left from an interrupted attempt is kept aside, not reused, not deleted
    for pth,suf in ((part,'.prev'),(fq,'.prev'),(od,'.failed')):
        k=keep_aside(pth,suf)
        if k: log(f'  {run}: kept earlier incomplete {os.path.basename(pth)} as {os.path.basename(k)}')
    resumable=CFG['max_resumes']>0 or CFG['job_max_s']>0; src=None
    with download_slot():
      manifest_add(label,sub,run,'DOWNLOAD_STARTED',reads_target=N,note=url)
      n=0; seen=0; eof=False; md5=hashlib.md5(); err=''
      try:
        src=(ResumableStream(url,CFG['read_timeout'],CFG['max_resumes'],time.time()+CFG['job_max_s'] if CFG['job_max_s'] else 0,f'{label} {sub} {run}')
             if resumable else get(url,CFG['read_timeout']))
        with src as resp, gzip.GzipFile(fileobj=resp) as gz, gzip.open(part,'wb',compresslevel=6) as out:
            while n<N:
                rec=[gz.readline() for _ in range(4)]
                if not rec[0]: eof=True; break
                if not rec[3]: raise ValueError('truncated FASTQ record in source stream')
                if resumable and not (rec[0][:1]==b'@' and rec[2][:1]==b'+' and len(rec[1].rstrip(b'\r\n'))==len(rec[3].rstrip(b'\r\n'))>0):
                    raise ValueError(f'malformed FASTQ record {seen+1} in source stream')
                seen+=1
                if seen<=skip: continue
                b=b''.join(rec); out.write(b); md5.update(b); n+=1
      except Exception as e: err=str(e)[:200]
    rn=(f'resumes={src.resumes}' if resumable and src is not None else '')
    if err or (n<N and not eof) or n==0:
        manifest_add(label,sub,run,'DOWNLOAD_FAILED',reads_target=N,reads_written=n,fastq_eof=eof,note=((err or 'incomplete')+(f'; {rn}' if rn else '')))
        log(f'  {label} {sub} {run} DOWNLOAD_FAILED reads={n} ({err or "incomplete"}); partial kept as .part'); return
    rc_n=count_reads_strict(part)
    if rc_n!=n:
        manifest_add(label,sub,run,'READCOUNT_FAILED',reads_written=n,reads_recounted=rc_n); log(f'  {run} READCOUNT_FAILED'); return
    safe_rename(part,fq); fsha=sha256_file(fq); fbytes=os.path.getsize(fq)
    manifest_add(label,sub,run,'FASTQ_VALIDATED',reads_target=N,reads_written=n,fastq_eof=eof,fastq_sha256=fsha,fastq_bytes=fbytes,reads_recounted=rc_n,note=rn)
    return quantify(K,idx,summ,label,sub,r,run,od,fq,url,n,eof,rc_n,fsha,fbytes,md5.hexdigest())
def last_validated_fastq(label,sub,run):
    if not os.path.exists(manifest_path()): return None
    with _LOCK:
        with open(manifest_path(),encoding='utf-8') as f: rows=[l.rstrip('\n').split('\t') for l in f]
    h=rows[0]; last=None
    for x in rows[1:]:
        d=dict(zip(h,x))
        if (d.get('study'),d.get('subset'),d.get('run'))==(label,sub,run):
            if d.get('state')=='FASTQ_VALIDATED': last=d
            elif d.get('state')=='DOWNLOAD_STARTED': last=None
    return last
def read_run_info(p):
    """Parse kallisto run_info.json. kallisto 0.51.1 on Windows writes the "call" field with raw, unescaped backslashes
    (e.g. "C:\\Users\\..."), which is invalid JSON. Standard parsing is tried first; on failure the numeric fields are
    extracted individually by exact key. Returns (dict, how)."""
    with open(p,encoding='utf-8',errors='replace') as f: txt=f.read()
    try: return json.loads(txt),'json'
    except ValueError: pass
    import re
    info={}
    for k in ('n_targets','n_bootstraps','n_processed','n_pseudoaligned','n_unique','index_version'):
        m=re.search(r'"%s"\s*:\s*(-?\d+)\s*[,}\n]'%re.escape(k),txt)
        if m: info[k]=int(m.group(1))
    for k in ('p_pseudoaligned','p_unique'):
        m=re.search(r'"%s"\s*:\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*[,}\n]'%re.escape(k),txt)
        if m: info[k]=float(m.group(1))
    m=re.search(r'"kallisto_version"\s*:\s*"([^"\\]*)"',txt)
    if m: info['kallisto_version']=m.group(1)
    return info,'regex-fallback (invalid JSON from kallisto)'
def quantify(K,idx,summ,label,sub,r,run,od,fq,url,n,eof,rc_n,fsha,fbytes,md5hex):
    tmpod=od+'.inprogress'; keep_aside(tmpod,'.failed'); os.makedirs(tmpod)
    cmd=[K,'quant','-i',idx,'-o',tmpod,'--single','-l','200','-s','30','-t','2',fq]
    with QUANT_LOCK: q=subprocess.run(cmd,capture_output=True,text=True,errors='replace')   # one kallisto at a time
    with open(os.path.join(tmpod,'kallisto_stderr.txt'),'w',encoding='utf-8') as f: f.write(q.stderr or '')
    info={}; ab=os.path.join(tmpod,'abundance.tsv'); rows=0
    how='missing'
    try: info,how=read_run_info(os.path.join(tmpod,'run_info.json'))
    except FileNotFoundError: pass
    except Exception as e: how='unreadable: '+str(e)[:100]
    if os.path.exists(ab):
        with open(ab,encoding='utf-8') as f: rows=sum(1 for _ in f)-1
    ri=os.path.join(tmpod,'run_info.json'); good=(q.returncode==0 and info.get('n_processed')==n and isinstance(info.get('n_pseudoaligned'),int) and rows>0 and os.path.exists(ri) and os.path.getsize(ri)>0)
    if not good:
        manifest_add(label,sub,run,'QUANT_FAILED',reads_written=n,fastq_sha256=fsha,kallisto_rc=q.returncode,n_processed=info.get('n_processed',''),abundance_rows=rows,note=f'run_info: {how}; '+(q.stderr or '')[-200:])
        log(f'  {label} {sub} {run} QUANT_FAILED rc={q.returncode}; output kept in {os.path.basename(tmpod)}'); return
    rec=dict(study=label,subset=sub,run=run,sample_alias=r.get('sample_alias',''),reads_target=N,reads_written=n,fastq_eof=eof,reads_recounted=rc_n,
             fastq=os.path.relpath(fq,ROOT),fastq_sha256=fsha,fastq_bytes=fbytes,source_url=url,kallisto_command=' '.join([os.path.basename(cmd[0])]+cmd[1:]),
             kallisto_rc=q.returncode,n_processed=info.get('n_processed'),n_pseudoaligned=info.get('n_pseudoaligned'),p_pseudoaligned=info.get('p_pseudoaligned'),
             abundance_rows=rows,run_info_parsed_with=how,completed=time.strftime('%Y-%m-%dT%H:%M:%S'))
    tmpc=os.path.join(tmpod,'_COMPLETE.json.tmp')
    with open(tmpc,'w',encoding='utf-8') as f: json.dump(rec,f,indent=1); f.flush(); os.fsync(f.fileno())
    safe_rename(tmpc,os.path.join(tmpod,'_COMPLETE.json'))
    os.makedirs(os.path.dirname(od),exist_ok=True); safe_rename(tmpod,od)   # od does not exist (kept aside above)
    manifest_add(label,sub,run,'COMPLETE',reads_target=N,reads_written=n,fastq_eof=eof,fastq_sha256=fsha,fastq_bytes=fbytes,reads_recounted=rc_n,kallisto_rc=0,n_processed=info.get('n_processed'),abundance_rows=rows,note='run_info: '+how)
    with _LOCK, open(summ,'a',encoding='utf-8') as f: f.write(f"{label}\t{sub}\t{run}\t{r.get('sample_alias','')}\t{n}\t{info.get('n_processed','')}\t{info.get('n_pseudoaligned','')}\t{info.get('p_pseudoaligned','')}\t{md5hex}\tTrue\tOK\t{fq}\t{url}\t{rec['kallisto_command']}\n")
    log(f"  {label} {sub} {run} COMPLETE reads={n} pseudoaligned={info.get('p_pseudoaligned','?')}%")
def preflight_gate():
    """Refuse to run unless a full Windows preflight (offline + network) PASSED for this exact script file.
    The marker preflight_scratch/PREFLIGHT_PASSED.json is written only by scripts/validation/windows_preflight.py."""
    if os.environ.get('FQ_TEST_HARNESS')=='1': return          # set only by the test harness (scratch folders)
    mk=os.path.join(ROOT,'preflight_scratch','PREFLIGHT_PASSED.json'); me=sha256_file(os.path.abspath(__file__))
    try:
        with open(mk,encoding='utf-8') as f: rec=json.load(f)
    except Exception: rec={}
    if rec.get('script_sha256')!=me or not rec.get('network_test_passed'):
        log(f'REFUSING TO START: no passing Windows preflight recorded for this script (sha256 {me[:16]}). Run scripts/validation/windows_preflight.py <project> --network first.')
        raise SystemExit(2)
    log(f"preflight gate OK: {rec.get('timestamp')} sha256 {me[:16]} python {rec.get('python')}")
def manifest_rows_all():
    p=manifest_path()
    if not os.path.exists(p): return []
    with _LOCK:
        with open(p,encoding='utf-8') as f: L=f.read().splitlines()
    if not L: return []
    h=L[0].split('\t'); return [dict(zip(h,l.split('\t'))) for l in L[1:]]
def load_exclusions():
    """data/metadata/run_exclusions.tsv (approved audit trail): runs listed there are NOT processed. Each row names the run kept
    instead. Absent file = no exclusions. Rows starting with '#' are comments. Returns {(study, run): row}."""
    p=os.path.join(ROOT,'data','metadata','run_exclusions.tsv')
    if not os.path.exists(p): return {}
    with open(p,encoding='utf-8') as f: L=[l for l in f.read().splitlines() if l.strip() and not l.startswith('#')]
    h=L[0].split('\t'); out={}
    for l in L[1:]:
        d=dict(zip(h,l.split('\t')))
        if d.get('study') and d.get('run_accession'): out[(d['study'],d['run_accession'])]=d
    return out
def apply_exclusions(jobs,excl,allrows):
    """Removes excluded runs from the queue. Stops (SystemExit) instead of guessing if the exclusion list is inconsistent:
    a kept run that is absent from its study's ENA table, or a kept run that is itself excluded."""
    for (study,run),d in excl.items():
        runs={r['run_accession'] for r in allrows.get(study,[])}
        if run not in runs: raise SystemExit(f'run_exclusions.tsv: {study} {run} is not in the ENA table - refusing to continue')
        k=d.get('kept_run','')
        if k and (k not in runs or (study,k) in excl): raise SystemExit(f'run_exclusions.tsv: kept run {k} for {run} is missing or itself excluded - refusing to continue')
    kept=[j for j in jobs if (j[0],j[2]['run_accession']) not in excl]
    return kept,[j for j in jobs if (j[0],j[2]['run_accession']) in excl]
def job_digest(jobs,study):
    """SHA-256 of the ordered (study, subset, run, fastq_ftp) list of one study - proves other studies' queues are unchanged."""
    h=hashlib.sha256()
    for label,sub,r in jobs:
        if label==study: h.update(f"{label}\t{sub}\t{r['run_accession']}\t{r.get('fastq_ftp','')}\n".encode())
    return h.hexdigest()
def load_selection(path,jobs):
    """Approved run-selection file (FQ_SELECTION). Columns: priority, study, subset, run, stage, authorized (YES/NO), ... .
    Every row must be a job of the approved queue (after exclusions); duplicates or unknown jobs stop the script.
    Returns the authorized jobs in priority order, the selection rows, and the file's sha256."""
    if not os.path.isabs(path): path=os.path.join(ROOT,path)
    with open(path,encoding='utf-8') as f: txt=f.read()
    L=[l for l in txt.splitlines() if l.strip() and not l.startswith('#')]; h=L[0].split('\t'); rows=[dict(zip(h,l.split('\t'))) for l in L[1:]]
    byk={(j[0],j[1],j[2]['run_accession']):j for j in jobs}; seen=set(); out=[]
    for d in rows:
        k=(d.get('study'),d.get('subset'),d.get('run'))
        if k not in byk: raise SystemExit(f'run selection: {k} is not a job of the approved queue (excluded duplicate, unknown run or wrong subset) - refusing to continue')
        if k in seen: raise SystemExit(f'run selection: {k} listed twice - refusing to continue')
        if d.get('authorized') not in ('YES','NO'): raise SystemExit(f"run selection: {k} authorized must be YES or NO, not {d.get('authorized')!r}")
        seen.add(k)
    for d in sorted(rows,key=lambda d:int(d['priority'])):
        if d['authorized']=='YES': out.append(byk[(d['study'],d['subset'],d['run'])])
    return out,rows,hashlib.sha256(txt.encode('utf-8')).hexdigest()
def _safe_process(K,idx,tmpd,summ,job):
    label,sub,r=job
    try: process_run(K,idx,tmpd,summ,label,sub,r)
    except Exception as e:
        manifest_add(label,sub,r['run_accession'],'ERROR',note=str(e)[:300]); log(f"  {label} {sub} {r['run_accession']} ERROR {e}")
def run_queue(K,idx,tmpd,summ,queue):
    if CFG['workers']==1:
        for job in queue: _safe_process(K,idx,tmpd,summ,job)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=CFG['workers']) as ex:
            list(ex.map(lambda j: _safe_process(K,idx,tmpd,summ,j),queue))
def build_jobs(allrows,exp):
    g78=allrows.get('GSE178742',[]); jobs=[]
    jobs+= [('GSE168740','A',r) for r in allrows.get('GSE168740',[])]
    jobs+= [('GSE178742','A',r) for r in g78 if exp.get(r.get('sample_alias',''),'').startswith('Apiary')]
    jobs+= [('GSE178742','A',r) for r in g78 if not exp.get(r.get('sample_alias',''),'').startswith('Apiary')]
    jobs+= [('GSE168740','B',r) for r in allrows.get('GSE168740',[])]
    jobs+= [('PRJNA1017469','A',r) for r in allrows.get('PRJNA1017469',[])]
    return jobs
def main():
    log(f'start N={N} platform={platform.platform()} python={sys.version.split()[0]}')
    preflight_gate()
    try: refseq_and_mapping()
    except Exception as e: log('mapping verification failed: '+str(e))
    K=kallisto(); v=subprocess.run([K,'version'],capture_output=True,text=True).stdout.strip(); log('kallisto: '+v); open(os.path.join(ROOT,'logs','tool_versions.txt'),'a').write(v+'\n')
    cds=os.path.join(ROOT,'data','raw','BeeBase','amel_OGSv3.2_cds.fa.gz'); idx=os.path.join(ROOT,'data','reference','amel_OGSv3.2_cds.kidx')
    os.makedirs(os.path.dirname(idx),exist_ok=True)
    if not os.path.exists(idx):
        if not os.path.exists(cds): raise SystemExit('missing '+cds+' - run the download step first')
        r=subprocess.run([K,'index','-i',idx,cds],capture_output=True,text=True); log('index cmd: kallisto index -i '+idx+' '+cds+' | '+(r.stderr or '')[-300:])
    tmpd=os.path.join(ROOT,'data','raw','reads_subsample'); os.makedirs(tmpd,exist_ok=True)
    summ=os.path.join(ROOT,'results','rawreads','subsample_log.tsv'); os.makedirs(os.path.dirname(summ),exist_ok=True)
    if not os.path.exists(summ): open(summ,'w').write('project\tsubset\trun\tsample_alias\treads_written\tkallisto_n_processed\tn_pseudoaligned\tp_pseudoaligned\tfastq_md5\tmd5_reverified\tstatus\tfastq_path\tsource_url\tkallisto_command\n')
    allrows={}
    for p,label in PROJECTS:
        rows=ena_runs(p); allrows[label]=rows
        if p=='PRJNA1017469':
            try: biosample_attrs(rows,p)
            except Exception as e: log('biosample attributes failed: '+str(e))
    meta=os.path.join(ROOT,'data','metadata','GSE178742_sample_metadata_verified.tsv'); exp={}
    if os.path.exists(meta):
        L=open(meta,encoding='utf-8').read().strip().split('\n'); h=L[0].split('\t')
        for l in L[1:]: d=dict(zip(h,l.split('\t'))); exp[d['gsm']]=d['experiment']
    jobs=build_jobs(allrows,exp)
    excl=load_exclusions(); jobs,skipped=apply_exclusions(jobs,excl,allrows)
    done_excl={(d['study'],d['run']) for d in manifest_rows_all() if d.get('state')=='EXCLUDED_DUPLICATE'}
    for label,sub,r in skipped:
        if (label,r['run_accession']) not in done_excl:   # one audit row per excluded run, not one per restart
            d=excl[(label,r['run_accession'])]
            manifest_add(label,sub,r['run_accession'],'EXCLUDED_DUPLICATE',note=f"kept {d.get('kept_run','')}; {d.get('reason','')}; {d.get('evidence','')}; approved: {d.get('decided_by','')}")
    for st in ('GSE168740','GSE178742','PRJNA1017469'): log(f'job-list sha256 {st}: {job_digest(jobs,st)}')
    log(f'job queue: {len(jobs)} quantifications after {len(skipped)} listed exclusions (priority: GSE168740 A, GSE178742 apiary, GSE178742 field, GSE168740 B, PRJNA1017469); N={N}')
    if any(CFG[k] for k in ('selection','max_resumes','job_max_s','retry_passes','slow403')) or CFG['workers']>1 or CFG['read_timeout']!=300:
        log('acceleration config: '+', '.join(f'{k}={v}' for k,v in CFG.items()))
    queue=jobs
    if CFG['selection']:
        queue,selrows,selsha=load_selection(CFG['selection'],jobs)
        log(f"run selection {CFG['selection']} sha256 {selsha}: {len(selrows)} rows, {len(queue)} authorized; "
            f"{len(jobs)-len(selrows)} approved jobs not in this selection are left untouched (not attempted, not marked)")
        h=hashlib.sha256()
        for label,sub,r in queue: h.update(f"{label}\t{sub}\t{r['run_accession']}\n".encode())
        log(f'authorized queue sha256: {h.hexdigest()}')
    run_queue(K,idx,tmpd,summ,queue)
    for p in range(CFG['retry_passes']):
        last={}
        for d in manifest_rows_all(): last[(d.get('study'),d.get('subset'),d.get('run'))]=d.get('state')
        again=[j for j in queue if last.get((j[0],j[1],j[2]['run_accession'])) in ('DOWNLOAD_FAILED','ERROR')]
        if not again: break
        log(f'retry pass {p+1}/{CFG["retry_passes"]}: {len(again)} jobs whose download failed in this run'); time.sleep(12*CFG['backoff'])
        run_queue(K,idx,tmpd,summ,again)
    log('DONE')
if __name__=='__main__': main()
