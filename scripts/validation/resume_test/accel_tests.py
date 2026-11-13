"""Tests for the acceleration options of fetch_and_quant.py (resume, parallel downloads, fail-fast, retry pass, run selection).
Portable (Windows + Linux). Everything is written under <base>; a project folder, if given, is only READ.
Usage: python accel_tests.py <empty_base_folder> [<project_dir>]     (exit code 1 on any failure)
The mock ENA server serves local multi-member gzip FASTQs (like ENA's), honours HTTP Range requests and can inject faults:
connection reset, silent early close, read timeout, HTTP 403, ignoring Range, lying about the resume offset, changed ETag/length."""
import os, sys, gzip, json, random, socket, hashlib, importlib.util, subprocess, time, threading, urllib.error, io, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness as H
def J(*p): return os.path.join(*p)
SCRIPT = os.environ.get('FQ_SCRIPT') or os.path.normpath(J(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'rawreads', 'fetch_and_quant.py'))
NREADS = 1500                      # N used by the pipeline in these tests (subset A = reads 1..N, B = N+1..2N)

# ---------------- synthetic multi-member gzip sources ----------------
def make_source(path, acc, nreads, member=400, seed=0):
    rnd = random.Random(f'{acc}{seed}'); recs = []
    for i in range(nreads):
        s = ''.join(rnd.choice('ACGT') for _ in range(76)); q = ''.join(rnd.choice('?@ABCDEFGHI') for _ in range(76))
        recs.append(f'@{acc}.{i+1} {i+1}/1\n{s}\n+\n{q}\n'.encode())
    with open(path, 'wb') as f:
        for k in range(0, nreads, member): f.write(gzip.compress(b''.join(recs[k:k+member]), compresslevel=6))
    return recs

# ---------------- mock ENA server ----------------
class Resp:
    def __init__(s, data, start, status, headers, fault, sleep=0.0):
        s.d = data; s.p = start; s.status = status; s.headers = headers; s.fault = fault; s.sent = 0; s.sleep = sleep
    def getcode(s): return s.status
    def read(s, n=-1):
        if s.sleep: time.sleep(s.sleep)
        n = 8192 if n is None or n < 0 else min(n, 8192)
        if s.fault:
            kind, at = s.fault
            if s.sent >= at:
                if kind == 'reset': raise ConnectionResetError('mock: connection reset by peer')
                if kind == 'timeout': raise socket.timeout('mock: The read operation timed out')
                if kind == 'close': return b''
            n = min(n, at - s.sent)
        out = s.d[s.p:s.p+n]; s.p += len(out); s.sent += len(out); return out
    def readable(s): return True
    def close(s): pass
    def __enter__(s): return s
    def __exit__(s, *a): pass
class Server:
    """plan[acc] = list of per-connection dicts (connection i uses plan[acc][i] if present):
       fault=('reset'|'close'|'timeout', bytes_after_which), http403=True, ignore_range=True, lie_shift=k, etag='x', length=+k, sleep=s"""
    def __init__(s, src_dir, plan=None):
        s.src = src_dir; s.plan = plan or {}; s.conn = {}; s.lock = threading.Lock(); s.log = []
    def _spec(s, acc):
        with s.lock:
            i = s.conn.get(acc, 0); s.conn[acc] = i + 1
        p = s.plan.get(acc, []); return (p[i] if i < len(p) else {}), i
    def _open(s, url, start, ranged=False):
        path = url.split('://LOCAL', 1)[1]; acc = os.path.basename(path).split('.')[0]
        with open(path, 'rb') as f: data = f.read()
        spec, i = s._spec(acc); s.log.append((acc, i, start))
        if spec.get('http403'): raise urllib.error.HTTPError(url, 403, 'Forbidden', {}, None)
        total = len(data) + spec.get('length', 0); hd = {'Content-Length': str(total - start), 'ETag': spec.get('etag', '"v1"'), 'Last-Modified': 'Mon, 01 Jan 2024 00:00:00 GMT'}
        if not ranged or spec.get('ignore_range'):
            hd['Content-Length'] = str(total); return Resp(data, 0, 200, hd, spec.get('fault'), spec.get('sleep', 0))
        served = start - spec.get('lie_shift', 0)
        hd['Content-Range'] = f'bytes {start}-{total-1}/{total}'
        return Resp(data, served, 206, hd, spec.get('fault'), spec.get('sleep', 0))
    def urlopen(s, req, timeout=None):
        """Stands in for urllib.request.urlopen: the pipeline's own get(), get_ranged() and Range checks run unchanged."""
        url = req.full_url if hasattr(req, 'full_url') else req
        if '://LOCAL' not in url: raise urllib.error.URLError('no network in mock test: ' + url)
        rg = req.get_header('Range') if hasattr(req, 'get_header') else None
        start = int(rg.split('=')[1].rstrip('-')) if rg else 0
        return s._open(url, start, rg is not None)
M_CUR = [None]

# ---------------- running the pipeline in a scratch project ----------------
def setup(base, runs, extra_rows=None):
    c = H.Ctx(base); R = c.R; src = c.src
    for d in ['data/metadata', 'data/raw/BeeBase', 'data/reference', 'logs', 'data/raw/reads_subsample']: os.makedirs(J(R, *d.split('/')), exist_ok=True)
    os.makedirs(src, exist_ok=True); recs = {}
    for acc, n in runs: recs[acc] = make_source(J(src, f'{acc}.fq.gz'), acc, n)
    hdr = 'run_accession\tsample_alias\tsample_accession\tfastq_ftp\n'
    gse = [a for a, _ in runs if a.startswith('G')]; prj = [a for a, _ in runs if a.startswith('P')]
    with open(J(R, 'data/metadata/ena_PRJNA713798_runs.tsv'), 'w', encoding='utf-8') as f: f.write(hdr + ''.join(f'{a}\tGSM{a}\tS\tLOCAL{J(src, a + ".fq.gz")}\n' for a in gse))
    with open(J(R, 'data/metadata/ena_PRJNA740269_runs.tsv'), 'w', encoding='utf-8') as f: f.write(hdr)
    with open(J(R, 'data/metadata/ena_PRJNA1017469_runs.tsv'), 'w', encoding='utf-8') as f: f.write(hdr + ''.join(f'{a}\tS{a}\tS\tLOCAL{J(src, a + ".fq.gz")}\n' for a in prj))
    for p in ['data/metadata/ena_PRJNA1017469_sample_attributes.tsv', 'data/raw/BeeBase/amel_OGSv3.2_cds.fa.gz', 'data/reference/amel_OGSv3.2_cds.kidx']:
        with open(J(R, p), 'w') as f: f.write('x\n')
    with open(J(c.base, 'fake_kallisto.py'), 'w') as f: f.write(FAKE)
    return c, recs
FAKE = H.FAKE_SRC.replace("od=a[a.index('-o')+1]; fq=a[-1]; mode=os.environ.get('FAKEK','ok')",
    "od=a[a.index('-o')+1]; fq=a[-1]; mode=os.environ.get('FAKEK','ok')\nlg=os.environ.get('FAKEK_LOG')\nif lg: open(lg,'a').write('start %.6f %s\\n'%(time.time(),os.getpid()))\ntime.sleep(float(os.environ.get('FAKEK_SLEEP','0')))") \
    + "\nif os.environ.get('FAKEK_LOG'): open(os.environ['FAKEK_LOG'],'a').write('end %.6f %s\\n'%(time.time(),os.getpid()))\n"
ENVK = ['FQ_WORKERS', 'FQ_MAX_RESUMES', 'FQ_READ_TIMEOUT', 'FQ_JOB_MAX_SECONDS', 'FQ_RETRY_PASSES', 'FQ_SELECTION', 'FQ_SLOW_403', 'FQ_BACKOFF']
def load(R, env):
    for k in ENVK: os.environ.pop(k, None)
    os.environ.update({k: str(v) for k, v in env.items()}); os.environ['FQ_TEST_HARNESS'] = '1'
    spec = importlib.util.spec_from_file_location('fq_accel_%d' % random.randrange(1 << 30), SCRIPT); m = importlib.util.module_from_spec(spec)
    old = sys.argv; sys.argv = ['x', R]
    try: spec.loader.exec_module(m)
    finally: sys.argv = old
    return m
def run(c, env=None, plan=None, N=NREADS, fakek='ok'):
    env = dict(env or {}); env.setdefault('FQ_BACKOFF', '0.01'); os.environ['FAKEK'] = fakek
    m = load(c.R, env); M_CUR[0] = m; m.N = N; m.refseq_and_mapping = lambda: None; m.kallisto = lambda: H.FAKE_ID; m.subprocess = H.SubprocShim(c)
    srv = Server(c.src, plan); real = m.urllib.request.urlopen; m.urllib.request.urlopen = srv.urlopen
    try: m.main()
    finally: m.urllib.request.urlopen = real
    return srv, m
def rows(c): return H.manifest_rows(c)
def last(c):
    out = {}
    for d in rows(c): out[(d['study'], d['subset'], d['run'])] = d
    return out
def fq_content(c, acc, sub='A'):
    p = J(c.R, 'data/raw/reads_subsample', f'{acc}_R1_subset{sub}.fastq.gz')
    if not os.path.exists(p): return None
    with gzip.open(p, 'rb') as f: return f.read()
def expect(recs, acc, sub='A', N=NREADS): return b''.join(recs[acc][(0 if sub == 'A' else N):(N if sub == 'A' else 2 * N)])
def logtxt(c):
    with open(J(c.R, 'logs/rawreads.log'), encoding='utf-8') as f: return f.read()

def run_p4(base, runs, RS, wait_s=60, throttle=0.25):
    """P4 hard kill of the whole process tree WHILE >= 3 downloads are verifiably in flight (partial files growing, none validated).
    Returns [(name, ok, detail)]. If the in-flight state is never reached, P4a fails (the test never passes by default)."""
    out = []
    c, recs = setup(base, runs); names = [a for a, _ in runs]; FQ = J(c.R, 'data/raw/reads_subsample')
    code = (f"import sys,os; sys.path.insert(0,{os.path.dirname(os.path.abspath(__file__))!r}); os.environ['FQ_SCRIPT']={SCRIPT!r}; import accel_tests as A, harness as H; "
            f"c=H.Ctx({c.base!r}); A.run(c,{{'FQ_MAX_RESUMES':3,'FQ_WORKERS':3}},{{a:[{{'sleep':{throttle!r}}}] for a in {names!r}}})")   # ~0.25 s per 8 KB read -> ~9 s per download
    kw = {'creationflags': subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == 'nt' else {'start_new_session': True}
    p = subprocess.Popen([sys.executable, '-c', code], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **kw); t0 = time.time(); inflight = {}
    while time.time() - t0 < wait_s and p.poll() is None:
        try: st = last(c) if os.path.exists(J(c.R, 'results/rawreads/run_manifest.tsv')) else {}
        except (KeyError, IndexError, ValueError): st = {}          # a line being written at this instant
        started = [k for k, d in st.items() if d['state'] == 'DOWNLOAD_STARTED']
        parts = {k: J(FQ, f'{k[2]}_R1_subset{k[1]}.fastq.gz.part') for k in started}
        live = {k: pth for k, pth in parts.items() if os.path.exists(pth) and os.path.getsize(pth) > 0}
        if len(live) >= 3: inflight = live; break
        time.sleep(0.1)
    if inflight:
        sizes0 = {k: os.path.getsize(pth) for k, pth in inflight.items()}
    finished_early = p.poll() is not None                 # the run ended before any kill: the precondition cannot be met
    if not finished_early:
        try: H.kill_tree(p)
        except ProcessLookupError: finished_early = True    # exited between the check and the kill
    time.sleep(1.5); alive = p.poll() is None
    before = last(c) if os.path.exists(J(c.R, 'results/rawreads/run_manifest.tsv')) else {}
    still = [k for k in inflight if before.get(k, {}).get('state') == 'DOWNLOAD_STARTED']
    out.append(('P4a precondition: the kill happened while >= 3 parallel downloads were in flight (partial files present and growing, none validated)',
                len(still) >= 3 and not finished_early, f'{len(still)} downloads in flight at the kill ({", ".join(k[2] + "/" + k[1] for k in still)}); waited {time.time()-t0:.1f} s' + ('; the process had already finished - nothing to kill' if finished_early else '')))
    snap = {k: (os.path.getsize(inflight[k]), hashlib.sha256(open(inflight[k], 'rb').read()).hexdigest()) for k in still if os.path.exists(inflight[k])}
    time.sleep(1.0); frozen = all(os.path.getsize(inflight[k]) == snap[k][0] for k in snap)
    falsec = [k for k, d in before.items() if d['state'] == 'COMPLETE' and not os.path.exists(J(c.R, 'results/rawreads/kallisto', *k, '_COMPLETE.json'))]
    grew = all(snap[k][0] >= sizes0[k] for k in snap) if inflight else False
    out.append(('P4b the hard kill reached the process tree: process exited, partial files stopped changing, no false COMPLETE, no FASTQ accepted for interrupted jobs',
                bool(snap) and not alive and frozen and grew and not falsec and not any(os.path.exists(J(FQ, f'{k[2]}_R1_subset{k[1]}.fastq.gz')) for k in still),
                f'process exited={not alive}; partial sizes frozen={frozen}; false COMPLETE={len(falsec)}'))
    run(c, {**RS, 'FQ_WORKERS': 3}); L = last(c); R_ = rows(c)
    kept = {k: J(FQ, f'{k[2]}_R1_subset{k[1]}.fastq.gz.part.prev') for k in snap}
    kept_ok = all(os.path.exists(pth) and hashlib.sha256(open(pth, 'rb').read()).hexdigest() == snap[k][1] for k, pth in kept.items())
    redl = all(sum(1 for d in R_ if (d['study'], d['subset'], d['run']) == k and d['state'] == 'DOWNLOAD_STARTED') == 2 for k in snap)
    ident = all(fq_content(c, a, s_) == expect(recs, a, s_) for (_, s_, a) in L)
    out.append(('P4c restart: every job COMPLETE with reads identical to the source (subsets A and B); each interrupted partial preserved byte-identical as .part.prev; each interrupted job downloaded again exactly once',
                bool(snap) and all(L[k]['state'] == 'COMPLETE' for k in L) and len(L) == 2 * len(runs) and ident and kept_ok and redl,
                f'{len(L)} jobs; identical reads={ident}; {len(kept)} partials preserved byte-identical={kept_ok}; re-downloaded once={redl}'))
    return out
def run_all(base, proj=None, report=print):
    res = []
    def check(name, ok, detail=''): res.append((name, bool(ok), detail)); report(('PASS ' if ok else 'FAIL ') + name + (f'   [{detail}]' if detail else ''))
    RS = {'FQ_MAX_RESUMES': 3, 'FQ_READ_TIMEOUT': 60}
    G = [('GRUN1', 4000)]; size = lambda c, a: os.path.getsize(J(c.src, a + '.fq.gz'))
    # R0/R1 resume after reset
    c, recs = setup(J(base, 'r1'), G); s = size(c, 'GRUN1')
    run(c, RS, {'GRUN1': [{'fault': ('reset', s // 3)}]}); L = last(c)
    fv = [d for d in rows(c) if d['state'] == 'FASTQ_VALIDATED']
    check('R1 connection reset at 1/3 of the stream -> resumed at that byte, COMPLETE', L[('GSE168740', 'A', 'GRUN1')]['state'] == 'COMPLETE' and fv and fv[0]['note'] == 'resumes=1', fv[0]['note'] if fv else '')
    check('R1 written reads byte-identical to reads 1..N of the source (resume adds/loses/duplicates nothing)', fq_content(c, 'GRUN1') == expect(recs, 'GRUN1'))
    # R2 three different faults in one download
    c, recs = setup(J(base, 'r2'), G); s = size(c, 'GRUN1')
    run(c, RS, {'GRUN1': [{'fault': ('reset', s // 5)}, {'fault': ('close', s // 6)}, {'fault': ('timeout', s // 7)}]})
    check('R2 reset + silent early close + read timeout (3 resumes) -> COMPLETE, reads identical', last(c)[('GSE168740', 'A', 'GRUN1')]['state'] == 'COMPLETE' and fq_content(c, 'GRUN1') == expect(recs, 'GRUN1'))
    # R3 more faults than allowed -> fail fast; partial kept; no COMPLETE
    c, recs = setup(J(base, 'r3'), G); s = size(c, 'GRUN1')
    run(c, RS, {'GRUN1': [{'fault': ('reset', s // 20)}] * 4}); L = last(c)
    check('R3 4 faults with FQ_MAX_RESUMES=3 -> DOWNLOAD_FAILED (resumes=3), .part kept, no _COMPLETE.json',
          L[('GSE168740', 'A', 'GRUN1')]['state'] == 'DOWNLOAD_FAILED' and 'resumes=3' in L[('GSE168740', 'A', 'GRUN1')]['note']
          and os.path.exists(J(c.R, 'data/raw/reads_subsample/GRUN1_R1_subsetA.fastq.gz.part')) and not os.path.exists(J(c.R, 'results/rawreads/kallisto/GSE168740/A/GRUN1/_COMPLETE.json')))
    # R3b retry pass in the same run
    c, recs = setup(J(base, 'r3b'), G); s = size(c, 'GRUN1')
    run(c, {**RS, 'FQ_RETRY_PASSES': 1}, {'GRUN1': [{'fault': ('reset', s // 20)}] * 4}); st = [d['state'] for d in rows(c) if d['run'] == 'GRUN1']
    check('R3b retry pass: failed download retried at the end of the same run -> COMPLETE; earlier partial kept as .part.prev',
          st[-1] == 'COMPLETE' and 'DOWNLOAD_FAILED' in st and 'retry pass 1/1: 1 jobs' in logtxt(c) and os.path.exists(J(c.R, 'data/raw/reads_subsample/GRUN1_R1_subsetA.fastq.gz.part.prev'))
          and fq_content(c, 'GRUN1') == expect(recs, 'GRUN1'), ' -> '.join(st))
    # R4-R7 the server does not continue correctly -> refused, never spliced
    for tag, spec, want in [('R4 server ignores Range (HTTP 200 from byte 0)', {'ignore_range': True}, 'did not continue'),
                            ('R5 server claims the right offset but sends shifted bytes', {'lie_shift': 777}, 'overlapping bytes'),
                            ('R6 ETag changed between requests', {'etag': '"v2"'}, 'ETag changed'),
                            ('R7 file length changed between requests', {'length': 10}, 'length changed')]:
        c, recs = setup(J(base, 'r' + tag[1]), G); s = size(c, 'GRUN1')
        run(c, RS, {'GRUN1': [{'fault': ('reset', s // 3)}, spec]}); d = last(c)[('GSE168740', 'A', 'GRUN1')]
        check(f'{tag} -> DOWNLOAD_FAILED, no FASTQ accepted, no _COMPLETE.json',
              d['state'] == 'DOWNLOAD_FAILED' and (want in d['note']) and fq_content(c, 'GRUN1') is None and not os.path.exists(J(c.R, 'results/rawreads/kallisto/GSE168740/A/GRUN1')), d['note'][:90])
    # R8 subset B with a resume
    c, recs = setup(J(base, 'r8'), [('GRUN1', 4000)]); s = size(c, 'GRUN1')
    plan = {'GRUN1': [{}, {'fault': ('reset', s // 2)}]}             # connection 0 = subset A (clean), 1 = subset B (reset), 2 = B resume
    run(c, RS, plan)
    check('R8 subset B (reads N+1..2N) with a reset mid-way -> resumed, reads identical to source records N+1..2N',
          last(c)[('GSE168740', 'B', 'GRUN1')]['state'] == 'COMPLETE' and fq_content(c, 'GRUN1', 'B') == expect(recs, 'GRUN1', 'B'))
    # R9 job time limit -> fail fast and continue with the next job
    c, recs = setup(J(base, 'r9'), [('GRUN1', 4000), ('GRUN2', 2000)]); t0 = time.time()
    run(c, {**RS, 'FQ_JOB_MAX_SECONDS': 2}, {'GRUN1': [{'sleep': 0.3}] * 2}); L = last(c)
    check('R9 slow download exceeds FQ_JOB_MAX_SECONDS -> DOWNLOAD_FAILED (time limit), next job still COMPLETE',
          L[('GSE168740', 'A', 'GRUN1')]['state'] == 'DOWNLOAD_FAILED' and 'time limit' in L[('GSE168740', 'A', 'GRUN1')]['note'] and L[('GSE168740', 'A', 'GRUN2')]['state'] == 'COMPLETE',
          f"{L[('GSE168740','A','GRUN1')]['note'][:60]}; {time.time()-t0:.0f} s")
    # R10 malformed record in the source (sequence/quality length mismatch) with the checks of the resumable path
    c, recs = setup(J(base, 'r10'), G); src = J(c.src, 'GRUN1.fq.gz')
    bad = list(recs['GRUN1']); x = bad[99].split(b'\n'); x[3] = x[3][:-5]; bad[99] = b'\n'.join(x)
    with open(src, 'wb') as f: f.write(gzip.compress(b''.join(bad)))
    run(c, RS); d = last(c)[('GSE168740', 'A', 'GRUN1')]
    check('R10 malformed FASTQ record (quality shorter than sequence) -> DOWNLOAD_FAILED, not accepted', d['state'] == 'DOWNLOAD_FAILED' and 'malformed FASTQ record 100' in d['note'], d['note'][:60])
    # P1-P2 three workers: identical results, one COMPLETE per job, kallisto never concurrent, manifest well-formed
    runs = [(f'GR{i:02d}', 2500) for i in range(6)] + [(f'PR{i:02d}', 1800) for i in range(6)]
    cs, recs = setup(J(base, 'p_seq'), runs); run(cs, {})
    c, recs = setup(J(base, 'p_par'), runs); klog = J(c.base, 'kallisto_times.log'); os.environ['FAKEK_LOG'] = klog; os.environ['FAKEK_SLEEP'] = '0.3'
    plan = {a: [{'fault': ('reset', 20000 + 3000 * i)}] for i, (a, _) in enumerate(runs) if i % 2 == 0}
    try: run(c, {**RS, 'FQ_WORKERS': 3}, plan)
    finally: os.environ.pop('FAKEK_LOG', None); os.environ.pop('FAKEK_SLEEP', None)
    L = last(c); allc = all(L[k]['state'] == 'COMPLETE' for k in L); R_ = rows(c)
    ident = all(fq_content(c, a, sub) == fq_content(cs, a, sub) for (st, sub, a) in L)
    same_ab = all(open(J(c.R, 'results/rawreads/kallisto', *k, 'abundance.tsv')).read() == open(J(cs.R, 'results/rawreads/kallisto', *k, 'abundance.tsv')).read() for k in L)
    check('P1 3 parallel downloads (half with resets): every job COMPLETE; FASTQs and abundance identical to the sequential run', allc and ident and same_ab, f'{len(L)} jobs')
    check('P1 exactly one DOWNLOAD_STARTED per job without a failure and one COMPLETE per job (no job run twice)',
          all(sum(1 for d in R_ if (d['study'], d['subset'], d['run']) == k and d['state'] == 'COMPLETE') == 1 for k in L)
          and all(sum(1 for d in R_ if (d['study'], d['subset'], d['run']) == k and d['state'] == 'DOWNLOAD_STARTED') == 1 for k in L))
    with open(J(c.R, 'results/rawreads/run_manifest.tsv'), encoding='utf-8') as f: lines = f.read().splitlines()
    check('P1 manifest well-formed under concurrency (every row has 15 columns)', all(len(l.split('\t')) == 15 for l in lines), f'{len(lines)-1} rows')
    with open(klog) as f: ev = [l.split() for l in f]
    starts = sorted(float(t) for k, t, _ in ev if k == 'start'); ends = sorted(float(t) for k, t, _ in ev if k == 'end')
    overlap = any(starts[i+1] < ends[i] for i in range(len(starts) - 1))
    check('P2 kallisto never ran concurrently (one quantification at a time)', not overlap and len(starts) == len(L), f'{len(starts)} runs')
    lg = logtxt(c); par = 'acceleration config:' in lg and 'workers=3' in lg
    check('P2 acceleration settings recorded in the log', par)
    # P3 HTTP 403 storm -> throttled to one download, jobs still complete
    c, recs = setup(J(base, 'p403'), runs[:6])
    plan = {a: [{'http403': True}] for a, _ in runs[:6]}
    run(c, {**RS, 'FQ_WORKERS': 3, 'FQ_SLOW_403': 1}, plan); L = last(c)
    check('P3 repeated HTTP 403 -> THROTTLE to 1 download logged; all jobs still COMPLETE', 'THROTTLE' in logtxt(c) and all(L[k]['state'] == 'COMPLETE' for k in L))
    for name, ok, detail in run_p4(J(base, 'p_kill'), runs[:6], RS): check(name, ok, detail)
    # S1-S4 run selection on the REAL project tables (read-only copies) + the real 48-run selection file
    if proj:
        import exclusion_tests as ET
        sel = os.environ.get('ACCEL_SELECTION') or J(proj, 'data', 'metadata', 'run_selection_vg48_2026-11-12.tsv')   # override only for pre-deployment testing
        def rq(tag, env):
            base2 = J(base, 's_' + tag); R = ET._scratch_root(base2, proj); q = []
            m = load(R, env); m.refseq_and_mapping = lambda: None; m.kallisto = lambda: H.FAKE_ID; m.subprocess = H.SubprocShim(H.Ctx(base2))
            m.get = lambda *a, **k: (_ for _ in ()).throw(RuntimeError('no network'))
            m.process_run = lambda K, idx, tmpd, summ, label, sub, r: q.append((label, sub, r['run_accession']))
            m.main(); return q, R, m
        with open(sel, encoding='utf-8') as f: SR = [d for d in csv.DictReader((l for l in f if not l.startswith('#')), delimiter='\t')]
        q, R, m = rq('sel', {'FQ_SELECTION': sel}); auth = [(d['study'], d['subset'], d['run']) for d in sorted(SR, key=lambda d: int(d['priority'])) if d['authorized'] == 'YES']
        with open(J(R, 'logs', 'rawreads.log'), encoding='utf-8') as f: lg = f.read()
        check('S1 real selection file: 48 rows; queue = exactly the 6 authorized jobs (2 retries, 4 probes) in priority order', len(SR) == 48 and q == auth and len(q) == 6, ' '.join(r for _, _, r in q))
        check('S1 full approved job lists still logged with the approved digests (GSE168740 0f7880dc, GSE178742 b823a47e, PRJNA1017469 632d1fae)',
              all(x in lg for x in ('job-list sha256 GSE168740: 0f7880dc9d973f25', 'job-list sha256 GSE178742: b823a47e989a7d62', 'job-list sha256 PRJNA1017469: 632d1faea7743813')) and 'authorized queue sha256' in lg)
        tmp = J(base, 'sel_all_yes.tsv')
        with open(sel, encoding='utf-8') as f, open(tmp, 'w', encoding='utf-8') as g: g.write(f.read().replace('\tmain\tNO\t', '\tmain\tYES\t'))
        q2, _, _ = rq('sel_all', {'FQ_SELECTION': tmp})
        check('S2 after authorizing stage "main": queue = all 48 in priority order, none outside the selection', len(q2) == 48 and q2[:6] == auth and {x[0] for x in q2} <= {'GSE168740', 'PRJNA1017469'})
        hdr = 'priority\tstudy\tsubset\trun\tstage\tauthorized\n'
        for why, txt in [('excluded duplicate listing', hdr + '1\tPRJNA1017469\tA\tSRR26062734\tx\tYES\n'), ('unknown run', hdr + '1\tGSE168740\tA\tSRR99999999\tx\tYES\n'),
                         ('wrong subset', hdr + '1\tPRJNA1017469\tB\tSRR26061248\tx\tYES\n'), ('listed twice', hdr + '1\tGSE168740\tA\tSRR13938902\tx\tYES\n2\tGSE168740\tA\tSRR13938902\tx\tYES\n'),
                         ('authorized not YES/NO', hdr + '1\tGSE168740\tA\tSRR13938902\tx\tY\n')]:
            bp = J(base, 'bad_sel.tsv'); open(bp, 'w', encoding='utf-8').write(txt)
            try: qb, _, _ = rq('bad', {'FQ_SELECTION': bp}); r_ = f'continued ({len(qb)} jobs)'; qb = qb
            except SystemExit as e: r_ = 'stopped: ' + str(e)[:80]; qb = []
            check(f'S3 invalid selection ({why}) -> stops before any job', r_.startswith('stopped') and not qb, r_)
        # S4 the 48 against the real run log (read-only): composition and no overlap with COMPLETE jobs
        man = {}
        with open(J(proj, 'results', 'rawreads', 'run_manifest.tsv'), encoding='utf-8') as f:
            for d in csv.DictReader(f, delimiter='\t'): man[(d['study'], d['subset'], d['run'])] = d['state']
        excl = {l.split('\t')[1] for l in open(J(proj, 'data', 'metadata', 'run_exclusions.tsv'), encoding='utf-8') if l.startswith('PRJNA')}
        with open(J(proj, 'data', 'metadata', 'ena_PRJNA1017469_runs.tsv'), encoding='utf-8') as f: prj = [d['run_accession'] for d in csv.DictReader(f, delimiter='\t') if d['run_accession'] not in excl]
        keys = [(d['study'], d['subset'], d['run']) for d in SR]
        check('S4 none of the 48 is already COMPLETE in the real run log (no completed result redone)', not any(man.get(k) == 'COMPLETE' for k in keys), f"{sum(1 for k in keys if k in man)} of 48 seen in log")
        check('S4 retries = the 2 failed GSE168740 subset-A jobs (SRR13938902 C1_10d, SRR13938905 T1_10d)',
              sorted(d['run'] for d in SR if d['stage'] == 'retry') == ['SRR13938902', 'SRR13938905'] and all(man.get(('GSE168740', 'A', r)) == 'DOWNLOAD_FAILED' for r in ('SRR13938902', 'SRR13938905')))
        pr = [d for d in SR if d['stage'] == 'probe']
        check('S4 probes = 4 unflagged control libraries, one per route x time (1h/4h CTL-IG, 1h/4h CTL-CONT)',
              len(pr) == 4 and sorted('-'.join(d['sample'].split(' ')[1].split('-')[:-1]) for d in pr) == ['1h-CTL-CONT', '1h-CTL-IG', '4h-CTL-CONT', '4h-CTL-IG'] and all(d['approved_qc_decision'].startswith('none') for d in pr))
        mn = {d['run'] for d in SR if d['stage'] == 'main'}
        check('S4 main = all 42 remaining kept PRJNA1017469 libraries (47 - library_10 - 4 probes), authorized=NO, QC-flagged ones included',
              len(mn) == 42 and mn == set(prj) - {'SRR26061202'} - {d['run'] for d in pr} and all(d['authorized'] == 'NO' for d in SR if d['stage'] == 'main')
              and {'SRR26061245', 'SRR26061222', 'SRR26061209', 'SRR26061228', 'SRR26061240'} <= mn)
        # S5 eligibility columns agree with the approved library QC decisions file used by the R analysis
        with open(J(proj, 'data', 'metadata', 'PRJNA1017469_library_qc_decisions.tsv'), encoding='utf-8') as f:
            QD = {d['run_accession']: d['decision'] for d in csv.DictReader((l for l in f if not l.startswith('#')), delimiter='\t')}
        want = {'EXCLUDE_PRIMARY': 'NOT ELIGIBLE', 'PENDING': 'PENDING'}
        bad = []
        for d in SR:
            if d['study'] != 'PRJNA1017469': continue
            q = QD.get(d['run'])
            ok = (d['primary_analysis_eligibility'].startswith(want[q]) if q in want else
                  (q is None and d['primary_analysis_eligibility'].startswith('ELIGIBLE')))
            if not ok or d['sensitivity_analysis_eligibility'] != 'YES': bad.append(d['run'])
        check('S5 every PRJNA1017469 row: primary eligibility matches the approved QC decisions file (12, 33 NOT ELIGIBLE; 17, 45, 28 PENDING; others ELIGIBLE); all sensitivity YES',
              not bad and sum(1 for d in SR if d['primary_analysis_eligibility'].startswith('NOT')) == 2 and sum(1 for d in SR if d['primary_analysis_eligibility'].startswith('PENDING')) == 3
              and {d['run'] for d in SR if d['primary_analysis_eligibility'].startswith('PENDING')} == {'SRR26061240', 'SRR26061209', 'SRR26061228'}, f'mismatches: {bad}')
        check('S5 the 42 held runs are not authorized and nothing in the selection authorizes them automatically (all main rows NO; only priorities 1-6 YES)',
              sorted(int(d['priority']) for d in SR if d['authorized'] == 'YES') == [1, 2, 3, 4, 5, 6] and sum(1 for d in SR if d['stage'] == 'main' and d['authorized'] == 'NO') == 42)
    for k in ENVK: os.environ.pop(k, None)
    return res

if __name__ == '__main__':
    base = os.path.abspath(sys.argv[1]); os.makedirs(base, exist_ok=False)
    if '--only-p4' in sys.argv:          # isolated re-run of the hard-kill test only
        r = run_p4(J(base, 'p_kill'), [(f'GR{i:02d}', 2500) for i in range(6)], {'FQ_MAX_RESUMES': 3, 'FQ_READ_TIMEOUT': 60})
        for n_, ok, d in r: print(('PASS ' if ok else 'FAIL ') + n_ + (f'   [{d}]' if d else ''))
    else:
        proj = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else None
        r = run_all(base, proj)
    print(f'{sum(ok for _, ok, _ in r)}/{len(r)} passed'); sys.exit(0 if all(ok for _, ok, _ in r) else 1)
