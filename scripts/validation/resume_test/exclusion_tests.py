"""Tests for the run-exclusion change (FINAL PROJECT APPROVAL 2026-11-12). Portable (Windows + Linux).
Every write goes into <base>; the project folder is only READ (ENA tables, metadata, exclusion list, manifest, pilot outputs).

Usage: python exclusion_tests.py <empty_base_folder> <project_dir> <original_script.py>   (exit code 1 on any failure)
  original_script.py = the archived pre-change fetch_and_quant.py (reference for the byte-identity test)
The script under test is FQ_SCRIPT (default: <project>/scripts/rawreads/fetch_and_quant.py)."""
import os, sys, shutil, gzip, json, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harness as H
def J(*p): return os.path.join(*p)
STUDIES = ('GSE168740', 'GSE178742', 'PRJNA1017469')

def _load(script, root):
    os.environ['FQ_TEST_HARNESS'] = '1'
    spec = importlib.util.spec_from_file_location('fq_' + str(abs(hash((script, root)))), script); m = importlib.util.module_from_spec(spec)
    old = sys.argv; sys.argv = ['x', root]
    try: spec.loader.exec_module(m)
    finally: sys.argv = old
    return m

def _scratch_root(base, proj, with_exclusions=True, excl_text=None):
    """Copy of the REAL inputs that decide the job list (read from the project, written to scratch)."""
    R = J(base, 'R')
    for d in ['data/metadata', 'data/raw/BeeBase', 'data/reference', 'logs', 'data/raw/reads_subsample']: os.makedirs(J(R, *d.split('/')), exist_ok=True)
    for f in ['ena_PRJNA713798_runs.tsv', 'ena_PRJNA740269_runs.tsv', 'ena_PRJNA1017469_runs.tsv', 'ena_PRJNA1017469_sample_attributes.tsv', 'GSE178742_sample_metadata_verified.tsv']:
        shutil.copy2(J(proj, 'data', 'metadata', f), J(R, 'data', 'metadata', f))
    if excl_text is not None:
        with open(J(R, 'data', 'metadata', 'run_exclusions.tsv'), 'w', encoding='utf-8') as f: f.write(excl_text)
    elif with_exclusions:
        shutil.copy2(J(proj, 'data', 'metadata', 'run_exclusions.tsv'), J(R, 'data', 'metadata', 'run_exclusions.tsv'))
    for p in ['data/raw/BeeBase/amel_OGSv3.2_cds.fa.gz', 'data/reference/amel_OGSv3.2_cds.kidx']:
        with open(J(R, *p.split('/')), 'w') as f: f.write('x\n')
    with open(J(base, 'fake_kallisto.py'), 'w') as f: f.write(H.FAKE_SRC)
    return R

def record_queue(script, base, proj, **kw):
    """Runs main() of `script` with process_run replaced by a recorder: returns the exact ordered queue it would process."""
    R = _scratch_root(base, proj, **kw); c = H.Ctx(base); m = _load(script, R); q = []
    m.refseq_and_mapping = lambda: None; m.kallisto = lambda: H.FAKE_ID; m.subprocess = H.SubprocShim(c)
    def nonet(url, *a, **k): raise RuntimeError('no network in this test: ' + url)
    m.get = nonet
    m.process_run = lambda K, idx, tmpd, summ, label, sub, r: q.append((label, sub, r['run_accession'], r.get('fastq_ftp', '')))
    m.main(); return q, m, R

def digest(q, study):
    import hashlib; h = hashlib.sha256()
    for label, sub, run, ftp in q:
        if label == study: h.update(f'{label}\t{sub}\t{run}\t{ftp}\n'.encode())
    return h.hexdigest()

def read_excl(proj):
    with open(J(proj, 'data', 'metadata', 'run_exclusions.tsv'), encoding='utf-8') as f: L = [l for l in f.read().splitlines() if l.strip() and not l.startswith('#')]
    h = L[0].split('\t'); return [dict(zip(h, l.split('\t'))) for l in L[1:]]

def manifest(R):
    p = J(R, 'results', 'rawreads', 'run_manifest.tsv')
    if not os.path.exists(p): return []
    with open(p, encoding='utf-8') as f: L = f.read().splitlines()
    h = L[0].split('\t'); return [dict(zip(h, l.split('\t'))) for l in L[1:]]

def run_all(base, proj, original, report=print):
    res = []
    def check(name, cond, detail=''): res.append((name, bool(cond), detail)); report(('PASS ' if cond else 'FAIL ') + name + (f'   [{detail}]' if detail else ''))
    script = os.environ.get('FQ_SCRIPT') or J(proj, 'scripts', 'rawreads', 'fetch_and_quant.py')
    excl = read_excl(proj); ex_runs = {d['run_accession'] for d in excl}; kept_runs = {d['kept_run'] for d in excl}
    # ---- E0 the approved exclusion list itself
    check('E0.1 exclusion list has 47 rows, all PRJNA1017469, 47 distinct excluded and 47 distinct kept runs, disjoint',
          len(excl) == 47 and all(d['study'] == 'PRJNA1017469' for d in excl) and len(ex_runs) == 47 and len(kept_runs) == 47 and not (ex_runs & kept_runs),
          f'{len(excl)} rows, {len(ex_runs)} excluded, {len(kept_runs)} kept')
    check('E0.2 every row carries reason, evidence and approval', all(d.get('reason') and d.get('evidence') and d.get('decided_by') for d in excl))
    # ---- E1 queue on the REAL ENA tables: original script vs changed script
    q0, _, _ = record_queue(original, J(base, 'e1_original'), proj, with_exclusions=False)
    q1, m1, R1 = record_queue(script, J(base, 'e1_new'), proj)
    for st in ('GSE168740', 'GSE178742'):
        a = [x for x in q0 if x[0] == st]; b = [x for x in q1 if x[0] == st]
        check(f'E1 {st} job list byte-identical to the original script (order, subset, run, URL)', a == b and digest(q0, st) == digest(q1, st),
              f'{len(a)} jobs; sha256 {digest(q1, st)[:16]}')
    p0 = [x for x in q0 if x[0] == 'PRJNA1017469']; p1 = [x for x in q1 if x[0] == 'PRJNA1017469']
    check('E1 PRJNA1017469: original queue 94 runs, new queue 47 runs', len(p0) == 94 and len(p1) == 47, f'{len(p0)} -> {len(p1)}')
    check('E1 PRJNA1017469: new queue = original queue minus the excluded runs, same order', p1 == [x for x in p0 if x[2] not in ex_runs])
    check('E1 PRJNA1017469: the 47 queued runs are exactly the 47 kept runs (one per library)', {x[2] for x in p1} == kept_runs)
    check('E1 no excluded run anywhere in the new queue', not any(x[2] in ex_runs for x in q1))
    check('E1 total queue = original total - 47', len(q1) == len(q0) - 47, f'{len(q0)} -> {len(q1)}')
    ms = manifest(R1); ed = [d for d in ms if d['state'] == 'EXCLUDED_DUPLICATE']
    check('E2 first start writes one EXCLUDED_DUPLICATE manifest row per excluded run (47) naming the kept run',
          len(ed) == 47 and {d['run'] for d in ed} == ex_runs and all('kept SRR' in d['note'] and 'approved' in d['note'] for d in ed), f'{len(ed)} rows')
    m2 = _load(script, R1); q2 = []
    m2.refseq_and_mapping = lambda: None; m2.kallisto = lambda: H.FAKE_ID; m2.subprocess = H.SubprocShim(H.Ctx(J(base, 'e1_new')))
    m2.get = lambda *a, **k: (_ for _ in ()).throw(RuntimeError('no network'))
    m2.process_run = lambda K, idx, tmpd, summ, label, sub, r: q2.append((label, sub, r['run_accession'], r.get('fastq_ftp', '')))
    m2.main(); ed2 = [d for d in manifest(R1) if d['state'] == 'EXCLUDED_DUPLICATE']
    check('E2 restart does not duplicate audit rows (still 47) and gives the same queue', len(ed2) == 47 and q2 == q1, f'{len(ed2)} rows')
    with open(J(R1, 'logs', 'rawreads.log'), encoding='utf-8') as f: lg = f.read()
    check('E2 log records the per-study job-list sha256', all(f'job-list sha256 {st}: {digest(q1, st)}' in lg for st in STUDIES))
    # ---- E3 an inconsistent exclusion list stops the script before anything is processed
    hdr = 'study\trun_accession\treason\tkept_run\tevidence\tdecided_by\tdate\n'
    bad = {'kept run itself excluded': hdr + 'PRJNA1017469\tSRR26062734\tx\tSRR26061248\te\tt\td\nPRJNA1017469\tSRR26061248\tx\tSRR26062734\te\tt\td\n',
           'excluded run not in ENA table': hdr + 'PRJNA1017469\tSRR99999999\tx\tSRR26061248\te\tt\td\n',
           'kept run not in ENA table': hdr + 'PRJNA1017469\tSRR26062734\tx\tSRR99999999\te\tt\td\n'}
    for i, (why, txt) in enumerate(bad.items()):
        try: q, _, Rb = record_queue(script, J(base, f'e3_{i}'), proj, excl_text=txt); r = f'continued ({len(q)} jobs)'
        except SystemExit as e: r = 'stopped: ' + str(e)[:90]; q = []
        check(f'E3 inconsistent list ({why}) -> stops, nothing processed', r.startswith('stopped') and not q, r)
    # ---- E4 full mock pipeline: an excluded run is never downloaded or quantified; its kept twin completes
    c = H.setup(J(base, 'e4_mock'))
    H.add_run(c, 'ena_PRJNA1017469_runs.tsv', 'KEEP1', 'RUN1'); H.add_run(c, 'ena_PRJNA1017469_runs.tsv', 'DUP1', 'RUN1')
    with open(J(c.R, 'data', 'metadata', 'run_exclusions.tsv'), 'w', encoding='utf-8') as f:
        f.write('# test\n' + hdr + 'PRJNA1017469\tDUP1\tDUPLICATE_LISTING of library_x\tKEEP1\ttest evidence\ttest approval\t2026-11-12\n')
    os.environ['FQ_SCRIPT'] = script; H.run(c); H.run(c)
    rows = H.manifest_rows(c); st, _ = H.state(c)
    check('E4 kept run COMPLETE; excluded run only ever EXCLUDED_DUPLICATE (once), never DOWNLOAD_STARTED',
          st.get(('PRJNA1017469', 'A', 'KEEP1')) == 'COMPLETE' and [d['state'] for d in rows if d['run'] == 'DUP1'] == ['EXCLUDED_DUPLICATE'])
    check('E4 no FASTQ and no output folder for the excluded run', not any(f.startswith('DUP1') for f in os.listdir(J(c.R, 'data/raw/reads_subsample')))
          and not os.path.exists(J(c.R, 'results/rawreads/kallisto/PRJNA1017469/A/DUP1')))
    check('E4 other studies unaffected in the mock (RUN1-3, SHORT complete)', all(st.get(('GSE168740', 'A', r)) == 'COMPLETE' for r in ('RUN1', 'RUN2', 'RUN3')) and st.get(('GSE178742', 'A', 'SHORT')) == 'COMPLETE')
    # ---- E5 pilot outputs in the REAL project are reused (read-only check: is_complete() == True means process_run returns at once)
    with open(J(proj, 'results', 'rawreads', 'run_manifest.tsv'), encoding='utf-8') as f: L = f.read().splitlines()
    h = L[0].split('\t'); prow = [dict(zip(h, l.split('\t'))) for l in L[1:]]
    lastst = {}
    for d in prow: lastst[(d['study'], d['subset'], d['run'])] = d['state']
    done = sorted(k for k, v in lastst.items() if v == 'COMPLETE')                       # pilot runs + production runs completed so far
    PILOTS = {('GSE168740', 'A', 'SRR13938890'), ('GSE168740', 'A', 'SRR13938893'), ('GSE178742', 'A', 'SRR14887980'), ('GSE178742', 'A', 'SRR14887967'), ('PRJNA1017469', 'A', 'SRR26061202')}
    mp = _load(script, proj); qset = {(x[0], x[1], x[2]) for x in q1}; bad = []
    for label, sub, run in done:
        od = J(proj, 'results', 'rawreads', 'kallisto', label, sub, run); fq = J(proj, 'data', 'raw', 'reads_subsample', f'{run}_R1_subset{sub}.fastq.gz')
        ok, why = mp.is_complete(od, fq)
        if not (ok and (label, sub, run) in qset): bad.append(f'{run}:{why}')
    for label, sub, run in sorted(PILOTS):
        od = J(proj, 'results', 'rawreads', 'kallisto', label, sub, run); fq = J(proj, 'data', 'raw', 'reads_subsample', f'{run}_R1_subset{sub}.fastq.gz')
        ok, why = mp.is_complete(od, fq)
        check(f'E5 pilot {label} {sub} {run}: in the new queue and re-validates as complete -> skipped, not re-downloaded or re-quantified', ok and (label, sub, run) in qset, why)
    check('E5 every job COMPLETE in the run log (pilot + production) re-validates and is in the queue -> skipped, never redone', bool(done) and not bad, f'{len(done)} complete; problems: {bad[:3]}')
    return res

if __name__ == '__main__':
    base, proj, orig = (os.path.abspath(a) for a in sys.argv[1:4]); os.makedirs(base, exist_ok=False)
    r = run_all(base, proj, orig)
    print(f'{sum(ok for _, ok, _ in r)}/{len(r)} passed'); sys.exit(0 if all(ok for _, ok, _ in r) else 1)
