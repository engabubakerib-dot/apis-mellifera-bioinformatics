r"""INDEPENDENT verification of the raw-read pilot (read-only on project data; Python >= 3.8, stdlib only).
Deliberately does NOT import scripts/rawreads/fetch_and_quant.py or pilot_runs.py: every check is re-implemented here
from the files themselves, and the pilot's own report is then compared against these results.

Usage: python verify_pilot_independent.py <project_dir> <output_dir> [--negative-control <scratch_dir>]
  writes <output_dir>/independent_verification_<ts>.tsv and .txt
  --negative-control: copies ONE completed pilot run (FASTQ + output folder + manifest) into <scratch_dir> (must be outside
  the project data folders), alters the COPY of the FASTQ (a) by appending one byte and (b) by flipping one byte (same size),
  and runs the deployed pipeline's process_run() on the copy. Expected: REVALIDATION_FAILED, copy left untouched, real
  project files unchanged (their SHA-256 is compared before/after)."""
import os, sys, re, json, gzip, hashlib, time, shutil, glob
PROJ = os.path.abspath(sys.argv[1]); OUT = os.path.abspath(sys.argv[2]); os.makedirs(OUT, exist_ok=True)
TS = time.strftime('%Y%m%d_%H%M%S'); TARGET = 500000; NTARGETS = 15314
PD = os.path.join(PROJ, 'results', 'rawreads', 'pilot'); MAN = os.path.join(PROJ, 'results', 'rawreads', 'run_manifest.tsv')
lines_out = []
def say(s=''): print(s, flush=True); lines_out.append(s)
def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
def fastq_count(p):
    """strict: gzip must decode to EOF; 4-line records; '@' header and '+' separator; seq/qual same length."""
    n = 0
    with gzip.open(p, 'rb') as f:
        while True:
            h = f.readline()
            if not h: break
            s, plus, q = f.readline(), f.readline(), f.readline()
            if not q or not h.startswith(b'@') or not plus.startswith(b'+') or len(s.rstrip(b'\r\n')) != len(q.rstrip(b'\r\n')):
                raise ValueError(f'malformed record after {n} reads')
            n += 1
    return n
def run_info(p):
    t = open(p, encoding='utf-8', errors='replace').read()
    try: json.loads(t); valid = 'valid JSON'
    except ValueError: valid = 'INVALID JSON (kallisto Windows backslashes)' if re.search(r'"call"\s*:\s*"[^"]*\\[^\\"/bfnrtu]', t) else 'INVALID JSON (other)'
    g = lambda k, pat=r'-?\d+': (lambda m: m.group(1) if m else None)(re.search(r'"%s"\s*:\s*(%s)' % (re.escape(k), pat), t))
    return dict(n_processed=int(g('n_processed')) if g('n_processed') else None, n_pseudoaligned=int(g('n_pseudoaligned')) if g('n_pseudoaligned') else None,
                p_pseudoaligned=float(g('p_pseudoaligned', r'-?\d+(?:\.\d+)?')) if g('p_pseudoaligned', r'-?\d+(?:\.\d+)?') else None,
                n_bootstraps=g('n_bootstraps'), version=(lambda m: m.group(1) if m else None)(re.search(r'"kallisto_version"\s*:\s*"([^"]*)"', t)), json=valid)
def manifest():
    if not os.path.exists(MAN): return []
    L = open(MAN, encoding='utf-8').read().splitlines(); h = L[0].split('\t'); return [dict(zip(h, l.split('\t'))) for l in L[1:]]

say(f'INDEPENDENT PILOT VERIFICATION {TS}  project={PROJ}')
sel_f = os.path.join(PD, 'pilot_selection.json')
if not os.path.exists(sel_f): say('PILOT NOT RUN: results/rawreads/pilot/pilot_selection.json missing'); sys.exit(2)
sel = json.load(open(sel_f, encoding='utf-8'))
snap_f = os.path.join(PD, 'snapshot_after_pass1.json'); snap = json.load(open(snap_f, encoding='utf-8')) if os.path.exists(snap_f) else None
reports = sorted(glob.glob(os.path.join(PD, 'pilot_report_*.txt'))); rep_txt = open(reports[-1], encoding='utf-8').read() if reports else ''
plog = open(os.path.join(PD, 'pilot_log.txt'), encoding='utf-8').read() if os.path.exists(os.path.join(PD, 'pilot_log.txt')) else ''
M = manifest(); pilot_runs = {r['run_accession'] for _, r, _ in sel}
say(f'selection: {len(sel)} runs; latest pilot report: {os.path.basename(reports[-1]) if reports else "NONE"}; snapshot after pass 1: {"present" if snap else "MISSING"}; manifest rows: {len(M)}')
others = sorted({d['run'] for d in M} - pilot_runs)
say(f'manifest rows for NON-pilot runs (must be none: full pipeline not started): {others or "none"}')
pass2_start = None
for l in plog.splitlines():
    if 'PILOT pass 2 (restart' in l: pass2_start = l[:19]
rows = []
for label, r, why in sel:
    run = r['run_accession']; od = os.path.join(PROJ, 'results', 'rawreads', 'kallisto', label, 'A', run)
    fq = os.path.join(PROJ, 'data', 'raw', 'reads_subsample', f'{run}_R1_subsetA.fastq.gz'); v = {'study': label, 'run': run, 'sample': f"{r.get('sample_alias','')} {r.get('sample_title','')}"}
    issues = []
    if not os.path.exists(fq): issues.append('FASTQ missing')
    else:
        try: v['fastq_reads'] = fastq_count(fq)
        except Exception as e: v['fastq_reads'] = f'ERROR {e}'; issues.append('FASTQ unreadable/malformed')
        v['fastq_sha256'] = sha256(fq); v['fastq_bytes'] = os.path.getsize(fq)
    cj = os.path.join(od, '_COMPLETE.json')
    try: rec = json.load(open(cj, encoding='utf-8'))
    except Exception as e: rec = {}; issues.append('_COMPLETE.json missing/unreadable')
    v['target'] = rec.get('reads_target')
    if rec:
        if rec.get('reads_target') != TARGET: issues.append(f"target {rec.get('reads_target')} != {TARGET}")
        if v.get('fastq_reads') != rec.get('reads_written') or rec.get('reads_written') != rec.get('reads_recounted'): issues.append('read count != _COMPLETE.json')
        if v.get('fastq_reads') != TARGET and rec.get('fastq_eof') is not True: issues.append('fewer reads than target without source EOF')
        v['sha_ok'] = v.get('fastq_sha256') == rec.get('fastq_sha256') and v.get('fastq_bytes') == rec.get('fastq_bytes')
        if not v['sha_ok']: issues.append('FASTQ SHA-256/size != _COMPLETE.json')
        v['kallisto_rc'] = rec.get('kallisto_rc')
        if rec.get('kallisto_rc') != 0: issues.append('kallisto exit code != 0')
    ri = os.path.join(od, 'run_info.json')
    if os.path.exists(ri) and os.path.getsize(ri) > 0:
        info = run_info(ri); v.update(n_processed=info['n_processed'], n_pseudoaligned=info['n_pseudoaligned'], p_pseudoaligned=info['p_pseudoaligned'], run_info_json=info['json'], kallisto_version=info['version'], n_bootstraps=info['n_bootstraps'])
        if info['n_processed'] != v.get('fastq_reads'): issues.append('kallisto n_processed != FASTQ reads')
        if info['n_processed'] and info['n_pseudoaligned'] is not None and abs(100 * info['n_pseudoaligned'] / info['n_processed'] - (info['p_pseudoaligned'] or -1)) > 0.1: issues.append('p_pseudoaligned inconsistent')
        if rec and (rec.get('n_processed') != info['n_processed'] or rec.get('p_pseudoaligned') != info['p_pseudoaligned']): issues.append('run_info.json != _COMPLETE.json')
        if info['version'] != '0.51.1': issues.append(f"kallisto version {info['version']}")
    else: issues.append('run_info.json missing/empty')
    ab = os.path.join(od, 'abundance.tsv')
    if os.path.exists(ab):
        L = open(ab, encoding='utf-8').read().splitlines(); v['abundance_targets'] = len(L) - 1
        if L[0].split('\t') != ['target_id', 'length', 'eff_length', 'est_counts', 'tpm']: issues.append('abundance.tsv header unexpected')
        ids = [l.split('\t')[0] for l in L[1:]]; est = [float(l.split('\t')[3]) for l in L[1:]]
        if len(set(ids)) != len(ids) or v['abundance_targets'] != NTARGETS: issues.append(f"abundance targets {v['abundance_targets']} (expected {NTARGETS}, unique)")
        if min(est) < 0: issues.append('negative est_counts')
        v['sum_est_counts'] = round(sum(est), 1)
        if v.get('n_pseudoaligned') and abs(sum(est) - v['n_pseudoaligned']) > 0.01 * v['n_pseudoaligned']: issues.append('sum(est_counts) differs >1% from n_pseudoaligned')
    else: issues.append('abundance.tsv missing')
    h5 = os.path.join(od, 'abundance.h5')
    v['h5_ok'] = os.path.exists(h5) and open(h5, 'rb').read(8) == b'\x89HDF\r\n\x1a\n'
    if not v['h5_ok']: issues.append('abundance.h5 missing or not HDF5')
    se = os.path.join(od, 'kallisto_stderr.txt')
    if os.path.exists(se):
        m = re.search(r'processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned', open(se, encoding='utf-8', errors='replace').read())
        if not m or int(m.group(1).replace(',', '')) != v.get('n_processed'): issues.append('kallisto stderr read count disagrees')
    else: issues.append('kallisto_stderr.txt missing')
    mr = [d for d in M if (d['study'], d['subset'], d['run']) == (label, 'A', run)]
    v['manifest_states'] = '>'.join(d['state'] for d in mr); v['manifest_final'] = mr[-1]['state'] if mr else 'NONE'
    if v['manifest_final'] != 'COMPLETE': issues.append('manifest final state not COMPLETE')
    fv = [d for d in mr if d['state'] == 'FASTQ_VALIDATED']
    if not fv or fv[-1]['fastq_sha256'] != v.get('fastq_sha256'): issues.append('manifest FASTQ_VALIDATED sha256 != FASTQ')
    v['marker_valid'] = bool(rec) and not any(x for x in issues if 'COMPLETE.json' in x or 'SHA-256' in x)
    # restart: files identical to the pass-1 snapshot; no manifest row for this run after pass 2 started; pass-2 report says PASS
    if snap is None or pass2_start is None: v['restart'] = 'NOT TESTED (pass 2 not run)'; issues.append('restart not tested')
    else:
        keys = [k for k in snap if run in k]
        P = lambda k: os.path.join(PROJ, *re.split(r'[\\/]', k))   # snapshot keys may carry Windows separators
        same = bool(keys) and all(os.path.exists(P(k)) and sha256(P(k)) == snap[k][2] and os.path.getsize(P(k)) == snap[k][1] for k in keys)
        mrc = re.search(r'manifest rows before restart (\d+), after restart (\d+)', rep_txt)
        later = [d['state'] for d in M[int(mrc.group(1)):] if (d['study'], d['subset'], d['run']) == (label, 'A', run)] if mrc else ['UNKNOWN (row count not in report)']
        rep_ok = bool(re.search(re.escape(run) + r'.*?restart: PASS', rep_txt, re.S)) and ('restart: FAIL' not in rep_txt.split(run, 1)[1].split('\n\n', 1)[0] if run in rep_txt else False)
        v['restart'] = 'SKIPPED CORRECTLY' if same and not later and rep_ok else f'PROBLEM: files_identical={same} rows_after_pass2={later} report_says_pass={rep_ok}'
        if not v['restart'].startswith('SKIPPED'): issues.append('restart check failed')
    rep_line = next((l for l in rep_txt.splitlines() if l.startswith(f'{label} {run} ')), '')
    v['pilot_report_verdict'] = rep_line.split(': ', 1)[1] if ': ' in rep_line else 'NOT IN REPORT'
    v['independent_verdict'] = 'PASS' if not issues else 'FAIL: ' + '; '.join(issues)
    if v['pilot_report_verdict'].startswith('PASS') != (not issues): v['independent_verdict'] += '  [DISAGREES WITH PILOT REPORT]'
    rows.append(v)
prev = os.path.join(PROJ, 'data', 'raw', 'reads_subsample', 'SRR13938890_R1_subsetA.fastq.gz.prev')
say(f"legacy partial SRR13938890 file preserved as .prev: {os.path.exists(prev)}{' (' + str(os.path.getsize(prev)) + ' bytes; original 20,316,182)' if os.path.exists(prev) else ''}")
cols = ['study', 'run', 'sample', 'fastq_reads', 'target', 'sha_ok', 'fastq_sha256', 'kallisto_rc', 'n_processed', 'p_pseudoaligned', 'abundance_targets', 'sum_est_counts', 'n_pseudoaligned',
        'h5_ok', 'run_info_json', 'kallisto_version', 'n_bootstraps', 'marker_valid', 'manifest_final', 'manifest_states', 'restart', 'pilot_report_verdict', 'independent_verdict']
with open(os.path.join(OUT, f'independent_verification_{TS}.tsv'), 'w', encoding='utf-8') as f:
    f.write('\t'.join(cols) + '\n'); [f.write('\t'.join(str(v.get(c, '')) for c in cols) + '\n') for v in rows]
for v in rows:
    say(''); say(f"{v['study']} {v['run']} ({v['sample']}): {v['independent_verdict']}")
    for c in cols[3:-1]: say(f'  {c}: {v.get(c, "")}')

# ------------------------------------------------------------------ negative control on a COPY
if '--negative-control' in sys.argv:
    S = os.path.abspath(sys.argv[sys.argv.index('--negative-control') + 1])
    if S.startswith(os.path.join(PROJ, 'data')) or S.startswith(os.path.join(PROJ, 'results')): raise SystemExit('scratch must not be inside data/ or results/')
    done = [v for v in rows if v['independent_verdict'].startswith('PASS')]
    say(''); say('NEGATIVE CONTROL (altered FASTQ, on a COPY outside the project data)')
    if not done: say('  not run: no fully verified pilot run to copy')
    else:
        v = done[0]; label, run = v['study'], v['run']
        real = [os.path.join(PROJ, 'data', 'raw', 'reads_subsample', f'{run}_R1_subsetA.fastq.gz')] + glob.glob(os.path.join(PROJ, 'results', 'rawreads', 'kallisto', label, 'A', run, '*'))
        real_before = {p: sha256(p) for p in real}
        import importlib.util
        for variant in ('append_one_byte', 'flip_one_byte_same_size'):
            R = os.path.join(S, f'negctl_{TS}_{variant}'); os.makedirs(os.path.join(R, 'data', 'raw', 'reads_subsample')); os.makedirs(os.path.join(R, 'results', 'rawreads', 'kallisto', label, 'A'))
            fq = os.path.join(R, 'data', 'raw', 'reads_subsample', f'{run}_R1_subsetA.fastq.gz'); shutil.copy2(real[0], fq)
            od = os.path.join(R, 'results', 'rawreads', 'kallisto', label, 'A', run); shutil.copytree(os.path.join(PROJ, 'results', 'rawreads', 'kallisto', label, 'A', run), od)
            shutil.copy2(MAN, os.path.join(R, 'results', 'rawreads', 'run_manifest.tsv'))
            if variant == 'append_one_byte':
                with open(fq, 'ab') as f: f.write(b'x')
            else:
                with open(fq, 'r+b') as f: f.seek(os.path.getsize(fq) // 2); b = f.read(1); f.seek(-1, 1); f.write(bytes([b[0] ^ 0xFF]))
            out_before = {p: sha256(os.path.join(od, p)) for p in os.listdir(od)}; fq_before = sha256(fq)
            sys.argv = ['fetch_and_quant.py', R]
            spec = importlib.util.spec_from_file_location('fq_neg', os.path.join(PROJ, 'scripts', 'rawreads', 'fetch_and_quant.py')); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
            rr = next(x for _, x, _ in sel if x['run_accession'] == run)
            m.process_run(None, None, os.path.join(R, 'data', 'raw', 'reads_subsample'), os.path.join(R, 'subsample_log.tsv'), label, 'A', rr)   # no kallisto/network needed if it correctly refuses
            L = open(os.path.join(R, 'results', 'rawreads', 'run_manifest.tsv'), encoding='utf-8').read().splitlines()
            st = L[-1].split('\t')[4]; note = L[-1].split('\t')[-1]
            untouched = out_before == {p: sha256(os.path.join(od, p)) for p in os.listdir(od)} and fq_before == sha256(fq) and not os.path.exists(fq + '.prev')
            say(f"  {variant}: pipeline state {st} ({note}); copy untouched: {untouched} -> {'PASS (fails as expected)' if st == 'REVALIDATION_FAILED' and untouched else 'FAIL'}")
        say(f"  real project files unchanged by the negative control: {real_before == {p: sha256(p) for p in real}}")
allp = all(v['independent_verdict'] == 'PASS' for v in rows) and len(rows) == 5
say(''); say('INDEPENDENT RESULT: ' + ('ALL 5 RUNS PASS every technical check' if allp else 'NOT PASSED - see entries above'))
open(os.path.join(OUT, f'independent_verification_{TS}.txt'), 'w', encoding='utf-8').write('\n'.join(lines_out) + '\n')
