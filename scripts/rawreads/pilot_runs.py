r"""Controlled raw-read PILOT (5 runs, subset A). Uses the DEPLOYED, preflight-approved scripts/rawreads/fetch_and_quant.py
UNCHANGED: it is imported as a module, its own preflight gate is called, and every run goes through its own process_run()
(same download, read-count, SHA-256, kallisto and _COMPLETE.json logic as the full pipeline). Outputs land in the normal
results/rawreads/... locations, so the full run later skips these runs as completed. Nothing is deleted.

Usage (from the project folder, launcher's interpreter):
  python scripts\rawreads\pilot_runs.py <project> run       pass 1: process the 5 pilot runs, verify, snapshot
  python scripts\rawreads\pilot_runs.py <project> restart   pass 2 (NEW process = real restart): re-run the same 5 jobs;
                                                            they must be skipped after SHA-256 re-verification, with no file
                                                            changed and no manifest row added; writes the final report
Report: results/rawreads/pilot/pilot_report_<timestamp>.tsv/.txt (new files every restart pass; plus pass snapshots and pilot_log.txt)."""
import os, sys, json, time, hashlib, importlib.util, gzip
PROJ = os.path.abspath(sys.argv[1]); MODE = sys.argv[2]
FQS = os.path.join(PROJ, 'scripts', 'rawreads', 'fetch_and_quant.py')
sys.argv = [FQS, PROJ]
spec = importlib.util.spec_from_file_location('fetch_and_quant', FQS); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
PD = os.path.join(PROJ, 'results', 'rawreads', 'pilot'); os.makedirs(PD, exist_ok=True)
# Pilot selection (fixed in advance; subset A only; the experimental design is not changed)
SELECT = [('GSE168740', 'PRJNA713798', 'sample_alias', 'GSM5163793', 'C1_1d control; also tests recovery of the earlier interrupted partial FASTQ'),
          ('GSE168740', 'PRJNA713798', 'sample_alias', 'GSM5163796', 'T1_1d dinotefuran-exposed'),
          ('GSE178742', 'PRJNA740269', 'sample_alias', 'GSM5396748', 'C21-F1 apiary forager, exposed colony'),
          ('GSE178742', 'PRJNA740269', 'sample_alias', 'GSM5396735', '1J7-F field forager, unexposed site'),
          ('PRJNA1017469', 'PRJNA1017469', 'first', '', 'first run by accession (design still to be established)')]
def plog(msg):
    s = time.strftime('%Y-%m-%dT%H:%M:%S') + '  ' + msg; print(s, flush=True)
    with open(os.path.join(PD, 'pilot_log.txt'), 'a', encoding='utf-8') as f: f.write(s + '\n')
def paths(label, run):
    od = os.path.join(PROJ, 'results', 'rawreads', 'kallisto', label, 'A', run)
    fq = os.path.join(PROJ, 'data', 'raw', 'reads_subsample', f'{run}_R1_subsetA.fastq.gz')
    return od, fq
def manifest_rows():
    p = m.manifest_path()
    if not os.path.exists(p): return []
    with open(p, encoding='utf-8') as f: L = f.read().splitlines()
    h = L[0].split('\t'); return [dict(zip(h, l.split('\t'))) for l in L[1:]]
def snapshot(jobs):
    s = {}
    for label, r, _ in jobs:
        od, fq = paths(label, r['run_accession'])
        for p in ([fq] + ([os.path.join(od, f) for f in sorted(os.listdir(od))] if os.path.isdir(od) else [])):
            if os.path.isfile(p): s[os.path.relpath(p, PROJ)] = [os.path.getmtime(p), os.path.getsize(p), m.sha256_file(p)]
    return s
def setup():
    m.preflight_gate()                                   # identical gate to the full pipeline (sha256 of fetch_and_quant.py)
    K = m.kallisto(); idx = os.path.join(PROJ, 'data', 'reference', 'amel_OGSv3.2_cds.kidx')
    if not os.path.exists(idx): raise SystemExit('kallisto index missing: ' + idx)
    tmpd = os.path.join(PROJ, 'data', 'raw', 'reads_subsample'); os.makedirs(tmpd, exist_ok=True)
    summ = os.path.join(PROJ, 'results', 'rawreads', 'subsample_log.tsv')
    rows = {p: m.ena_runs(p) for p in ('PRJNA713798', 'PRJNA740269', 'PRJNA1017469')}
    try: m.biosample_attrs(rows['PRJNA1017469'], 'PRJNA1017469')
    except Exception as e: plog('biosample attributes for PRJNA1017469 failed (not needed for the pilot): ' + str(e))
    jobs = []
    for label, proj, how, key, why in SELECT:
        rr = sorted(rows[proj], key=lambda r: r['run_accession'])
        r = rr[0] if how == 'first' else next((x for x in rr if x.get('sample_alias') == key), None)
        if r is None: plog(f'SELECTION FAILED: {label} {key} not found in ena_{proj}_runs.tsv'); continue
        jobs.append((label, r, why))
    return K, idx, tmpd, summ, jobs
def run_job(K, idx, tmpd, summ, label, r):
    try: m.process_run(K, idx, tmpd, summ, label, 'A', r)
    except Exception as e:                               # same handling as the full pipeline's main loop
        m.manifest_add(label, 'A', r['run_accession'], 'ERROR', note=str(e)[:300]); plog(f"  {label} A {r['run_accession']} ERROR {e}")
def verify(label, r):
    run = r['run_accession']; od, fq = paths(label, run); v = {'study': label, 'run': run, 'sample_alias': r.get('sample_alias', ''), 'title': r.get('sample_title', '')}
    states = [d['state'] for d in manifest_rows() if (d['study'], d['subset'], d['run']) == (label, 'A', run)]
    v['manifest_states'] = '>'.join(states); v['manifest_last'] = states[-1] if states else ''
    cj = os.path.join(od, '_COMPLETE.json'); v['complete_marker'] = os.path.exists(cj)
    if not v['complete_marker']: v['verdict'] = 'NOT COMPLETE'; return v
    with open(cj, encoding='utf-8') as f: rec = json.load(f)
    v['reads_target'] = rec['reads_target']; v['reads_written'] = rec['reads_written']; v['fastq_eof'] = rec['fastq_eof']
    try: v['reads_recounted_now'] = m.count_reads_strict(fq)
    except Exception as e: v['reads_recounted_now'] = 'ERROR ' + type(e).__name__ + ': ' + str(e)[:80]
    v['sha256_recorded'] = rec['fastq_sha256']
    v['sha256_now'] = m.sha256_file(fq) if os.path.exists(fq) else 'FASTQ MISSING'
    v['bytes_ok'] = os.path.exists(fq) and os.path.getsize(fq) == rec['fastq_bytes']
    v['kallisto_rc'] = rec['kallisto_rc']; v['n_processed'] = rec['n_processed']; v['p_pseudoaligned'] = rec['p_pseudoaligned']; v['run_info_parsed_with'] = rec.get('run_info_parsed_with', '')
    files = {f: os.path.getsize(os.path.join(od, f)) if os.path.exists(os.path.join(od, f)) else 0 for f in ('abundance.tsv', 'abundance.h5', 'run_info.json')}
    try:
        with open(os.path.join(od, 'abundance.tsv'), encoding='utf-8') as f: v['abundance_rows'] = sum(1 for _ in f) - 1
    except Exception as e: v['abundance_rows'] = 'ERROR ' + str(e)[:80]
    mv = [d for d in manifest_rows() if (d['study'], d['subset'], d['run'], d['state']) == (label, 'A', run, 'FASTQ_VALIDATED')]
    checks = {'reads == target (or source ended, intact)': rec['reads_written'] == rec['reads_target'] or rec['fastq_eof'] is True,
              'strict recount == reads written == recorded recount': v['reads_recounted_now'] == rec['reads_written'] == rec['reads_recounted'],
              'SHA-256 now == _COMPLETE.json == manifest FASTQ_VALIDATED': v['sha256_now'] == rec['fastq_sha256'] and bool(mv) and mv[-1]['fastq_sha256'] == rec['fastq_sha256'] and v['bytes_ok'],
              'kallisto exit code 0': rec['kallisto_rc'] == 0,
              'kallisto n_processed == reads written': rec['n_processed'] == rec['reads_written'],
              'outputs non-empty (abundance.tsv, abundance.h5, run_info.json)': all(sz > 0 for sz in files.values()),
              'abundance.tsv has 15,314 targets': v['abundance_rows'] == 15314,
              'manifest last state COMPLETE': v['manifest_last'] == 'COMPLETE'}
    v['checks'] = checks; v['verdict'] = 'PASS' if all(checks.values()) else 'FAIL: ' + '; '.join(k for k, ok in checks.items() if not ok)
    return v
if MODE == 'run':
    plog(f'PILOT pass 1 (run) - fetch_and_quant.py sha256 {m.sha256_file(FQS)[:16]}, N={m.N}')
    K, idx, tmpd, summ, jobs = setup()
    with open(os.path.join(PD, 'pilot_selection.json'), 'w', encoding='utf-8') as f: json.dump([(l, r, w) for l, r, w in jobs], f, indent=1)
    for label, r, why in jobs:
        t0 = time.time(); plog(f"start {label} A {r['run_accession']} ({r.get('sample_title', '')}; {why})")
        run_job(K, idx, tmpd, summ, label, r)
        plog(f"end   {label} A {r['run_accession']} after {(time.time() - t0) / 60:.1f} min; manifest: {verify(label, r)['manifest_states']}")
    snap = snapshot(jobs)
    with open(os.path.join(PD, 'snapshot_after_pass1.json'), 'w', encoding='utf-8') as f: json.dump(snap, f, indent=1)
    plog(f'pass 1 finished; {len(snap)} files snapshotted; manifest rows {len(manifest_rows())}')
elif MODE == 'restart':
    plog('PILOT pass 2 (restart in a new process)')
    K, idx, tmpd, summ, jobs = setup()
    with open(os.path.join(PD, 'pilot_selection.json'), encoding='utf-8') as f: sel = [r['run_accession'] for _, r, _ in json.load(f)]
    jobs = [j for j in jobs if j[1]['run_accession'] in sel]
    with open(os.path.join(PD, 'snapshot_after_pass1.json'), encoding='utf-8') as f: snap1 = json.load(f)
    before = snapshot(jobs); n0 = len(manifest_rows())
    calls = []; real_sha = m.sha256_file
    def counting_sha(p): calls.append(os.path.basename(p)); return real_sha(p)
    m.sha256_file = counting_sha                          # instrumentation only: proves the skip re-verified each checksum
    restart = {}
    for label, r, why in jobs:
        k = len(calls); run_job(K, idx, tmpd, summ, label, r); restart[r['run_accession']] = [c for c in calls[k:]]
    m.sha256_file = real_sha
    after = snapshot(jobs); n1 = len(manifest_rows())
    rows = []
    for label, r, why in jobs:
        run = r['run_accession']; fqn = f'{run}_R1_subsetA.fastq.gz'
        try: v = verify(label, r)
        except Exception as e: v = {'study': label, 'run': run, 'verdict': 'FAIL: verification error ' + repr(e)[:150]}
        mine = [k for k in after if run in k]
        unchanged = all(before.get(k) == after.get(k) == snap1.get(k) for k in mine) and len(mine) > 0
        new_rows = [d for d in manifest_rows()[n0:] if d['run'] == run]
        if v.get('complete_marker'):
            ok = unchanged and not new_rows and fqn in restart[run]
            v['restart'] = ('PASS: skipped after SHA-256 re-verification; files unchanged; no manifest row' if ok else
                            f'FAIL: unchanged={unchanged} new_manifest_rows={[d["state"] for d in new_rows]} sha256_checked={fqn in restart[run]}')
        else:
            v['restart'] = f'not complete in pass 1 -> retried in pass 2; new states: {[d["state"] for d in new_rows]}'
        rows.append((why, v))
    cols = ['study', 'run', 'sample_alias', 'title', 'reads_target', 'reads_written', 'reads_recounted_now', 'p_pseudoaligned', 'kallisto_rc', 'n_processed',
            'abundance_rows', 'sha256_recorded', 'sha256_now', 'complete_marker', 'manifest_states', 'verdict', 'restart', 'run_info_parsed_with']
    TS = time.strftime('%Y%m%d_%H%M%S'); RT = os.path.join(PD, f'pilot_report_{TS}.tsv'); RX = os.path.join(PD, f'pilot_report_{TS}.txt')
    with open(RT, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\tselection_reason\n')
        for why, v in rows: f.write('\t'.join(str(v.get(c, '')) for c in cols) + '\t' + why + '\n')
    with open(RX, 'w', encoding='utf-8') as f:
        f.write(f'PILOT REPORT {time.strftime("%Y-%m-%d %H:%M")} | fetch_and_quant.py sha256 {real_sha(FQS)} | N = {m.N} reads (Read 1, subset A) per run\n')
        f.write(f'manifest rows before restart {n0}, after restart {n1}\n\n')
        for why, v in rows:
            f.write(f"{v['study']} {v['run']} ({v.get('sample_alias')}, {v.get('title')}): {v['verdict']}\n  reason selected: {why}\n")
            for k in ('reads_target', 'reads_written', 'reads_recounted_now', 'p_pseudoaligned', 'kallisto_rc', 'n_processed', 'abundance_rows', 'sha256_now', 'complete_marker', 'manifest_states', 'restart'):
                f.write(f'  {k}: {v.get(k, "")}\n')
            for k, ok in v.get('checks', {}).items(): f.write(f"  [{'PASS' if ok else 'FAIL'}] {k}\n")
            f.write('\n')
        allok = all(v['verdict'] == 'PASS' and v['restart'].startswith('PASS') for _, v in rows) and len(rows) == len(SELECT)
        f.write('PILOT RESULT: ' + ('PASS - all runs complete, verified, and correctly skipped on restart' if allok else 'NOT PASSED - see entries above') + '\n')
    plog('pass 2 finished; report: ' + RX + ' | ' + ('PILOT RESULT: PASS' if allok else 'PILOT RESULT: NOT PASSED'))
