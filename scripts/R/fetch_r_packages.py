r"""Fallback downloader for R packages when R itself cannot reach the repositories (stdlib only).

Used by scripts/R/setup_r_packages.R ONLY if R's own HTTPS access fails. It reads the official repository indexes
(PACKAGES files) for the R/Bioconductor versions given, computes the dependency closure of the requested packages
(Depends, Imports, LinkingTo; minus base packages and packages already installed in a sufficient version), downloads
each package file from the official repository, verifies its MD5 against the index, and writes an install plan in
dependency order. R then installs the local files with install.packages(<file>, repos = NULL, lib = <project library>).

Usage: python fetch_r_packages.py <project_dir> <bioc_version> <r_minor e.g. 4.6> <installed.tsv> <pkg1,pkg2,...> [--mirrors URL1,URL2,...] [--base-url URL]
  --mirrors: official Bioconductor mirrors tried after bioconductor.org; when two sources are reachable, a package is used only if
             both indexes list the same version with the same MD5 (otherwise rejected).
  installed.tsv: Package<TAB>Version (everything already in .libPaths(), written by R)
Outputs (all under <project>/envs/R-pkgcache/<bioc_version>/):
  <repo>/PACKAGES (downloaded index), files/<package file>, install_plan.tsv, fetch_log.txt
Nothing outside the project folder is written. Nothing is deleted."""
import os, sys, re, hashlib, urllib.request, urllib.error, time, json
BASE = ['R','base','compiler','datasets','graphics','grDevices','grid','methods','parallel','splines','stats','stats4','tcltk','tools','utils']
def log(m):
    s = time.strftime('%Y-%m-%dT%H:%M:%S') + '  ' + m; print(s, flush=True)
    with open(LOG, 'a', encoding='utf-8') as f: f.write(s + '\n')
def get(url, timeout=120, tries=4):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Apis-project-R-fetch/1.0'}), timeout=timeout) as r: return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 410): raise RuntimeError(f'failed {url}: HTTP {e.code}')
            last = e; log(f'  retry {i+1} {url}: {e}'); time.sleep(3 * (i + 1))
        except Exception as e:
            last = e; log(f'  retry {i+1} {url}: {e}'); time.sleep(3 * (i + 1))
    raise RuntimeError(f'failed {url}: {last}')
def parse_dcf(txt):
    recs, cur, key = [], {}, None
    for line in txt.splitlines():
        if not line.strip():
            if cur: recs.append(cur); cur, key = {}, None
            continue
        if line[0] in ' \t' and key: cur[key] += ' ' + line.strip(); continue
        k, _, v = line.partition(':'); key = k.strip(); cur[key] = v.strip()
    if cur: recs.append(cur)
    return recs
def vkey(v): return [int(x) for x in re.split(r'[.-]', v) if x.isdigit()]
def vcmp(a, b):
    a, b = vkey(a), vkey(b); n = max(len(a), len(b)); a += [0] * (n - len(a)); b += [0] * (n - len(b))
    return (a > b) - (a < b)
def deps(rec):
    out = []
    for f in ('Depends', 'Imports', 'LinkingTo'):
        for item in rec.get(f, '').split(','):
            item = item.strip()
            if not item: continue
            m = re.match(r'^([A-Za-z0-9.]+)\s*(?:\(\s*([<>=]+)\s*([^)\s]+)\s*\))?', item)
            if m: out.append((m.group(1), m.group(2), m.group(3)))
    return out
def satisfies(have, op, need):
    if not op: return True
    c = vcmp(have, need)
    return {'>=': c >= 0, '>': c > 0, '==': c == 0, '<=': c <= 0, '<': c < 0}.get(op, True)
def main():
    global LOG
    a = sys.argv[1:]
    proj, bioc, rminor, inst_tsv, pkgs = a[0], a[1], a[2], a[3], [p for p in a[4].split(',') if p]
    base_url = a[a.index('--base-url') + 1].rstrip('/') if '--base-url' in a else None
    cache = os.path.join(proj, 'envs', 'R-pkgcache', bioc); os.makedirs(os.path.join(cache, 'files'), exist_ok=True)
    LOG = os.path.join(cache, 'fetch_log.txt')
    # Bioconductor sources: official site first, then official mirrors (https://bioconductor.org/BioC_mirrors.csv) given with --mirrors
    bases = [base_url] if base_url else ['https://bioconductor.org']
    if '--mirrors' in a: bases += [m.rstrip('/') for m in a[a.index('--mirrors') + 1].split(',') if m and m.rstrip('/') not in bases]
    C = base_url + '/cran' if base_url else 'https://cloud.r-project.org'
    bioc_repos = [('BioCsoft', f'/packages/{bioc}/bioc/bin/windows/contrib/{rminor}', 'win.binary'),
                  ('BioCann', f'/packages/{bioc}/data/annotation/src/contrib', 'source'),
                  ('BioCexp', f'/packages/{bioc}/data/experiment/src/contrib', 'source')]
    cran_repos = [('CRAN', f'{C}/bin/windows/contrib/{rminor}', 'win.binary')]
    bioc_src = [('BioCsoftsrc', f'/packages/{bioc}/bioc/src/contrib', 'source')]   # Bioconductor packages without a Windows binary
    cran_src = [('CRANsrc', f'{C}/src/contrib', 'source')]
    log(f'start: bioc {bioc}, R {rminor}, requested {pkgs}; Bioconductor sources tried in order: {bases}')
    index = {}
    def load_index(name, url):
        txt = get(url + '/PACKAGES', tries=2).decode('utf-8', 'replace'); return {r['Package']: r for r in parse_dcf(txt) if r.get('Package')}, txt
    def add(name, typ, primary_url, recs, cross=None, cross_base=None):
        n = 0
        for p, r in recs.items():
            if p in index: continue
            x = cross.get(p) if cross else None
            xs = ('agree: ' + cross_base) if x and x.get('Version') == r.get('Version') and x.get('MD5sum') == r.get('MD5sum') else \
                 ('DISAGREE: ' + cross_base) if x and x.get('Version') == r.get('Version') else 'single source'
            index[p] = dict(r, _repo=name, _url=primary_url, _type=typ, _cross=xs); n += 1
        log(f'index {name}: {primary_url} ({n} new packages; cross-check {"with " + cross_base if cross else "none"})')
    def bioc_group(repos):
        for name, path, typ in repos:
            got = []
            for b in bases:
                try: recs, txt = load_index(name, b + path)
                except Exception as e: log(f'  {name} at {b}: {e}'); continue
                got.append((b, recs, txt))
                if len(got) == 2: break
            if not got: log(f'index {name} unavailable at every source'); continue
            d = os.path.join(cache, name); os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, 'PACKAGES'), 'w', encoding='utf-8') as f: f.write(got[0][2])
            add(name, typ, got[0][0] + path, got[0][1], got[1][1] if len(got) > 1 else None, got[1][0] if len(got) > 1 else None)
    def cran_group(repos):
        for name, url, typ in repos:
            try: recs, txt = load_index(name, url)
            except Exception as e: log(f'index {name} unavailable: {e}'); continue
            add(name, typ, url, recs)
    bioc_group(bioc_repos); cran_group(cran_repos); bioc_group(bioc_src); cran_group(cran_src)
    if not index: log('FATAL: no repository index could be downloaded'); sys.exit(3)
    have = {}
    with open(inst_tsv, encoding='utf-8') as f:
        for line in f.read().splitlines()[1:]:
            p, v = (line.split('\t') + [''])[:2]
            if p and (p not in have or vcmp(v, have[p]) > 0): have[p] = v
    need, order, missing, seen = {}, [], [], set()
    def visit(p, op=None, ver=None, forced=False):
        # requested (forced) packages are always (re)installed and their dependencies always resolved,
        # so an installed-but-broken package (e.g. clusterProfiler without GO.db) gets its missing dependencies
        if p in BASE: return
        if p in seen: return
        if p in have and satisfies(have[p], op, ver) and not forced:
            # installed and sufficient: not reinstalled, but ITS dependencies are still checked, because an installed
            # package can lack one of its own dependencies (observed: stringr present, stringi missing)
            seen.add(p); r = index.get(p)
            if r is not None:
                for q, o, v in deps(r): visit(q, o, v)
            return
        seen.add(p)
        r = index.get(p)
        if r is None: missing.append(p); return
        if op and not satisfies(r['Version'], op, ver): missing.append(f'{p} ({op} {ver}; repo has {r["Version"]})'); return
        for q, o, v in deps(r): visit(q, o, v)
        need[p] = r; order.append(p)
    for p in pkgs: visit(p, forced=True)
    log(f'install plan: {len(order)} packages in dependency order; unresolved: {missing or "none"}')
    with open(os.path.join(cache, 'unresolved.txt'), 'w', encoding='utf-8') as f: f.write('')   # never leave a stale list from an earlier run
    plan = []
    for p in order:
        r = need[p]; ext = '.zip' if r['_type'] == 'win.binary' else '.tar.gz'
        fn = f"{p}_{r['Version']}{ext}"; dest = os.path.join(cache, 'files', fn); md5 = r.get('MD5sum', '')
        ok = os.path.exists(dest) and md5 and hashlib.md5(open(dest, 'rb').read()).hexdigest() == md5
        if r['_cross'].startswith('DISAGREE'):
            log(f'REJECTED {fn}: MD5 differs between {r["_url"]} and {r["_cross"]}'); missing.append(p + ' (mirror MD5 disagreement)'); continue
        if not ok:
            if os.path.exists(dest): os.replace(dest, dest + '.bad_' + time.strftime('%Y%m%d%H%M%S'))
            b = None
            rel = r['_url'][len(next((x for x in bases if r['_url'].startswith(x)), '')):] if any(r['_url'].startswith(x) for x in bases) else None
            srcs = [r['_url']] + ([x + rel for x in bases if x + rel != r['_url']] if rel else [])   # same file on the official mirrors
            used = None
            for su in srcs:
                try: b = get(f"{su}/{fn}", timeout=1800); used = su; break
                except Exception as e: log(f'  download failed from {su}: {e}')
            if b is None: missing.append(p + ' (download failed at every source)'); continue
            got = hashlib.md5(b).hexdigest()
            if md5 and got != md5: log(f'MD5 MISMATCH {fn}: index {md5} got {got} - not used'); missing.append(p + ' (md5 mismatch)'); continue
            with open(dest + '.part', 'wb') as f: f.write(b)
            os.replace(dest + '.part', dest)
            log(f'downloaded {fn} ({len(b)/1e6:.1f} MB) md5 {"verified" if md5 else "NOT IN INDEX"} from {used}; index cross-check: {r["_cross"]}')
        else: log(f'cached {fn} md5 verified')
        plan.append((p, r['Version'], r['_repo'], r['_type'], dest, md5, r['_url'], r['_cross']))
    with open(os.path.join(cache, 'install_plan.tsv'), 'w', encoding='utf-8') as f:
        f.write('package\tversion\trepo\ttype\tfile\tmd5\tsource\tcross_check\n')
        for row in plan: f.write('\t'.join(row) + '\n')
    with open(os.path.join(cache, 'unresolved.txt'), 'w', encoding='utf-8') as f: f.write('\n'.join(missing) + ('\n' if missing else ''))
    log(f'done: plan {len(plan)} files -> {os.path.join(cache, "install_plan.tsv")}; unresolved {len(missing)}')
    sys.exit(0 if not missing else 4)
if __name__ == '__main__': main()
