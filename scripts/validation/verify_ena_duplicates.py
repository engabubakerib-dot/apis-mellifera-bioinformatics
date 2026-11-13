r"""Verify that two ENA runs contain the SAME reads (duplicate listing) - Python stdlib only; nothing is written to disk.

Compares decompressed FASTQ records, not gzip files: two different gzip containers of identical reads give different
file MD5s but identical record digests. Per record it compares the sequence, the quality string and the ORIGINAL Illumina
read ID (the header without the run-specific '@SRRxxxx.N ' prefix). Streams both files in parallel, record by record.

Usage:
  python verify_ena_duplicates.py <pairs.tsv> <out.tsv> [--max-records N]
    pairs.tsv: columns library, run_a, url_a, run_b, url_b (url may be https://... or a local path)
    --max-records N: compare only the first N records of each pair (0 or absent = whole files)
Output per pair: records compared, records identical, first differing record (or -), end-of-file agreement,
SHA-256 over all compared (illumina_id, sequence, quality) records for each run, verdict IDENTICAL_CONTENT / DIFFERENT / INCOMPLETE."""
import sys, gzip, hashlib, urllib.request, os, io, time
def open_fastq(src):
    raw = open(src, 'rb') if os.path.exists(src) else urllib.request.urlopen(urllib.request.Request(src, headers={'User-Agent': 'Apis-dup-check/1.0'}), timeout=300)
    return gzip.GzipFile(fileobj=raw)          # multi-member gzip handled by GzipFile
def records(f):
    while True:
        h = f.readline()
        if not h: return
        s, p, q = f.readline(), f.readline(), f.readline()
        if not q: raise ValueError('truncated FASTQ record')
        yield h.rstrip(b'\r\n').split(b' ', 1)[-1] if h.startswith(b'@SRR') or h.startswith(b'@ERR') or h.startswith(b'@DRR') else h.rstrip(b'\r\n'), s.rstrip(b'\r\n'), q.rstrip(b'\r\n')
def compare(a, b, maxn):
    ha, hb = hashlib.sha256(), hashlib.sha256(); n = same = 0; first = '-'; eof_a = eof_b = False
    with open_fastq(a) as fa, open_fastq(b) as fb:
        ia, ib = records(fa), records(fb)
        while not maxn or n < maxn:
            ra, rb = next(ia, None), next(ib, None)
            if ra is None or rb is None: eof_a, eof_b = ra is None, rb is None; break
            n += 1; ha.update(b'\t'.join(ra) + b'\n'); hb.update(b'\t'.join(rb) + b'\n')
            if ra == rb: same += 1
            elif first == '-': first = str(n)
    complete = (eof_a and eof_b) or (maxn and n == maxn)
    verdict = 'IDENTICAL_CONTENT' if complete and same == n and n > 0 and eof_a == eof_b else ('DIFFERENT' if same != n or eof_a != eof_b else 'INCOMPLETE')
    return n, same, first, f'{eof_a}/{eof_b}', ha.hexdigest(), hb.hexdigest(), verdict
if __name__ == '__main__':
    pairs, out = sys.argv[1], sys.argv[2]; maxn = int(sys.argv[sys.argv.index('--max-records') + 1]) if '--max-records' in sys.argv else 0
    L = open(pairs, encoding='utf-8').read().strip().split('\n'); h = L[0].split('\t')
    with open(out, 'w', encoding='utf-8') as f:
        f.write('library\trun_a\trun_b\tmode\trecords_compared\trecords_identical\tfirst_difference\teof_a/eof_b\tsha256_records_a\tsha256_records_b\tverdict\ttime\n')
        for line in L[1:]:
            d = dict(zip(h, line.split('\t')))
            try: r = compare(d['url_a'], d['url_b'], maxn)
            except Exception as e: r = (0, 0, '-', '-', '-', '-', 'ERROR ' + str(e)[:80])
            f.write('\t'.join([d['library'], d['run_a'], d['run_b'], f'first {maxn}' if maxn else 'whole files'] + [str(x) for x in r] + [time.strftime('%Y-%m-%dT%H:%M:%S')]) + '\n'); f.flush()
            print(d['library'], r[-1], r[0], 'records', flush=True)
