"""Mock resume/checkpoint tests. Usage: python tests.py <empty_base_folder>   (returns exit code 1 on any failure)"""
import os,sys,json,subprocess,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from harness import *
def J(*p): return os.path.join(*p)
def run_all(base,report=print):
    res=[]
    def check(name,cond): res.append((name,bool(cond))); report(('PASS ' if cond else 'FAIL ')+name)
    c=setup(base); R=c.R; K=J(R,'results/rawreads/kallisto'); FQ=J(R,'data/raw/reads_subsample')
    # legacy partial FASTQ under the final name (SRR13938890 situation) + stale output folder without _COMPLETE.json
    import gzip
    with open(J(FQ,'RUN2_R1_subsetA.fastq.gz'),'wb') as f: f.write(gzip.compress(b'@x\nAC\n+\nII\n')[:-5])
    os.makedirs(J(K,'GSE168740/A/RUN3'),exist_ok=True)
    with open(J(K,'GSE168740/A/RUN3/abundance.tsv'),'w') as f: f.write('stale partial')
    sz=os.path.getsize(J(c.src,'RUN1.fq.gz')); szs=os.path.getsize(J(c.src,'SHORT.fq.gz'))
    run(c,breaks={'RUN1':sz//3,'SHORT':-(szs//2)}); st,_=state(c)
    check('T1 connection reset mid-stream -> DOWNLOAD_FAILED, no _COMPLETE.json', st[('GSE168740','A','RUN1')]=='DOWNLOAD_FAILED' and not os.path.exists(J(K,'GSE168740/A/RUN1/_COMPLETE.json')))
    w=[d for d in manifest_rows(c) if d['run']=='RUN1' and d['subset']=='A' and d['state']=='DOWNLOAD_FAILED']
    check('T1 the reset happened after reads were written (true mid-stream cut)', w and int(w[0]['reads_written'])>0)
    check('T1 partial download kept as .part', os.path.exists(J(FQ,'RUN1_R1_subsetA.fastq.gz.part')))
    check('T1 legacy partial FASTQ kept aside as .prev; run completed', os.path.exists(J(FQ,'RUN2_R1_subsetA.fastq.gz.prev')) and st[('GSE168740','A','RUN2')]=='COMPLETE')
    check('T1 stale folder without _COMPLETE.json kept as .failed; run redone', os.path.exists(J(K,'GSE168740/A/RUN3.failed/abundance.tsv')) and st[('GSE168740','A','RUN3')]=='COMPLETE')
    check('T1 server closes early without error -> DOWNLOAD_FAILED', st[('GSE178742','A','SHORT')]=='DOWNLOAD_FAILED')
    run(c); st,_=state(c)
    with open(J(K,'GSE178742/A/SHORT/_COMPLETE.json'),encoding='utf-8') as f: shortrec=json.load(f)
    check('T1b genuinely short intact source completes with fastq_eof=True', st[('GSE178742','A','SHORT')]=='COMPLETE' and shortrec['fastq_eof'] is True)
    check('T1b RUN1 resumes after reset and completes (old partial kept as .part.prev)', st[('GSE168740','A','RUN1')]=='COMPLETE' and os.path.exists(J(FQ,'RUN1_R1_subsetA.fastq.gz.part.prev')))
    add_run(c,'ena_PRJNA713798_runs.tsv','RUN4','RUN1')
    done=lambda k: any(x in k for x in ('/RUN1/','/RUN2/','/RUN3/','/SHORT/')) or k in ('fq/RUN1_R1_subsetA.fastq.gz','fq/RUN2_R1_subsetA.fastq.gz','fq/RUN3_R1_subsetA.fastq.gz')
    before=snap(c); run(c,fakek='fail'); st,_=state(c); after=snap(c)
    check('T2 validated runs skipped: files byte-identical, mtimes unchanged', all(before[k]==after[k] for k in before if done(k)))
    check('T2 kallisto exit 1 -> QUANT_FAILED, no final folder', st[('GSE168740','A','RUN4')]=='QUANT_FAILED' and not os.path.exists(J(K,'GSE168740/A/RUN4')))
    check('T2 failed output kept as .inprogress', os.path.isdir(J(K,'GSE168740/A/RUN4.inprogress')))
    run(c,fakek='partial'); st,_=state(c)
    check('T3 kallisto killed mid-write (rc 137, partial abundance) -> QUANT_FAILED', st[('GSE168740','A','RUN4')]=='QUANT_FAILED')
    check('T3 earlier .inprogress kept as .inprogress.failed', os.path.isdir(J(K,'GSE168740/A/RUN4.inprogress.failed')))
    before=snap(c); run(c); st,_=state(c); after=snap(c)
    check('T4 resume completes the incomplete run', st[('GSE168740','A','RUN4')]=='COMPLETE')
    check('T4 previously validated outputs unchanged', all(before[k]==after[k] for k in before if k in after and done(k)))
    fq=J(FQ,'RUN2_R1_subsetA.fastq.gz'); ab=J(K,'GSE168740/A/RUN2/abundance.tsv'); h0=sha(ab)
    with open(fq,'ab') as f: f.write(b'tamper')
    run(c); st,_=state(c)
    check('T5 altered FASTQ of a completed run -> REVALIDATION_FAILED, outputs not overwritten', st[('GSE168740','A','RUN2')]=='REVALIDATION_FAILED' and sha(ab)==h0)
    _,a=state(c); run(c); new=manifest_rows(c)[a:]
    check('T6 clean rerun touches nothing except the run flagged for review', all(d['run']=='RUN2' for d in new))
    add_run(c,'ena_PRJNA713798_runs.tsv','RUN5','RUN3')
    run(c,fakek='fail'); fqp=J(FQ,'RUN5_R1_subsetA.fastq.gz'); m0=os.path.getmtime(fqp)
    _,a=state(c); run(c); st,_=state(c); new=manifest_rows(c)[a:]
    check('T7 validated FASTQ reused after QUANT_FAILED (no re-download; file untouched)', st[('GSE168740','A','RUN5')]=='COMPLETE' and not any(d['run']=='RUN5' and d['state']=='DOWNLOAD_STARTED' for d in new) and os.path.getmtime(fqp)==m0)
    add_run(c,'ena_PRJNA713798_runs.tsv','RUN6','RUN3')
    run(c,fakek='fail'); fqp=J(FQ,'RUN6_R1_subsetA.fastq.gz')
    with open(fqp,'ab') as f: f.write(b'x')
    run(c); st,_=state(c)
    check('T8 altered validated FASTQ not reused: kept as .prev, downloaded again, completes', st[('GSE168740','A','RUN6')]=='COMPLETE' and os.path.exists(fqp+'.prev'))
    return res
def run_kill(base,report=print):
    """Hard-kill (power-loss simulation) of the whole process tree while kallisto is running."""
    res=[]
    def check(name,cond): res.append((name,bool(cond))); report(('PASS ' if cond else 'FAIL ')+name)
    c=setup(base); K=J(c.R,'results/rawreads/kallisto')
    code=f"import sys; sys.path.insert(0,{os.path.dirname(os.path.abspath(__file__))!r}); from harness import *; c=Ctx({c.base!r}); run(c,fakek='hang')"
    kw={'creationflags':subprocess.CREATE_NEW_PROCESS_GROUP} if os.name=='nt' else {'start_new_session':True}
    p=subprocess.Popen([sys.executable,'-c',code],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,**kw)
    t0=time.time()
    while time.time()-t0<90 and not os.path.exists(J(K,'GSE168740/A/RUN1.inprogress/abundance.tsv')): time.sleep(0.5)
    time.sleep(1); kill_tree(p); time.sleep(2)
    st,_=state(c)
    check('K1 killed during kallisto: last state FASTQ_VALIDATED, no final folder (no false COMPLETE)', st.get(('GSE168740','A','RUN1'))=='FASTQ_VALIDATED' and not os.path.exists(J(K,'GSE168740/A/RUN1')))
    run(c); st,_=state(c)
    check('K2 restart completes the run; killed attempt kept as .inprogress.failed', st[('GSE168740','A','RUN1')]=='COMPLETE' and os.path.isdir(J(K,'GSE168740/A/RUN1.inprogress.failed')))
    check('K3 FASTQ from the killed attempt reused (no second download)', sum(1 for d in manifest_rows(c) if d['run']=='RUN1' and d['subset']=='A' and d['state']=='DOWNLOAD_STARTED')==1)
    return res
if __name__=='__main__':
    base=os.path.abspath(sys.argv[1]); os.makedirs(base,exist_ok=True)
    r=run_all(J(base,'mock'))+run_kill(J(base,'kill'))
    print(f"{sum(ok for _,ok in r)}/{len(r)} passed"); sys.exit(0 if all(ok for _,ok in r) else 1)
