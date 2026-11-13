from harness import *
import subprocess,signal,time,sys,os
setup()
p=subprocess.Popen([sys.executable,'-c','from harness import *; run(fakek="hang")'],start_new_session=True)
time.sleep(4); os.killpg(p.pid,signal.SIGKILL); p.wait()
st,_=state(); print('after hard kill, last states:',st)
print('final dir exists:',os.path.exists(os.path.join(R,'results/rawreads/kallisto/GSE168740/A/RUN1')),' inprogress exists:',os.path.isdir(os.path.join(R,'results/rawreads/kallisto/GSE168740/A/RUN1.inprogress')))
run(); st,_=state(); print('after restart:',st[('GSE168740','A','RUN1')])
print('killed attempt kept as:',sorted(f for f in os.listdir(os.path.join(R,'results/rawreads/kallisto/GSE168740/A')) if 'RUN1' in f))
