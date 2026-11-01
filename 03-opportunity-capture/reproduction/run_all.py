"""Reproduce the rewritten paper in a working copy; preserve delivered archives."""
from pathlib import Path
import sys,os,subprocess,shutil,json,hashlib,time,platform
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'run_output';OUT.mkdir(exist_ok=True)
WORK=OUT/'material';WORK.mkdir(exist_ok=True)
shutil.copytree(ROOT/'experiments',WORK/'experiments',dirs_exist_ok=True)
for name in ['study.py','check_simplified_stopping.py']:shutil.copy2(ROOT/name,WORK/name)
env=dict(os.environ,PYTHONUTF8='1',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
jobs=['experiments/second_paper_bounded_revision.py','experiments/second_paper_bounded_independent.py','experiments/second_paper_bounded_deployment.py','experiments/application_contract_bridge.py','study.py','check_simplified_stopping.py']
results=[]
for job in jobs:
    print('Running '+job,flush=True);start=time.perf_counter()
    script=WORK/job
    p=subprocess.run([sys.executable,str(script)],cwd=WORK,capture_output=True,text=True,encoding='utf-8',env=env)
    log=OUT/(Path(job).stem+'.log');log.write_text(p.stdout+'\n'+p.stderr,'utf-8')
    row={'program':job,'exit_code':p.returncode,'seconds':time.perf_counter()-start,'sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
    results.append(row)
    (OUT/'run_report.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'jobs':results},indent=2),'utf-8')
    if p.returncode:raise RuntimeError(job+' failed; inspect '+str(log))
print('All programs completed. New results and figures: '+str(WORK/'focused_results'))
