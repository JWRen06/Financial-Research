"""Replay the paper's finite numerical certificates in a disposable copy."""
from pathlib import Path
import hashlib, json, os, shutil, subprocess, sys, tempfile, time

def main():
    here=Path(__file__).resolve().parent
    plan=json.loads((here/'plan.json').read_text(encoding='utf-8'))
    output=here/'run_output'
    output.mkdir(exist_ok=True)
    env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
    rows=[]
    with tempfile.TemporaryDirectory(prefix='certificate_',dir=output) as tmp:
        work=Path(tmp)/'inputs'
        shutil.copytree(here/'inputs',work)
        for job in plan['jobs']:
            script=work/job
            print('Running',script.name,flush=True)
            start=time.perf_counter()
            proc=subprocess.run([sys.executable,'-X','utf8',str(script)],cwd=script.parent,env=env,capture_output=True,text=True,encoding='utf-8',timeout=900)
            (output/(script.stem+'.log')).write_text(proc.stdout+'\n'+proc.stderr,encoding='utf-8')
            rows.append({'program':job,'source_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
                         'exit_code':proc.returncode,'seconds':time.perf_counter()-start})
            (output/'verification.json').write_text(json.dumps({'python':sys.version,'jobs':rows,'all_passed':all(r['exit_code']==0 for r in rows)},ensure_ascii=False,indent=2),encoding='utf-8')
            if proc.returncode:
                raise RuntimeError('Certificate failed; inspect run_output/'+script.stem+'.log')
    print('Passed',len(rows),'finite certificate checks; no optimizer search rerun.',flush=True)

if __name__=='__main__': main()
