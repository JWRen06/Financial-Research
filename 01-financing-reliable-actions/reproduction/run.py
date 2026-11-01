from pathlib import Path
import os,sys,subprocess,json,time,hashlib,platform
HERE=Path(__file__).resolve().parent
out=HERE/'run_output';out.mkdir(exist_ok=True)
env=dict(os.environ,FINANCING_OUTPUT=str(out),PYTHONUTF8='1',PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
report={'python':sys.version,'platform':platform.platform(),'jobs':[]}
for name in ['financing_case.py','verify_messages.py']:
    p=HERE/'experiments'/name;t=time.perf_counter()
    r=subprocess.run([sys.executable,str(p)],cwd=HERE,env=env,capture_output=True,text=True,encoding='utf-8')
    (out/(name+'.log')).write_text(r.stdout+'\n'+r.stderr,'utf-8')
    report['jobs'].append({'script':name,'exit_code':r.returncode,'seconds':time.perf_counter()-t,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (out/'run_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),'utf-8')
    print(name,r.returncode,flush=True)
    if r.returncode:raise RuntimeError(name)
for name in ['case_results.json','独立消息生成器核验.json','估计器归一化核验.json']:
    assert (out/name).read_bytes()==(HERE/'results'/name).read_bytes(),name
print('All runs passed and JSON outputs match this release.',flush=True)
