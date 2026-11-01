from pathlib import Path
from decimal import Decimal as D
import json, csv, platform, hashlib, sys, time
from hidden_system import HiddenSystem
from audit_interval import audit
from fractions import Fraction as Q

ROOT=Path(__file__).resolve().parent
rows=[]
def evaluate(group, label, system, costs, eps, **meta):
    # Centers are deterministic forward responses. Error boxes represent
    # measurement resolution, rather than invented repeated random samples.
    values=[system.reset_and_report(k,6) for k in costs]
    req={'rho':'1','costs':costs,'reports':[[str(v) for v in ch] for ch in values],
         'error_bounds':str(eps),'max_beta_width':'.02','threshold':'.7'}
    res=audit(req)
    if group!='fault':
        assert res['status']!='MODEL_INCOMPATIBLE'
        low,high=map(Q,res['beta_outer']['exact'])
        assert low<=Q(4,5)<=high
    row=dict(group=group,label=label,costs=costs,epsilon=str(eps),
             peak=float(max(abs(v) for ch in values for v in ch)),
             status=res['status'],width=res['beta_width'],
             interval=res['beta_outer']['outward_float'] if res['beta_outer'] else None,**meta)
    rows.append(row)
    (ROOT/'continuous_inputs').mkdir(exist_ok=True)
    (ROOT/'continuous_inputs'/f'{len(rows):03d}.json').write_text(json.dumps({'input':req,'output':res},indent=2),encoding='utf-8')
    return row

start=time.time()
for j in range(25):
    gap=10**(-3+j/8)
    cost=str(D(1)+D(str(gap)))
    for eps in ('1e-10','1e-8','1e-6'):
        evaluate('gap',f'{j}:{eps}',HiddenSystem(),['1',cost],eps,gap=gap)
for s in ('.8','.95','.99','.999','.9999'):
    system=HiddenSystem(nodes=(s,),weights=('1',))
    peak=max(abs(v) for k in ['1','2'] for v in system.reset_and_report(k,6))
    for mode,eps in [('absolute',D('1e-8')),('relative',D('1e-6')*peak)]:
        evaluate('persistence',s+':'+mode,system,['1','2'],eps,persistence=float(s),error_mode=mode)
for eps in ('1e-12','1e-10','1e-8','1e-6','1e-4'):
    evaluate('noise',eps,HiddenSystem(),['1','2'],eps)
for fault in ('different_poles','two_poles'):
    for eps in ('1e-12','1e-8','1e-4'):
        evaluate('fault',fault+':'+eps,HiddenSystem(misspec=fault),['1','2'],eps,fault=fault)
(ROOT/'continuous_results.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
with (ROOT/'continuous_results.csv').open('w',encoding='utf-8-sig',newline='') as f:
    keys=sorted(set().union(*(r.keys() for r in rows)))
    writer=csv.DictWriter(f,keys);writer.writeheader();writer.writerows(rows)
summary={'count':len(rows),'legal_boxes_checked':sum(r['group']!='fault' for r in rows),
         'fault_rejected':sum(r['status']=='MODEL_INCOMPATIBLE' for r in rows),
         'seconds':time.time()-start,'python':sys.version,'platform':platform.platform(),
         'center_design':'deterministic exact-response approximations; no stochastic sampling',
         'hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'audit_interval.py',ROOT/'hidden_system.py',Path(__file__)]}}
(ROOT/'continuous_verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
for r in rows:
    if r['group']!='gap': print(r['label'],r['status'],r['width'],r['peak'])
