"""Independent Decimal reports, physical Stein energies, noise splits and hashes.
Does not import SciPy, NumPy, proposal code, constructor or interval engine.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json,hashlib
HERE=Path(__file__).resolve().parent
def dec(x):
    f=F(x);return Decimal(f.numerator)/Decimal(f.denominator)
def inside(x,b):assert dec(b[0])<=x<=dec(b[1]),(x,b)
def close(x,y):assert abs(x-y)<Decimal('1e-160'),(x,y)
def root(beta,k):
    v=k*(1+beta)+beta/15;return 2*k/(v+(v*v-4*k*k*beta).sqrt())
def az(e,s,k):return (1-e)/(1+15*k*(1-e)*(1-s))
def pair_energy(e,e0,s,t,b):
    aa=az(e,s,1);bb=az(e,t,1);cs=az(e,s,1)-az(e0,s,1);ct=az(e,t,1)-az(e0,t,1)
    y11=1/(1-b*s*t);y01=b*aa*t*y11/(1-b*e*t);y10=b*bb*s*y11/(1-b*e*s)
    y00=b*(e*bb*y01+e*aa*y10+aa*bb*y11)/(1-b*e*e)
    return ((e-e0)**2*y00+(e-e0)*ct*y01+cs*(e-e0)*y10+cs*ct*y11)/(2*e0)
def original_objective(e,nodes,weights,b):
    n=len(nodes)+1;zero=Decimal(0);A=[[zero for _ in range(n)] for _ in range(n)];A[0][0]=e
    for j,(s,w) in enumerate(zip(nodes,weights)):A[0][j+1]=az(e,s,1)*w;A[j+1][j+1]=s
    Y=[[zero for _ in range(n)] for _ in range(n)]
    for i in range(n-1,-1,-1):
        for j in range(n-1,-1,-1):
            rhs=Decimal(i>0 and j>0)
            for k in range(i,n):
                for l in range(j,n):
                    if (k,l)!=(i,j):rhs+=b*A[i][k]*A[j][l]*Y[k][l]
            Y[i][j]=rhs/(1-b*A[i][i]*A[j][j])
    track=[e]+[(az(e,s,1)-1)*w for s,w in zip(nodes,weights)];adjust=[e-1]+[az(e,s,1)*w for s,w in zip(nodes,weights)]
    return sum(((b/15)*track[i]*track[j]+adjust[i]*adjust[j])*Y[i][j] for i in range(n) for j in range(n))/2
def main():
    raw=(HERE/'standard_sqp_comparison.json').read_bytes();d=json.loads(raw);pr=(HERE/'sqp_comparison_proposals.json').read_bytes();proposal=json.loads(pr);assert hashlib.sha256(pr).hexdigest()==d['proposal_sha256'];sources={}
    for key,z in d['sources'].items():
        p=HERE.parent/z['path'];rr=p.read_bytes();assert hashlib.sha256(rr).hexdigest()==z['sha256'];sources[key]=json.loads(rr)
        if key in ['B26','B28']:
            ck=json.loads((p.parent/'artifact_check.json').read_text(encoding='utf-8'));assert ck['source_sha256']==z['sha256']
    count=0;best=Decimal(0);bestindex=None
    with localcontext() as ctx:
        ctx.prec=180;cal=sources['reports'];b=dec('9999/10000');e0=root(b,Decimal(1));upper=dec(sources['B28']['absolute_cost_upper']);old=dec(sources['B26']['genuine_compatible_loss_lower'])
        for w in d['certified_worlds']:
            i=w['index'];assert w['method_run']==proposal['runs'][i];p=w['parameters'];beta=dec(p['beta']);r=dec(p['report_pole']);g=dec(p['gain']);nodes=list(map(dec,w['nodes']));weights=list(map(dec,w['nu_weights']));initials=list(map(dec,w['initials']))
            assert 0<beta<1 and 0<r<1 and g>0 and all(0<=s<1 for s in nodes) and min(weights+initials)>=0 and len(nodes)<=5
            assert all(F(w['parameters'][k])==F(str(proposal['candidates'][i][kk])) for k,kk in [('beta','beta'),('report_pole','r'),('gain','g')])
            assert list(map(F,w['nodes']))==[F(str(x)) for x in proposal['candidates'][i]['nodes']]
            assert list(map(F,w['nu_weights']))==[F(str(x)) for x in proposal['candidates'][i]['weights']]
            close(sum(ww/(1-s) for s,ww in zip(nodes,weights)),dec(w['finite_mu_mass']));pred=[]
            for ch,k in enumerate([Decimal(1),Decimal(2)]):
                e=root(beta,k);inv,rep=initials[2*ch:2*ch+2]
                for n in range(4):inv=e*inv+sum(ww*az(e,s,k)*(s**n if n else Decimal(1)) for s,ww in zip(nodes,weights));rep=r*rep+g*inv;pred.append(rep)
            inv,rep=initials[4:]
            for n in range(4):rep=r*rep+g*(inv+(dec(cal['known_physical_step']) if n>=2 else 0));pred.append(rep)
            for j,(v,z) in enumerate(zip(pred,w['report_checks'])):
                center=dec(cal['main_centers'][j//4][j%4]) if j<8 else dec(cal['calibration_centers'][j-8]);rad=dec(cal['main_raw_radius'] if j<8 else cal['calibration_raw_radius']);close(center,dec(z['center']));close(rad,dec(z['single_radius']));inside(v,z['prediction_outer']);observed=dec(z['concrete_observation']);close(observed-center,dec(z['actual_error']));inside(v-observed,z['candidate_error_outer']);assert abs(observed-center)<rad and abs(v-observed)<rad;count+=1
            e=root(beta,Decimal(1));energy=sum(ww*vv*pair_energy(e,e0,s,t,b) for s,ww in zip(nodes,weights) for t,vv in zip(nodes,weights));gap=original_objective(e,nodes,weights,b)-original_objective(e0,nodes,weights,b);close(energy,gap);inside(energy,w['residual_energy_outer']);inside(gap,w['original_objective_gap_outer']);assert dec(w['certified_loss_lower'])<=energy<upper
            close(upper/dec(w['certified_loss_lower']),dec(w['unchanged_global_upper_to_lower_ratio']));assert abs(energy-dec(str(proposal['candidates'][i]['value'])))<Decimal('1e-19')
            if dec(w['certified_loss_lower'])>best:best=dec(w['certified_loss_lower']);bestindex=i
        assert bestindex==d['best_world_index'] and best==dec(d['best_loss_lower']);close(best/old,dec(d['improvement_factor']));close(upper/best,dec(d['best_upper_to_lower_ratio']));close(upper,dec(d['unchanged_B28_global_upper']));assert best>old and upper/best<dec('1033/1000')
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_decimal_digits=180,standard_optimizer_runs=5,
        rationally_qualified_worlds=5,independent_original_report_noise_splits=count,independent_pair_Stein_energy_evaluations=5,independent_original_objective_gap_evaluations=5,
        optimizer_success_flags=sum(z['success'] for z in proposal['runs']),same_original_twelve_report_union_preserved=True,
        complex_step_proposal_candidates_independently_qualified=True,local_node_cap_imposed_on_global_upper=False,
        independent_SDP_optimization_completed=False,full_unknown_parameter_global_optimum_claimed=False,
        standard_solver_success_flags_overridden=False,algorithm_superiority_or_equal_resource_comparison_claimed=False,
        representation_splits_counted_as_independent_physical_starts=False,B28_global_upper_changed=False,
        best_loss_lower=d['best_loss_lower'],best_upper_to_lower_ratio=d['best_upper_to_lower_ratio'])
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in out.items() if k not in ['best_loss_lower','best_upper_to_lower_ratio']}))
if __name__=='__main__':main()
