from pathlib import Path
import financing_case as m
import json,math,numpy as np
from scipy.stats import beta
OUT=m.OUT
data=json.loads((OUT/'case_results.json').read_text('utf-8'));rng=np.random.default_rng(72026013);raw=[];n=200000
for rep in [False,True]:
    plan=next(p for p in data['plans'] if p['rho']==1 and p['repeat']==rep);q=plan['q'];ker,_,_=m.kernel(q,rep)
    for theta in [-1,1]:
        w=m.W0+theta*m.D;types=rng.choice(len(m.Z),size=n,p=w/w.sum());lam=q*m.Z[types];first=rng.exponential(1/lam);hit=first<=1;bins=np.minimum((first[hit]*8).astype(int),7);category=np.zeros(n,int)
        if rep:
            extra=rng.poisson(lam[hit]*(1-first[hit]));category[hit]=1+3*bins+np.minimum(extra,2)
        else:category[hit]=1+bins
        counts=np.bincount(category,minlength=ker.shape[1]);truth=w@ker;alpha=.001/(2*(9+25));miss=0
        for num,p in zip(counts,truth):
            lo=0 if num==0 else float(beta.ppf(alpha/2,num,n-num+1));hi=1 if num==n else float(beta.ppf(1-alpha/2,num+1,n-num));miss+=int(p<lo or p>hi)
        raw.append({'repeat':rep,'theta':theta,'identities':n,'categories':len(counts),'simultaneous_interval_misses':miss,'max_frequency_error':float(max(abs(counts/n-truth)))})
assert all(r['simultaneous_interval_misses']==0 for r in raw)
(OUT/'独立消息生成器核验.json').write_text(json.dumps(raw,ensure_ascii=False,indent=2),'utf-8')
normalization=[]
for p in data['plans']:
    model=m.model(p['rho']);design=m.design(model,p['q'],p['repeat']);I=design['score_information_per_identity']
    for theta in [-1.,-.5,0.,.5,1.]:
        prob=(m.W0+theta*p['rho']*m.D)@design['matrix'];err=abs(float(prob@design['score'])-theta*I);assert err<1e-12
    ncap=design['cap'];lg=math.log(2*len(m.CHECKS)/m.STAT_ERROR)
    radius=(math.sqrt(2*design['variance_max']*lg/(m.K*ncap))+2*design['range_bound']*lg/(3*m.K*ncap))/I
    assert abs(radius-design['radius'][ncap])<1e-12
    normalization.append({'rho':p['rho'],'repeat':p['repeat'],'radius_at_cap':radius,'includes_inverse_information_slope':True})
(OUT/'估计器归一化核验.json').write_text(json.dumps(normalization,ensure_ascii=False,indent=2),'utf-8')
print(json.dumps({'raw_message_identities':sum(r['identities'] for r in raw),'normalization_checks':len(normalization)}))
