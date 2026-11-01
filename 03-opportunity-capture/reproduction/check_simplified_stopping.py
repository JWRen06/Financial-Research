"""Numerical cross-check of the simplified deterministic stopping bound in Appendix C."""
from pathlib import Path
import sys,itertools,json,math
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[2]
for name in ['_local_dependencies','_plot_dependencies']:
    p=BASE/'工作区/04_第二第三篇共同研究/原始工作目录/experiments'/name
    if p.exists():sys.path.insert(0,str(p))
import numpy as np
def params(p):
    a=p[0];b=2*p[0]-p[1];d=3*p[0]-3*p[1]+p[2]
    if not a>b>d>0 or a*d<=b*b:raise ValueError()
    r=a-b*b/d;s=d/b;A=b**3/d**2;x=(b-d)/(a-b);w=(a-b)/(1-x)
    return r,A,s,a-w,w/x,x
def curves(p,q):
    r,A,s,m,M,x=params(p)
    return r*q+A*(-np.expm1(q*np.log1p(-s))),m+M*(-np.expm1(q*np.log1p(-x)))
def optimum(p,c,Q):
    r,A,s,m,M,x=params(p);zs=-np.log1p(-s);zx=-np.log1p(-x)
    opt=[(0.,0.)]
    for lo,hi,index in [(0,1,0),(1,2,1),(2,3,0),(3,Q,1)]:
        if lo>Q:continue
        hi=min(hi,Q)
        if index==0:q=hi if c<=r else float(np.clip(np.log(A*zs/(c-r))/zs,lo,hi))
        else:q=float(np.clip(np.log(M*zx/c)/zx,lo,hi))
        opt.append((q,curves(p,q)[index]-c*q))
    return max(opt,key=lambda x:x[1])
def minorant(p,q):
    r,A,s,m,M,x=params(p);right=1<=q<=2 or q>=3;x=x if right else s
    g=-np.expm1(q*np.log1p(-x))/x
    gp=(q*x*np.exp((q-1)*np.log1p(-x))-1+(1-x)**q)/x**2
    if right:gamma=(1-g-gp*(1-x))/(1-x)**2;beta=gp-2*gamma*x;alpha=g-gp*x+gamma*x*x
    else:gamma=(x*gp-g+q)/x**2;beta=2*(g-q)/x-gp;alpha=q
    return np.array([alpha+2*beta+3*gamma,-beta-3*gamma,gamma])
rng=np.random.default_rng(20261031);checked=0;largest_violation=-1e100
for case in range(120):
    v=rng.uniform(.02,.95,7);w=rng.uniform(.1,4,7)
    p=np.array([sum(w*(-np.expm1(k*np.log1p(-v)))) for k in [1,2,3]])
    e=1e-5*max(p);c=float(rng.uniform(.03,.8)*p[0]);Q=float(rng.choice([4,10,50,500]))
    vertices=np.array(list(itertools.product(*zip(p-e,p+e))))
    for vp in vertices:params(vp)
    q,_=optimum(p,c,Q)
    if q>0:
        coef=minorant(p,q);W=float(coef@p-e*sum(abs(coef))-c*q)
    else:W=0.
    if W<0:q=0.;W=0.
    Z=max(optimum(vp,c,Q)[1] for vp in vertices)
    for _ in range(40):
        truth=p+rng.uniform(-e,e,3)
        actual=0 if q==0 else min(curves(truth,q))-c*q
        best=optimum(truth,c,Q)[1]
        violation=max(W-actual,best-Z,best-actual-(Z-W))
        largest_violation=max(largest_violation,violation)
        assert violation<2e-8,(case,q,W,Z,actual,best)
        checked+=1
result={'boxes':120,'interior_points':checked,'largest_violation':largest_violation,'passed':True,'scope':'Implementation cross-check, not a replacement for the convexity and interpolation proof.'}
(HERE/'简化停止规则核查.json').write_text(json.dumps(result,indent=2),'utf-8');print(result)
