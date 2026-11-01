"""B40 candidate-only Clarabel matrix-SOS support synthesis.
Same 8 report coefficients and original calibration/error budget as B35.
PSD constraints cover persistence continuously at three beta proposal nodes.
Full continuous beta/persistence qualification is a separate mandatory step.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib,math
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
sys.path.insert(0,str(BASE/'运行依赖_SDP'));sys.path.insert(0,str(HERE.parent/'第三十五轮研究_20261027'))
import standard_lp_support_cover as lp
import numpy as np,clarabel
from scipy import sparse
base=lp.base;support=lp.support;pp=lp.proposal

def add(a,b):return np.polynomial.polynomial.polyadd(a,b)
def mul(*args):
    z=np.array([1.])
    for a in args:z=np.polynomial.polynomial.polymul(z,a)
    return z
def pad(a,n=9):
    assert len(a)<=n
    return np.pad(a,(0,n-len(a)))
def features(d,reduced=False):
    if not reduced:return [(r,np.r_[np.zeros(k),1.]) for k in range(d+1) for r in range(4)]
    return [(r,np.r_[np.zeros(k),1.]) for k in range(d+1) for r in range(3)]+[(3,np.r_[np.zeros(k),1.,-1.]) for k in range(d)]
def pairs(n):return [(i,j) for j in range(n) for i in range(j+1)]

def arrow(beta,direction,T):
    a,b=map(float,direction);rho=1/15;ref=.9999;e0=pp.w.t.root(ref,rho,1);e,e2=[pp.w.t.root(beta,rho,k) for k in [1,2]]
    D1=np.array([rho+1-e,-(1-e)]);D2=np.array([rho+2*(1-e2),-2*(1-e2)]);D0=np.array([rho+1-e0,-(1-e0)]);EE=np.array([1,-ref*e])
    C=(1-e)*(1-e0)/rho;R=1-e+C;V=1-ref*e*e;TT=V+ref*R*R;BN=V*(e+R)+R*C;NN=np.array([TT,-ref*BN]);S=math.sqrt(2*e0*V*TT)
    basis=[]
    for group,z in enumerate([e,e2]):
        for n in range(4):
            pol=np.array([z**(n-j) for j in range(n+1)])
            basis.append(pad(S*rho*(1-z)*mul(D0,EE,D2 if group==0 else D1,pol)))
    XV=np.array(basis).T@T;XC=pad(-a*rho*rho*mul(D2,NN));t=np.array([1,0,-ref]);u=np.array([1,-1])
    mv=np.zeros((9,4,4,8));mc=np.zeros((9,4,4))
    for k in range(8):
        x=XV[:,k]
        # Trim structural zero high powers before polynomial multiplication.
        x=np.trim_zeros(x,'b')
        mv[:,0,0,k]=pad(mul(t,x));mv[:,1,1,k]=pad(x);mv[:,2,2,k]=pad(x);mv[:,3,3,k]=pad(mul(u,x))
    xc=np.trim_zeros(XC,'b');mc[:,0,0]=pad(mul(t,xc));mc[:,1,1]=pad(xc);mc[:,2,2]=pad(xc);mc[:,3,3]=pad(mul(u,xc))
    W=e*(2-e*(1+ref))+R*R;k0=(1-ref*e)**2*(1-ref);k1=2*ref*(1-ref*e)**2;k2=ref*ref*W
    fac=b*rho*rho*(1+C)*math.sqrt(ref*V);sD=mul([0,1],D2)
    for j,z in [(1,fac*math.sqrt(k0)*sD),(2,fac*math.sqrt(k2)*mul(u,sD)),(3,fac*math.sqrt(k1)*mul(u,sD))]:mc[:,0,j]=pad(z);mc[:,j,0]=pad(z)
    return mv,mc,dict(beta=beta,e=e,e2=e2,e0=e0,arrow_scale=float(np.max(np.abs(mv))),endpoint_u1_positive_terms=[k0,k1,k2],target_direction=[a,b])

def pole_nodes(band,cal,pieces,step,strips,T):
    nodes={}
    for i in band['allowed_pole_indices']:
        ll,rr=map(F,strips[i]['pole_interval'])
        for pp0 in pieces:
            l=max(ll,pp0['interval'][0]);r=min(rr,pp0['interval'][1])
            if l>=r:continue
            for j in range(17):
                z=l+(r-l)*F(j,16);gain=np.polynomial.polynomial.polyval(float(z),list(map(float,pp0['p'])))/(float(step)*np.polynomial.polynomial.polyval(float(z),list(map(float,pp0['den']))))
                key=(z,float(gain));nodes[key]=(float(z),float(gain))
    return list(nodes.values())

def program(band,cal,pieces,step,strips,direction,kind,scalar_guard=True):
    prepared=lp.prepare(band,cal,pieces,step,strips);Agrid,T,oldmatrix,coords=prepared
    poles=pole_nodes(band,cal,pieces,step,strips,T);nbase=9+8*len(poles);nvars=nbase
    betas=[float(F(band['beta_interval'][0])),float(sum(map(F,band['beta_interval']),F(0))/2),float(F(band['beta_interval'][1]))]
    qfeat=features(4,True);rfeat=features(3);graminfo=[]
    if kind=='SDP':
        for beta in betas:
            for feat in [qfeat,rfeat]:
                graminfo.append((nvars,feat));nvars+=len(pairs(len(feat)))
    eq=[];erhs=[];ineq=[];irhs=[];meta=[]
    def eqrow(z,rhs):eq.append(z);erhs.append(rhs)
    def le(z,rhs=0.):ineq.append(z);irhs.append(rhs)
    # Same epigraph budget for both methods. No 256 raw-error vertices needed.
    y=np.array([float(F(x)) for row in cal['main_centers'] for x in row]);rad=2*float(F(cal['main_raw_radius']))
    for j,(pole,gain) in enumerate(poles):
        E=np.eye(8)
        for k in [0,1,2,4,5,6]:E[k,k+1]-=pole
        KT=E@T;start=9+8*j
        for k in range(8):
            z={i:KT[k,i] for i in range(8)};z[start+k]=-1;le(z)
            z={i:-KT[k,i] for i in range(8)};z[start+k]=-1;le(z)
        z={i:(y@KT)[i] for i in range(8)};z[8]=-gain
        for k in range(8):z[start+k]=rad
        le(z)
    le({8:-1.})
    for beta in betas:
        es=[pp.w.t.root(beta,1/15,k) for k in [1,2]]
        for group,e in enumerate(es):
            v=np.zeros(8);v[group*4]=1
            le({i:-(v@T)[i] for i in range(8)})
            v[group*4:group*4+4]=[e**(k+1) for k in range(4)]
            le({i:-(v@T)[i] for i in range(8)})
    # Additional scalar samples keep conic residuals from exploiting small
    # endpoint eigenvalues. Their zero-buffer version follows from SOS;
    # the 1e-8 buffer is a shared tightening for both proposal methods.
    a,b=map(float,direction);target=np.concatenate([np.r_[a*q+b*h,np.zeros(4)] for q,h in coords]);AT=Agrid@T
    if kind=='LP' or scalar_guard:
        for row,rhs in zip(AT,target):le({i:-row[i] for i in range(8)},-rhs-1e-8)
    if kind=='SDP':
        for ib,beta in enumerate(betas):
            mv,mc,info=arrow(beta,direction,T);scale=max(info['arrow_scale'],1e-14);meta.append(info)
            # Coefficient identity after exact endpoint facial reduction:
            # first Gram feature for component 3 has factor (1-s).
            mapping={}
            for gg,(start,feat) in enumerate(graminfo[2*ib:2*ib+2]):
                weight=[1.] if gg==0 else [0.,1.,-1.]
                for v,(i,j) in enumerate(pairs(len(feat))):
                    ri,fi=feat[i];rj,fj=feat[j];r,c0=sorted([ri,rj]);pol=mul(fi,fj,weight)*(2 if ri==rj and i!=j else 1)
                    for degree,coef in enumerate(pol):
                        if coef:mapping.setdefault((degree,r,c0),{})[start+v]=-scale*coef
            for degree in range(9):
                for r in range(4):
                    for c0 in range(r,4):
                        if c0==3 and degree==8:continue
                        z={i:mv[degree,r,c0,i] for i in range(8) if mv[degree,r,c0,i]!=0};z.update(mapping.get((degree,r,c0),{}));eqrow(z,-mc[degree,r,c0])
    def matrix(rows,rhs,scaled=True):
        ii=[];jj=[];vv=[];bb=[]
        for i,(z,b0) in enumerate(zip(rows,rhs)):
            sc=max(max(map(abs,z.values()),default=0.),abs(b0),1e-12) if scaled else 1.
            for j,v in z.items():
                if v:ii.append(i);jj.append(j);vv.append(v/sc)
            bb.append(b0/sc)
        return sparse.csc_matrix((vv,(ii,jj)),shape=(len(rows),nvars)),np.array(bb)
    em,eb=matrix(eq,erhs);im,ib=matrix(ineq,irhs);blocks=[em,im];rhs=[eb,ib];cones=[]
    if len(eq):cones.append(clarabel.ZeroConeT(len(eq)))
    cones.append(clarabel.NonnegativeConeT(len(ineq)))
    for start,feat in graminfo:
        nn=len(feat);pq=pairs(nn);mat=sparse.csc_matrix(([-1 if i==j else -2**.5 for i,j in pq],(range(len(pq)),[start+k for k in range(len(pq))])),shape=(len(pq),nvars))
        blocks.append(mat);rhs.append(np.zeros(len(pq)));cones.append(clarabel.PSDTriangleConeT(nn))
    mat=sparse.vstack(blocks,format='csc');bb=np.concatenate(rhs);obj=np.zeros(nvars);obj[8]=1
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.tol_gap_abs=1e-9;settings.tol_gap_rel=1e-9;settings.tol_feas=1e-9;settings.max_iter=200;settings.time_limit=90.;settings.max_threads=1
    begin=time.monotonic();sol=clarabel.DefaultSolver(sparse.csc_matrix((nvars,nvars)),obj,mat,bb,cones,settings).solve();elapsed=time.monotonic()-begin
    out=dict(method=kind,scalar_support_samples_included=(kind=='LP' or scalar_guard),solver='Clarabel '+clarabel.__version__,status=str(sol.status),objective=float(sol.obj_val),dual_objective=float(sol.obj_val_dual),iterations=sol.iterations,seconds=elapsed,primal_residual=float(sol.r_prim),dual_residual=float(sol.r_dual),variables=nvars,equalities=len(eq),nonnegative_rows=len(ineq),pole_samples=len(poles),PSD_orders=[len(x[1]) for x in graminfo],beta_samples=betas,matrix_diagnostics=meta,qualification_not_yet_run=True)
    if sol.x is not None and np.isfinite(sol.x).all():
        c=T@np.array(sol.x[:8]);out['candidate_coefficients']=[str(F(round(float(z)*10**12),10**12)) for z in c];out['floating_coefficients']=c.tolist();out['transformed_weights']=list(sol.x[:8]);out['max_scaled_constraint_residual']=float(np.max(np.abs(mat@sol.x+np.array(sol.s)-bb)));out['minimum_original_scalar_sample_slack']=float(np.min(Agrid@c-target))
        out['grams']=[]
        for start,feat in graminfo:
            Q=np.zeros((len(feat),len(feat)))
            for k,(i,j) in enumerate(pairs(len(feat))):Q[i,j]=Q[j,i]=sol.x[start+k]
            out['grams'].append(dict(order=len(feat),smallest_eigenvalue=float(np.linalg.eigvalsh(Q)[0])))
    return out

def main():
    roots={};data={}
    for key,name in [('B39','第三十九轮研究_20261027/adaptive_direction_cover.json'),('B28','第二十八轮研究_20261026/polygon_performance_certificate.json'),('cal','第七轮研究_20261023/physical_step_absolute_audit.json'),('excl','第二十七轮研究_20261026/joint_compatibility_exclusion.json')]:
        path=HERE.parent/name;raw=path.read_bytes();roots[key]=dict(path=name,sha256=hashlib.sha256(raw).hexdigest());data[key]=json.loads(raw)
    b39=data['B39'];entry=max(b39['completed_subbands'],key=lambda z:F(z['hybrid_local_loss_upper']));part=json.loads((HERE.parent/'第三十九轮研究_20261027'/entry['path']).read_text(encoding='utf-8'))
    parent=entry['source_parent_index'];child=entry['source_child_index'];band=next(z for z in data['B28']['beta_records'] if z['source_beta_certificate_index']==parent)['subbands'][child]
    cal=data['cal'];strips=data['excl']['strips'];_,_,step,pieces=base.calmod.frontier(cal);direction=list(map(F,part['runs'][0]['direction']))
    out=dict(source_artifacts=roots,construction_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parent=parent,child=child,direction=list(map(str,direction)),proposal_scope='Local controlled-direction support proposals only. Same original 8 coefficients/calibration/error set. LP and SDP share compressed raw-budget epigraph, scalar support samples and beta/pole proposal samples; SDP adds continuous-persistence matrix SOS. Direction/partition inherited. No continuous beta qualification, global bound, full 51-variable SDP, equal-resource superiority or exact optimum claimed.',runs=[])
    for kind,with_grid in [('LP',True),('SDP',False),('SDP',True)]:
        z=program(band,cal,pieces,step,strips,direction,kind,with_grid);out['runs'].append(z)
        (HERE/'sdp_support_pilot.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in z.items() if k not in ['candidate_coefficients','floating_coefficients','transformed_weights','grams','matrix_diagnostics']}),flush=True)
if __name__=='__main__':main()
