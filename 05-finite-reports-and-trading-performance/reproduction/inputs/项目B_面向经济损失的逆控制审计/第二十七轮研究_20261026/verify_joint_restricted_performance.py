"""Independent 110-digit rational proof of affine beta exclusions and budgets.
Reuses the unchanged, previously independently verified B20 target envelopes
by hash; reconstructs every new exclusion and every restricted raw budget.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import json,hashlib,time
HERE=Path(__file__).resolve().parent;DEN=10**110;RHO=F(1,15)
def floor(v):return F((v*DEN).__floor__(),DEN)
def ceil(v):return F((v*DEN).__ceil__(),DEN)
class Q:
    def __init__(self,l,r=None):self.l=floor(F(l));self.r=ceil(F(l if r is None else r));assert self.l<=self.r
    def __add__(self,b):
        b=asq(b);return Q(self.l+b.l,self.r+b.r)
    __radd__=__add__
    def __neg__(self):return Q(-self.r,-self.l)
    def __sub__(self,b):return self+-asq(b)
    def __rsub__(self,b):return asq(b)+-self
    def __mul__(self,b):
        b=asq(b);v=[self.l*b.l,self.l*b.r,self.r*b.l,self.r*b.r];return Q(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=asq(b);assert b.l*b.r>0;return self*Q(1/b.r,1/b.l)
def asq(v):return v if isinstance(v,Q) else Q(v)
def power(v,n):
    answer=Q(1)
    for _ in range(n):answer=answer*v
    return answer
@lru_cache(None)
def root(beta,k):
    left=F(0);right=F(1)
    for _ in range(360):
        mid=(left+right)/2;sign=k*(1-mid)*(1-beta*mid)-RHO*beta*mid
        if sign>0:left=mid
        else:right=mid
    assert k*(1-left)*(1-beta*left)-RHO*beta*left>0
    assert k*(1-right)*(1-beta*right)-RHO*beta*right<0
    return Q(left,right)
def multiply(a,b):
    c=[Q(0) for _ in range(len(a)+len(b)-1)]
    for j,x in enumerate(a):
        for k,y in enumerate(b):c[j+k]=c[j+k]+x*y
    return c
def add(*args):
    return [sum((x[j] if j<len(x) else Q(0) for x in args),Q(0)) for j in range(max(map(len,args)))]
def rebuild(band,leftcoef,rightcoef):
    bl,br=band;bm=(bl+br)/2;midcoef=[(a+b)/2 for a,b in zip(leftcoef,rightcoef)]
    slope=[(b-a)/(br-bl) for a,b in zip(leftcoef,rightcoef)]
    cb=[Q(min(a,b),max(a,b)) for a,b in zip(leftcoef,rightcoef)];den=[];dend=[];num=[];numd=[];centerden=[];centernum=[];ini=[];inid=[]
    for fee,k in enumerate([F(1),F(2)]):
        em=root(bm,k);e=Q(root(br,k).l,root(bl,k).r);beta=Q(bl,br)
        ep=-e*(k*(1-e)+RHO)/(k*(1+beta-2*beta*e)+RHO*beta);assert ep.r<0
        den.append([RHO+k*(1-e),-k*(1-e)]);dend.append([-k*ep,k*ep]);centerden.append([RHO+k*(1-em),-k*(1-em)])
        nc=[];nd=[];nm=[]
        for j in range(4):
            base=sum((cb[4*fee+n]*power(e,n-j) for n in range(j,4)),Q(0))
            deriv=sum((slope[4*fee+n]*power(e,n-j)+(cb[4*fee+n]*(n-j)*power(e,n-j-1)*ep if n>j else Q(0)) for n in range(j,4)),Q(0))
            nc.append(RHO*(1-e)*base);nd.append(RHO*(-ep*base+(1-e)*deriv))
            nm.append(RHO*(1-em)*sum((midcoef[4*fee+n]*power(em,n-j) for n in range(j,4)),Q(0)))
        num.append(nc);numd.append(nd);centernum.append(nm)
        ini.append(sum((midcoef[4*fee+n]*power(em,n+1) for n in range(4)),Q(0)))
        inid.append(sum((slope[4*fee+n]*power(e,n+1)+cb[4*fee+n]*(n+1)*power(e,n)*ep for n in range(4)),Q(0)))
    center=add(multiply(centernum[0],centerden[1]),multiply(centernum[1],centerden[0]))
    deriv=add(multiply(numd[0],den[1]),multiply(num[0],dend[1]),multiply(numd[1],den[0]),multiply(num[1],dend[0]))
    radius=(br-bl)/2;assert min((a+b*Q(-radius,radius)).l for a,b in zip(ini,inid))>=0
    return center,deriv
def bernstein(poly,left,right):
    n=len(poly)-1
    shifted=[sum((poly[k]*comb(k,j)*left**(k-j)*(right-left)**j for k in range(j,n+1)),Q(0)) for j in range(n+1)]
    return [sum((shifted[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),Q(0)) for i in range(n+1)]
def quadratic(poly,left,right):
    a,b,c=poly;points=[left,right]
    if c>0 and left<-b/(2*c)<right:points.append(-b/(2*c))
    return min(a+b*x+c*x*x for x in points)
def main():
    started=time.monotonic();raw=(HERE/'joint_restricted_performance.json').read_bytes();d=json.loads(raw)
    eraw=(HERE/'joint_compatibility_exclusion.json').read_bytes();exc=json.loads(eraw)
    oraw=(HERE.parent/'第二十轮研究_20261025/global_directional_certificate.json').read_bytes();old=json.loads(oraw)
    ocraw=(HERE.parent/'第二十轮研究_20261025/artifact_check.json').read_bytes();oc=json.loads(ocraw)
    lraw=(HERE.parent/'第二十六轮研究_20261026/continuous_five_node_witness.json').read_bytes();low=json.loads(lraw)
    lcraw=(HERE.parent/'第二十六轮研究_20261026/artifact_check.json').read_bytes();lc=json.loads(lcraw)
    calraw=(HERE.parent/'第七轮研究_20261023/physical_step_absolute_audit.json').read_bytes();cal=json.loads(calraw)
    for data,key in [(eraw,'beta_pole_exclusion_sha256'),(oraw,'old_directional_source_sha256'),(ocraw,'old_directional_check_sha256'),(lraw,'lower_source_sha256'),(calraw,'calibration_source_sha256')]:assert hashlib.sha256(data).hexdigest()==d[key]
    assert hashlib.sha256(oraw).hexdigest()==oc['source_sha256'] and hashlib.sha256(lraw).hexdigest()==lc['source_sha256']
    assert hashlib.sha256(calraw).hexdigest()==exc['source_calibration_sha256']
    cy=list(map(F,cal['calibration_centers']));ce=2*F(cal['calibration_raw_radius']);l0,l1,l2,l3=[v-ce for v in cy];u0,u1,u2,u3=[v+ce for v in cy]
    rlo=(l3-u2)/(u2-l1);rhi=(u3-l2)/(l2-u1);assert exc['pole_outer']==list(map(str,[rlo,rhi]))
    y=[F(x) for row in cal['main_centers'] for x in row];eps=2*F(cal['main_raw_radius']);strips=exc['strips'];lastpole=rlo;cells=rectangles=0
    for si,strip in enumerate(strips):
        rl,rr=map(F,strip['pole_interval']);assert rl==lastpole<rr;lastpole=rr;lastbeta=F('.978')
        for row in strip['certificates']:
            bl,br=map(F,row['beta_interval']);assert bl==lastbeta<br;lastbeta=br
            leftcoef,rightcoef=[list(map(F,v)) for v in row['holding_coefficients_at_beta_endpoints']]
            assert min(leftcoef[0],leftcoef[4],rightcoef[0],rightcoef[4])>=0
            for coef,saved in zip([leftcoef,rightcoef],row['raw_report_checks_at_beta_endpoints']):
                for r,stored in zip([rl,rr],saved):
                    rawcoef=[coef[j]-r*coef[j+1] if j%4<3 else coef[j] for j in range(8)]
                    center=sum((a*b for a,b in zip(rawcoef,y)),F(0));error=eps*sum(map(abs,rawcoef))
                    assert rawcoef==list(map(F,stored['raw_coefficients'])) and center==F(stored['center']) and error==F(stored['error'])
                    assert center+error==F(stored['upper'])<0
            center,deriv=rebuild((bl,br),leftcoef,rightcoef);radius=(br-bl)/2;lasts=F(0)
            for cover in row['whole_support_cover']:
                sl,sr=map(F,cover['s_interval']);assert sl==lasts<sr;lasts=sr
                bc=bernstein(center,sl,sr);bd=bernstein(deriv,sl,sr)
                assert min((a+b*Q(-radius,radius)).l for a,b in zip(bc,bd))>=0
                cells+=1
            assert lasts==1;rectangles+=1
        assert lastbeta==F(strip['excluded_beta_interval'][1]) and F(strip['excluded_beta_interval'][0])==F('.978')
        print(json.dumps(dict(exclusion_poles_checked=si+1,rectangles=rectangles,elapsed=round(time.monotonic()-started,1))),flush=True)
    assert lastpole==rhi and F(exc['new_uniform_beta_lower'])==min(F(s['excluded_beta_interval'][1]) for s in strips)
    # Reconstruct the original physical-calibration lower-gain frontier.
    r0=(l3-u2)/(u2-u1);r1=(l3-l2)/(l2-u1);step=F(cal['known_physical_step'])
    pieces=[([rlo,r0],[l3-u2,l3-u2,l0-u2],[F(0),F(1)]),([r0,r1],[l3-u1,l0-u1,l0-u1],[F(1),F(1)]),([r1,rhi],[l2-u1,l0-u1,F(0)],[F(1),F(0)])]
    for actual,(interval,p,den) in zip(d['gain_frontier'],pieces):assert actual==dict(pole_interval=list(map(str,interval)),numerator=list(map(str,p)),denominator=list(map(str,den)))
    rawcells=0;newbudgets=0
    for entry in d['registry'].values():
        original=old['registry'][entry['source_budget_id']];coef=list(map(F,original['normalized_coefficients']));M=F(entry['raw_upper']);assert 0<=M<=F(original['normalized_raw_upper'])
        allowed=entry['allowed_pole_indices']
        if entry['inherited_full_pole_budget']:
            assert allowed==list(range(len(strips))) and M==F(original['normalized_raw_upper']);continue
        newbudgets+=1
        for index in allowed:
            last=F(strips[index]['pole_interval'][0]);target=F(strips[index]['pole_interval'][1])
            for cell in [v for v in entry['raw_proof'] if v['pole_index']==index]:
                left,right=map(F,cell['pole_interval']);assert left==last<right;last=right
                interval,p,den=pieces[cell['gain_piece']];assert interval[0]<=left<right<=interval[1]
                forms=[(coef[j],-coef[j+1] if j%4<3 else F(0)) for j in range(8)];A=B=F(0)
                for (a,b),sign,observed in zip(forms,cell['raw_signs'],y):
                    assert sign in [-1,1] and sign*(a+b*left)>=0 and sign*(a+b*right)>=0
                    A+=a*(observed+eps*sign);B+=b*(observed+eps*sign)
                assert [A,B]==list(map(F,cell['raw_numerator']))
                d0,d1=den;prod=[d0*A,d0*B+d1*A,d1*B];poly=[M*x-step*z for x,z in zip(p,prod)]
                assert poly==list(map(F,cell['quadratic'])) and quadratic(poly,left,right)==F(cell['minimum'])>=0
                assert quadratic(p,left,right)>0;rawcells+=1
            assert last==target
    assert len(d['beta_records'])==len(old['certificates']);uppers=[];retained=0
    for i,(row,previous) in enumerate(zip(d['beta_records'],old['certificates'])):
        assert row['source_beta_certificate_index']==i and row['beta_interval']==previous['beta_interval']
        right=F(row['beta_interval'][1]);allowed=[j for j,s in enumerate(strips) if right>F(s['excluded_beta_interval'][1])];assert allowed==row['allowed_pole_indices']
        if not allowed:continue
        retained+=1;M={}
        for kind,key in row['components'].items():
            entry=d['registry'][key];assert entry['allowed_pole_indices']==allowed and entry['source_budget_id']==previous['components'][kind]['budget_id'];M[kind]=F(entry['raw_upper'])
        upper=F(previous['root_gap_outer'][1])**2*(max(M['plus']**2,M['minus']**2)+M['perpendicular']**2)
        assert upper==F(row['absolute_cost_upper'])<=F(previous['absolute_cost_upper']);uppers.append(upper)
    upper=max(uppers);lower=F(low['genuine_compatible_loss_lower']);prior=F(old['absolute_cost_upper'])
    assert lower<upper<prior and upper==F(d['absolute_cost_upper']) and prior==F(d['previous_upper'])
    assert lower==F(d['genuine_compatible_lower']) and upper/lower==F(d['union_gap_ratio'])
    assert F(d['nominal_compatible_lower'])==F(old['nominal_compatible_lower']) and upper/F(d['nominal_compatible_lower'])==F(d['nominal_gap_ratio'])
    assert 1-upper/prior==F(d['relative_upper_reduction'])
    result=dict(source_sha256=hashlib.sha256(raw).hexdigest(),exclusion_sha256=hashlib.sha256(eraw).hexdigest(),
        lower_check_sha256=hashlib.sha256(lcraw).hexdigest(),independent_affine_beta_pole_rectangles=rectangles,
        independent_exclusion_persistence_cells=cells,independent_precision_digits=110,
        independent_restricted_raw_budgets=newbudgets,independent_raw_quadratic_cells=rawcells,
        old_beta_envelopes_accounted_for=len(d['beta_records']),retained_beta_envelopes=retained,
        original_actual_error_union_and_shared_calibration_retained=True,
        old_verified_target_proofs_reused_by_hash=True,unexcluded_rectangles_claimed_feasible=False,
        new_global_upper_verified=True,sharp_global_optimum_claimed=False,
        absolute_cost_upper=str(upper),genuine_compatible_lower=str(lower),union_gap_ratio=str(upper/lower))
    (HERE/'artifact_check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['absolute_cost_upper','genuine_compatible_lower','union_gap_ratio']}))
if __name__=='__main__':main()
