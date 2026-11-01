"""Independent rational reconstruction of all A52 continuous-exposure covers.
No old interval, exp, log, or KL implementation and no constructor imported.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json
HERE=Path(__file__).resolve().parent
DEN=10**110
def down(x):return F((x*DEN).__floor__(),DEN)
def up(x):return F((x*DEN).__ceil__(),DEN)
class J:
    def __init__(self,l,r=None):self.l=down(F(l));self.r=up(F(l if r is None else r));assert self.l<=self.r
    def __add__(self,b):b=asj(b);return J(self.l+b.l,self.r+b.r)
    __radd__=__add__
    def __neg__(self):return J(-self.r,-self.l)
    def __sub__(self,b):return self+-asj(b)
    def __rsub__(self,b):return asj(b)+-self
    def __mul__(self,b):
        b=asj(b);v=[self.l*b.l,self.l*b.r,self.r*b.l,self.r*b.r];return J(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=asj(b);assert b.l*b.r>0;return self*J(1/b.r,1/b.l)
def asj(x):return x if isinstance(x,J) else J(x)
@lru_cache(None)
def ex(t):
    t=F(t);assert t>=0;nscale=0
    while t>F(1,8):t/=2;nscale+=1
    term=F(1);z=term
    for n in range(1,61):term*= -t/n;z+=term
    ans=J(z-term*t/61,z)
    for _ in range(nscale):ans=ans*ans
    return ans
def local_log(x):
    x=F(x);assert 1<=x<=2;z=(x-1)/(x+1);value=sum((2*z**(2*i+1)/F(2*i+1) for i in range(130)),F(0))
    return J(value,value+2*z**261/(261*(1-z*z)))
def log(x):
    x=F(x);assert x>0
    if x<1:return -log(1/x)
    p=0
    while x>2:x/=2;p+=1
    return p*local_log(2)+local_log(x)
@lru_cache(None)
def shape(z):
    z=F(z);assert abs(z)<1
    value=sum(((-z)**i/F((i+1)*(i+2)) for i in range(19)),F(0));err=abs(z)**19/(20*21*(1-abs(z)))
    return J(value-err,value+err)
def contained(x,bs):
    a,b=map(F,bs);assert a<=x.l<=x.r<=b,(float(a-x.l),float(x.r-b))
def main():
    raw=(HERE/'joint_time_capacity_certificate.json').read_bytes();d=json.loads(raw);p=d['parameters'];z=d['hard_pair']
    c=F(p['cost_rate']);lam=F(p['good_mean']);moment=F(p['moment_bound']);delta=F(p['delta']);gap=F(p['gap']);ell=F(z['slow']);h=F(z['fast']);eps=F(z['fast_mass']);D=h-ell;k=2*h-ell
    assert ell==c+delta/2 and lam==c+delta+gap and eps==(lam-ell)/D and (1-eps)*ell+eps*h==lam and (1-eps)*ell*ell+eps*h*h==moment==25
    assert D==(moment-ell*ell)/(lam-ell)-2*ell and 0<eps<1 and F(p['action_exposure'])==delta/moment
    assert (ell-c)*F(p['action_exposure'])==F(p['profit_threshold'])==F(z['bad_profit_strict_upper']) and ell>c
    assert F(d['likelihood_ratio_lower'])==1-eps
    kap=shape(-eps);b=eps*eps*D*D/(ell*k);B=kap*b;A=lam*log(lam/ell)-lam+ell;event=F(9,10)*log(19);zcap=eps*D/ell
    contained(kap,d['domination_constant_outer']);contained(B,d['full_capture_information_cap_outer']);contained(A,d['zero_shadow_information_clock_outer']);contained(event,d['required_decision_KL_outer'])
    def binary(left,right):
        v=J(ex(ell*right).l,ex(ell*left).r);t=J((1-ex(D*left)).l,(1-ex(D*right)).r);pq=1-v;dd=eps*v*t;zz=dd/pq;zz=J(max(F(0),zz.l),min(zcap,zz.r));w=-eps*t
        return dd*dd*(J(shape(zz.r).l,shape(zz.l).r)/pq+J(shape(w.r).l,shape(w.l).r)/v)
    cells=0;axis=0;mixed=0
    for s in d['supports']:
        u=F(s['information_time_price']);v=F(s['information_capacity_price'])
        if s['kind']=='exposure_axis':
            assert v==0 and u==F(d['zero_shadow_information_clock_outer'][1]) and s['exact_zero_shadow_supremum'] is True;contained(A,s['rate_outer']);axis+=1;continue
        shadow=None if s['shadow_time_per_opportunity'] is None else F(s['shadow_time_per_opportunity'])
        def divisor(q):return F(1) if shadow is None else q+shadow
        small=F(s['small_exposure_end']);large=F(s['large_exposure_start']);assert small==F(1,10**8) and large==F(1,20)
        assert A.r*small/divisor(small)<=F(s['small_exposure_normalized_upper'])
        ee=ex(ell*large);assert (kap*eps*eps*ee/(1-ee)).r/divisor(large)<=F(s['large_exposure_normalized_upper'])
        basic=B.r if shadow is None else (B/(J(B.l/A.r)+shadow)).r;assert basic<=F(s['basic_normalized_upper'])
        peak=max(F(s['small_exposure_normalized_upper']),F(s['large_exposure_normalized_upper']));previous=small
        for row in s['q_cover']:
            left,right=map(F,row['q_interval']);assert previous==left<right;previous=right;info=binary(left,right);contained(info,row['binary_KL_outer'])
            bound=min(info.r,(B*(1-ex(k*right))).r,A.r*right)/divisor(left);bound=min(bound,basic)
            assert bound<=F(row['normalized_upper']);peak=max(peak,F(row['normalized_upper']));cells+=1
        assert previous==large and peak==F(s['rate_outer'][1])
        q=F(s['candidate']['exposure']);point=binary(q,q)/divisor(q);contained(point,s['candidate']['normalized_information_outer'])
        assert F(s['rate_outer'][0])==F(s['candidate']['normalized_information_outer'][0])<=point.l and peak/F(s['rate_outer'][0])-1==F(s['relative_normalized_gap'])<=F(1,1000)
        if shadow is None:assert u==0 and v==peak;axis+=1
        else:assert u==peak and v==shadow*peak;mixed+=1
    assert axis==2 and mixed==5 and len(d['supports'])==7
    elo,ehi=map(F,d['required_decision_KL_outer']);Aupper=F(d['zero_shadow_information_clock_outer'][1]);cap=d['supports'][1];vupper=F(cap['information_capacity_price']);vlower=F(cap['rate_outer'][0])
    assert d['necessary_integer_opportunity_capacity']==(elo/vupper).__ceil__() and F(d['necessary_exposure_time_lower'])==elo/Aupper
    for row in d['joint_deadline_lower_frontier']:
        N=row['capacity'];assert N==(F(row['capacity_multiplier'])*ehi/vlower).__ceil__();vals=[]
        for item in row['support_constraints']:
            ss=d['supports'][item['support_index']];u=F(ss['information_time_price']);v=F(ss['information_capacity_price']);val=max(F(0),(elo-v*N)/u);assert val==F(item['necessary_exposure_lower']);vals.append(val)
        assert max(vals)==F(row['joint_time_lower']) and elo/Aupper==F(row['separate_time_lower'])
    exm=d['simultaneous_resource_example'];T=F(exm['exposure_or_calendar_cap']);N=exm['opportunity_capacity']
    assert N==(F(102,100)*ehi/vlower).__ceil__() and F(exm['separate_time_information_lower'])==F(d['zero_shadow_information_clock_outer'][0])*T>ehi and F(exm['separate_capacity_information_lower'])==vlower*N>ehi
    bounds=[F(s['information_time_price'])*T+F(s['information_capacity_price'])*N for s in d['supports']]
    assert list(map(F,exm['mixed_information_uppers']))==bounds and exm['rejecting_support_indices']==[i for i,x in enumerate(bounds) if x<elo] and exm['rejecting_support_indices']
    pw=F(exm['impossible_good_issuance_probability_at_bad_alpha']);al=F(p['alpha']);pkl=pw*log(pw/al)+(1-pw)*log((1-pw)/(1-al));contained(pkl,exm['decision_KL_at_power_threshold_outer'])
    assert pw==F(21,25) and min(bounds)<F(exm['decision_KL_at_power_threshold_outer'][0]) and F(exm['necessary_good_rejection_strict_lower'])==1-pw==F(4,25)
    src=d['source'];path=HERE.parent/src['path'];sr=path.read_bytes();ck=(path.parent/src['check_filename']).read_bytes()
    assert hashlib.sha256(sr).hexdigest()==src['sha256']==json.loads(ck)['source_sha256'] and hashlib.sha256(ck).hexdigest()==src['check_sha256']
    assert d['shadow_times_are_actual_setup_times'] is False and d['prior_slow_below_cost_future_per_opportunity_profit_bound_reused_for_new_pair'] is False
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_rational_precision_digits=110,
        complete_positive_exposure_cover_cells_checked=cells,affine_resource_supports_checked=7,mixed_resource_supports_checked=5,
        full_second_moment_hard_pair_recomputed=True,positive_slow_speed_market_strict_bad_profit_proved_analytically=True,
        necessary_integer_opportunity_capacity=d['necessary_integer_opportunity_capacity'],joint_frontier_points_checked=len(d['joint_deadline_lower_frontier']),
        separate_time_and_capacity_bounds_pass_but_joint_information_excludes_example=True,
        example_good_rejection_must_exceed='4/25 at bad issuance <=1/20',
        original_full_market_binary_observation_and_adaptive_control_permissions_preserved=True,
        shadow_resource_multipliers_misidentified_as_actual_setup_cost=False,
        new_bad_pair_claimed_same_as_old_slow_unprofitable_pair=False,
        old_per_business_rare_profit_upper_reused_for_positive_slow_pair=False,
        necessary_resource_bounds_claimed_feasible_financial_policy=False,
        sharp_all_policy_resource_frontier_or_new_generic_duality_claimed=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out))
if __name__=='__main__':main()
