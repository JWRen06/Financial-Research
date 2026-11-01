"""Independent exact cash transition audit and rational-interval full NPV."""
from pathlib import Path
from fractions import Fraction as F
import numpy as np,json,hashlib,math,sys,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第五十二轮研究_20261028'))
import verify_joint_time_capacity as interval
interval.DEN=10**180;interval.ex.cache_clear()
J=interval.J;ex=interval.ex
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def eb(t):
    term=z=F(1)
    for n in range(1,13):term*= -t/n;z+=term
    return z-t**13/math.factorial(13),z
def main():
    start=time.monotonic();raw=(HERE/'positive_rate_point0875_certificate.json').read_bytes();d=json.loads(raw);s=d['source_A51'];source=HERE.parent/s['path'];cc=source.parent/'artifact_check.json';assert sha(source)==s['sha256']==json.loads(cc.read_text(encoding='utf-8'))['source_sha256'] and sha(cc)==s['check_sha256'];src=json.loads(source.read_text(encoding='utf-8'));assert d['parameters']==src['parameters'] and d['late_full_market_capture_learner']==src['late_full_market_capture_learner']
    rho=F(249,100000);rate=F(7,80);rc=rho*rate;c=F(3,2);lam=F(301,200);m2=F(25);a=F(1,2000000);fee=F(1,10**12);sigma=F(1,10**9);duration=a+sigma;k=c*a+fee;p=lam*a-m2*a*a/2
    assert F(d['valuation_rate'])==rho and F(d['cash_rate_fraction'])==rate and F(d['cash_rate'])==rc and F(d['optional_business_exposure'])==a<F(d['certification_action_exposure'])==F(1,25000)
    assert F(d['business_fixed_access_fee'])==fee and F(d['business_access_time'])==sigma and F(d['business_total_window'])==duration and F(d['business_capture_probability_lower'])==p
    D=d['cash_grid_denominator'];E=d['exponential_denominator'];assert D==d['value_denominator']==10**18 and E==10**24;vn=d['upward_bank_discount_numerator'];dn=d['downward_shareholder_discount_numerator'];assert eb(rc*duration)[1]<=F(vn,E) and F(dn,E)<=eb(rho*duration)[0]
    grid=HERE/d['cash_grid_file'];assert sha(grid)==d['cash_grid_sha256']
    with np.load(grid,allow_pickle=False) as npz:X=npz['cash'].copy();W=npz['value'].copy()
    U=d['nonzero_cash_states'];assert U==2933800 and len(X)==len(W)==U+1 and X[0]==W[0]==0 and int(W.min())>=0
    destinations=np.minimum(np.searchsorted(X,X[:-1]+D,side='right')-1,U);resmin=jumpmin=None
    for i in range(1,U+1):
        oldcash=int(X[i-1]);assert int(X[i])==int(k*D)+(vn*oldcash+E-1)//E and X[i]>X[i-1]
        j=int(destinations[i-1]);hitcash=oldcash+D;pay=hitcash-int(X[j]);assert i<=j<=U and pay>=0 and (j==U or hitcash<int(X[j+1]))
        before=int(W[i-1]);after=int(W[j])+pay;jump=after-before;assert jump>=0
        residual=dn*((p.denominator-p.numerator)*before+p.numerator*after)-int(W[i])*E*p.denominator;assert residual>=0
        resmin=residual if resmin is None else min(resmin,residual);jumpmin=jump if jumpmin is None else min(jumpmin,jump)
    assert str(resmin)==d['minimum_exact_Bellman_residual_numerator'] and F(jumpmin,D)==F(d['minimum_exact_capture_jump_gain'])
    ix=d['initial_state'];M=F(int(X[ix]),D)-F(99,100);tail=F(d['discounted_tail_allowance']);V=F(int(W[ix]),D)-tail;H=d['finite_business_steps'];assert M>=0 and M==F(d['initial_supplement']) and V==F(d['finite_business_payout_lower']) and F(int(X[ix]),D)==F(d['initial_business_cash']) and F(int(W[ix]),D)==F(d['initial_potential']) and F(int(W.max()),D)==F(d['maximum_potential']) and F(int(X[-1]),D)==F(d['barrier_cash'])
    assert H==20_000_000_000 and duration*H==F(d['business_calendar_upper']) and tail==F(1,10**8)
    term=el=F(1)
    for n in range(1,121):term*=rho*duration*H/n;el+=term
    assert F(int(W.max()),D)/el<=tail and lam*(V-M)-c-(rho-rc)*M==F(d['full_policy_first_order_coefficient'])>0
    late=d['late_full_market_capture_learner'];budget=F(late['all_path_prepaid_budget']);calendar=F(late['certificate_deadline_without_wait']);count=0;npvl=[]
    for row in d['positive_full_market_policies']:
        q=F(row['paid_probe_exposure']);good=lam*q-m2*q*q/2;bad=F(1501,1000)*q;deposit=q*q/100;B=M+deposit+c*q;wait=row['late_waiting_time']
        assert q in [F(1,10**8),F(1,10**9),F(1,10**10)] and B==F(row['initial_equity']) and deposit==F(row['late_initial_deposit']) and good==F(row['whole_good_probe_success_lower']) and bad==F(row['whole_bad_probe_success_upper']) and wait>q
        assert budget*ex(rc*wait).r<=deposit<budget*ex(rc*(wait-1)).l
        early=ex(rho*q)*(J(M+deposit)/ex(rc*q)-M+V);failed=M*ex((rho-rc)*q);value=good*early+(1-good)*failed-B
        vl,vh=map(F,row['whole_good_equity_lower_outer']);assert 0<vl<=value.l<=value.r<=vh and early.l>failed.r
        assert bad+(1-bad)*F(1,40)==F(row['full_bad_issuance_upper'])<=F(1,20) and (1-good)*F(1,20)==F(row['full_good_rejection_upper'])<=F(1,20)
        bf=F(row['early_failure_evidence_factor']);bh=F(row['early_capture_evidence_factor']);assert bf==F(1,2) and bh==(1+1/bad)/2 and (1-bad)*bf+bad*bh==1 and bh>=20
        assert F(row['late_certificate_deadline'])==wait+calendar and F(row['all_cash_terminated_by'])==max(q+duration*H,F(wait)+calendar+F(1,25000));npvl.append(float(vl));count+=1
    a54=HERE.parent/'第五十四轮研究_20261028/cash_rate_bracket_certificate.json';old=json.loads(a54.read_text(encoding='utf-8'));assert sha(a54)==json.loads((a54.parent/'artifact_check.json').read_text(encoding='utf-8'))['source_sha256'] and old['negative_cash_rate_fraction_upper']==d['negative_rate_bound_unchanged']=='329/4000'
    assert count==3 and d['candidate_positive_interval_by_virtual_account']==['7/80','1'] and not d['point_one_target_or_realistic_resource_claimed'] and not d['old_A54_A55_sources_overwritten']
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_rational_interval_digits=180,exact_cash_grid_transitions_checked=U,exact_Bellman_inequalities_checked=U,exact_capture_direction_inequalities_checked=U,finite_business_tail_and_positive_initial_capture_coefficient_checked=True,original_positive_access_fee_and_delay_preserved=True,original_A51_full_market_learner_and_cash_proof_reused_explicitly=True,full_market_positive_capture_policies_checked=count,rational_full_NPV_intervals_positive=True,display_NPV_lowers=npvl,positive_rate_fraction_lower='7/80',negative_rate_fraction_upper_unchanged='329/4000',critical_fraction_bracket=['329/4000','7/80'],higher_rate_virtual_account_extension_proved_analytically=True,A54_negative_source_sha256=sha(a54),original_point_one_target_result_or_realistic_resources_claimed=False,old_certificates_overwritten=False,seconds=time.monotonic()-start)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
