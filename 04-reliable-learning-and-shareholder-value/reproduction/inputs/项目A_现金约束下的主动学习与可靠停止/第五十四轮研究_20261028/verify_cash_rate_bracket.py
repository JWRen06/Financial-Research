"""Independent A54 negative potentials, exact cash grid and full NPVs."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_CEILING
import numpy as np,json,hashlib,math,time
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dd(x):x=F(x);return Decimal(x.numerator)/Decimal(x.denominator)
def eb(t):
    s=F(1);v=s
    for n in range(1,11):s*= -t/n;v+=s
    return v-t**11/math.factorial(11),v
def main():
    starttime=time.monotonic();raw=(HERE/'cash_rate_bracket_certificate.json').read_bytes();d=json.loads(raw);npath=HERE/d['negative']['path'];ppath=HERE/d['positive']['path'];assert sha(npath)==d['negative']['sha256'] and sha(ppath)==d['positive']['sha256'];neg=json.loads(npath.read_text(encoding='utf-8'));pos=json.loads(ppath.read_text(encoding='utf-8'));sources={}
    for k,s in neg['sources'].items():
        p=HERE.parent/s['path'];cc=p.parent/'artifact_check.json';assert sha(p)==s['sha256']==json.loads(cc.read_text(encoding='utf-8'))['source_sha256'] and sha(cc)==s['check_sha256'];sources[k]=json.loads(p.read_text(encoding='utf-8'))
    assert pos['source_A51']==neg['sources']['A51'] and pos['parameters']==sources['A51']['parameters'] and pos['late_full_market_capture_learner']==sources['A51']['late_full_market_capture_learner']
    assert sha(HERE/'cash_rate_refinement.json')==neg['proposal_source_sha256'];C=F('1.5');lam=F('1.505');rho=F('.00249');m2=F(25);ds=F('.01')-1+C/lam;checks=[]
    for cert in neg['certificates']:
        rate=F(cert['cash_rate_fraction']);rc=rho*rate;spread=rho-rc;h=cert['increment_numerators'];S=cert['grid_step_denominator'];D=cert['increment_denominator'];U=len(h)-1;assert (S,D,U)==(1000,10**14,3000) and h[0]==0 and all(0<=x<=y for x,y in zip(h,h[1:]));ma=mi=mf=None
        for i in range(U):
            slope=F((h[i+1]-h[i])*S,D)
            for j in [i,i+1]:
                x=F(j,S);hj=F(h[j],D);hit=F(h[min(j+S,U)]-h[j],D);av=-(lam-C-spread*x+(rc*x-C)*slope+lam*hit-rho*hj);iv=-(-spread*x+rc*x*slope-rho*hj);ma=av if ma is None else min(ma,av);mi=iv if mi is None else min(mi,iv)
        for j in range(990,U+1):
            val=-(lam*(F(h[j],D)-ds)-spread*(F(j,S)-F('.99')));mf=val if mf is None else min(mf,val)
        at=-(lam-C-spread*3-rho*F(h[-1],D));ft=-(lam*(F(h[-1],D)-ds)-spread*(3-F('.99')));eps=F(cert['strict_first_capture_loss']);assert min(ma,mi,at)>=0 and mf>=eps>0 and ft>=eps and eps==mf/2
        z=cert['exact_check'];assert [str(ma),str(mi),str(mf),str(at),str(ft)]==[z[k] for k in ['minimum_active_slack','minimum_idle_slack','minimum_first_capture_slack','active_tail_slack','first_tail_slack']];checks.append(rate)
    assert len(checks)==4 and max(checks)==F(329,4000)==F(d['negative_cash_rate_fraction_upper'])
    rate=F(pos['cash_rate_fraction']);rc=F(pos['cash_rate']);a=F(pos['optional_business_exposure']);fee=F(pos['business_fixed_access_fee']);sigma=F(pos['business_access_time']);duration=a+sigma;k=C*a+fee;p=lam*a-m2*a*a/2
    assert rate==F(9,100) and rc==rho*rate and a==F(1,10**6) and fee==F(1,10**12) and sigma==F(1,10**9) and F(pos['certification_action_exposure'])==F(1,25000) and duration==F(pos['business_total_window']) and p==F(pos['business_capture_probability_lower'])
    D=pos['cash_grid_denominator'];E=pos['exponential_denominator'];assert D==pos['value_denominator']==10**16 and E==10**18;vn=pos['upward_bank_discount_numerator'];dn=pos['downward_shareholder_discount_numerator'];assert eb(rc*duration)[1]<=F(vn,E) and F(dn,E)<=eb(rho*duration)[0]
    grid=HERE/d['cash_grid']['path'];assert sha(grid)==d['cash_grid']['sha256']==pos['cash_grid_sha256'] and grid.name==pos['cash_grid_file']
    with np.load(grid,allow_pickle=False) as f:X=f['cash'].copy();W=f['value'].copy()
    U=pos['nonzero_cash_states'];assert U==1466908 and len(X)==len(W)==U+1 and X[0]==W[0]==0 and int(W.min())>=0;dest=np.minimum(np.searchsorted(X,X[:-1]+D,side='right')-1,U);resmin=jumpmin=None
    for i in range(1,U+1):
        assert int(X[i])==int(k*D)+(vn*int(X[i-1])+E-1)//E and X[i]>X[i-1];j=int(dest[i-1]);y=int(X[i-1])+D;pay=y-int(X[j]);assert i<=j<=U and pay>=0 and (j==U or y<int(X[j+1]));before=int(W[i-1]);hit=int(W[j])+pay;assert hit>=before;res=dn*((p.denominator-p.numerator)*before+p.numerator*hit)-int(W[i])*E*p.denominator;assert res>=0;resmin=res if resmin is None else min(resmin,res);jumpmin=hit-before if jumpmin is None else min(jumpmin,hit-before)
    assert str(resmin)==pos['minimum_exact_Bellman_residual_numerator'] and F(jumpmin,D)==F(pos['minimum_exact_capture_jump_gain']);ix=pos['initial_state'];M=F(int(X[ix]),D)-F('.99');tail=F(pos['discounted_tail_allowance']);V=F(int(W[ix]),D)-tail;H=pos['finite_business_steps'];assert M==F(pos['initial_supplement'])>=0 and V==F(pos['finite_business_payout_lower']) and F(int(X[ix]),D)==F(pos['initial_business_cash']) and F(int(W[ix]),D)==F(pos['initial_potential']) and F(int(W.max()),D)==F(pos['maximum_potential']) and F(int(X[-1]),D)==F(pos['barrier_cash'])
    assert H==10**10 and duration*H==F(pos['business_calendar_upper']);term=F(1);el=term
    for n in range(1,121):term*=rho*duration*H/n;el+=term
    assert F(int(W.max()),D)/el<=tail and lam*(V-M)-C-(rho-rc)*M==F(pos['full_policy_first_order_coefficient'])>0
    with localcontext() as ctx:
        ctx.prec=180;late=pos['late_full_market_capture_learner'];qt=Decimal('.004')/(10*dd(m2));p0=Decimal('1.501')*qt;pg=dd(lam)*qt-dd(m2)*qt*qt/2;theta=(p0+dd(lam)*qt)/2
        def kl(x,y):return x*(x/y).ln()+(1-x)*((1-x)/(1-y)).ln()
        d0=kl(theta,p0);d1=kl(theta,pg);N=late['trial_count'];cut=late['success_count_cutoff'];assert N==int(max(Decimal(40).ln()/d0,Decimal(20).ln()/d1).to_integral_value(rounding=ROUND_CEILING)) and cut==int((N*theta).to_integral_value(rounding=ROUND_CEILING));K=N*dd(C)*qt+dd(C*F(1,25000));assert K==dd(late['all_path_prepaid_budget']) and (-N*d0).exp()<=Decimal('.025') and (-N*d1).exp()<=Decimal('.05')
        for z in pos['positive_full_market_policies']:
            q=dd(z['paid_probe_exposure']);good=dd(lam)*q-dd(m2)*q*q/2;bad=Decimal('1.501')*q;deposit=q*q/100;wait=z['late_waiting_time'];B=dd(M)+deposit+dd(C)*q;assert B==dd(z['initial_equity']) and deposit==dd(z['late_initial_deposit']) and good==dd(z['whole_good_probe_success_lower']) and bad==dd(z['whole_bad_probe_success_upper']) and wait==int(((K/deposit).ln()/dd(rc)).to_integral_value(rounding=ROUND_CEILING)) and deposit*(dd(rc)*wait).exp()>=K
            early=(-dd(rho)*q).exp()*((dd(M)+deposit)*(dd(rc)*q).exp()-dd(M)+dd(V));failed=dd(M)*(-dd(rho-rc)*q).exp();value=good*early+(1-good)*failed-B;assert 0<dd(z['whole_good_equity_lower_outer'][0])<=value<=dd(z['whole_good_equity_lower_outer'][1]) and early>failed
            assert bad+(1-bad)*Decimal('.025')==dd(z['full_bad_issuance_upper'])<=Decimal('.05') and (1-good)*Decimal('.05')==dd(z['full_good_rejection_upper'])<=Decimal('.05');bf=F(z['early_failure_evidence_factor']);bh=F(z['early_capture_evidence_factor']);assert bf==F(1,2) and bh==(1+1/F(bad))/2 and (1-F(bad))*bf+F(bad)*bh==1 and bh>=20
            assert dd(z['late_certificate_deadline'])==Decimal(wait)+N*qt and dd(z['all_cash_terminated_by'])==max(q+dd(duration)*H,Decimal(wait)+N*qt+dd(F(1,25000)))
    assert len(pos['positive_full_market_policies'])==3 and d['critical_fraction_bracket']==['329/4000','9/100'] and d['positive_interval_by_virtual_account']==['9/100','1']
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_decimal_digits=180,negative_piecewise_cash_certificates_checked=4,negative_cash_endpoint_inequalities_checked=4*(6000+6000+2011),negative_rate_fraction_upper='329/4000',positive_rate_fraction_lower='9/100',critical_fraction_bracket=['329/4000','9/100'],exact_cash_grid_transitions_checked=U,exact_Bellman_inequalities_checked=U,exact_capture_direction_inequalities_checked=U,full_market_positive_capture_policies_checked=3,late_learner_and_finite_business_tail_checked=True,
        cash_rate_monotonic_negative_extension_and_positive_virtual_account_extension_proved=True,once_equity_learning_fee_actual_capture_alpha_beta_kept=True,positive_access_fee_time_same_as_A51=True,
        failed_finite_LP_family_claimed_global_result=False,exact_frontier_point_one_profit_or_realistic_resources_claimed=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
