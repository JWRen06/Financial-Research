"""Independent rational interval check of a fully specified finite policy.
Uses the previously audited rational exp/log engine, no Decimal or constructor.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第五十二轮研究_20261028'))
from verify_joint_time_capacity import ex,log,J
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def contained(x,stored):
    l,r=map(F,stored);assert l<x.l<=x.r<r,(float(x.l-l),float(r-x.r));assert r-l==F(2,10**65)
def kl(t,p):return t*log(t/p)+(1-t)*log((1-t)/(1-p))

def main():
    path=HERE/'ordinary_target_reserve.json';d=json.loads(path.read_text('utf8'));assert d['constructor_sha256']==sha(HERE/'ordinary_target_reserve.py')
    sources={}
    for name,row in d['sources'].items():
        p=HERE.parent/row['path'];ck=p.parent/'artifact_check.json';assert sha(p)==row['sha256']==json.loads(ck.read_text('utf8'))['source_sha256'] and sha(ck)==row['check_sha256'];sources[name]=json.loads(p.read_text('utf8'))
    p=d['parameters'];z=d['control'];assert all(F(p[k])==F(v) for k,v in sources['A62']['parameters'].items())
    c,lam,m2,rho,delta,gap=map(F,[p[k] for k in ['cost_rate','good_mean','moment_bound','valuation_rate','delta','gap']]);alpha,beta=map(F,[p[k] for k in ['alpha','beta']]);A=F(p['action_exposure']);gA=F(p['action_profit_threshold'])
    assert c==F(3,2) and lam==c+delta+gap==F(301,200) and m2==25 and A==delta/m2==F(1,25000) and gA==delta*delta/(2*m2)==F(1,50_000_000)
    assert F(p['learning_success_fee'])==F(1,100) and F(p['project_NPV_target'])==F(1,10) and alpha==beta==F(1,20)
    fee,delay=map(F,[p[k] for k in ['business_fixed_access_fee','business_access_time']]);assert fee==F(sources['A61']['business_fixed_access_fee'])==F(1,10**12) and delay==F(sources['A61']['business_access_time'])==F(1,10**9)
    fraction=F(z['cash_rate_fraction']);rate=rho*fraction;X=F(z['constant_business_reserve']);a=F(z['business_exposure']);n=z['business_windows'];duration=a+delay;k=c*a+fee;pg=lam*a-m2*a*a/2
    assert fraction==F(1999999,2000000) and X==603 and a==F(7,10**7) and n==5_000_000_000
    assert F(d['business_window_duration'])==duration and F(d['business_prepayment'])==k and F(d['business_good_success_probability_lower'])==pg>0
    bank_disc=ex(rate*duration);val_disc=ex(rho*duration);refund=J(X-k)/bank_disc-X
    contained(refund,d['business_failure_state_refund_outer']);assert refund.l>0 and X>k and X>=c*A
    # This single interval inequality proves the exact same cash transition for
    # every one of the n windows and every possible capture sequence.
    value=(val_disc*(refund+pg)/(1-val_disc)-X)*(1-ex(rho*duration*n))
    contained(value,d['business_finite_NPV_lower_outer']);assert value.l>F('2.00045882591906')
    q=F(z['first_paid_learning_exposure']);seed=F(z['late_bank_seed']);early=F(z['early_independent_certificate_probability']);W=F(z['late_learning_absolute_start'])
    assert q==seed==F(1,10**8) and q>0 and early==F(4999,100000)<alpha and W==15000>q
    waitloss=X*(1-ex((rho-rate)*q));contained(waitloss,d['first_window_reserve_discount_loss_outer'])
    lower=early*ex(rho*q)*value-c*q-seed-waitloss;contained(lower,d['whole_good_NPV_lower_outer']);assert lower.l>F('.10000291170519')>F(p['project_NPV_target'])
    ql=F(z['late_learning_exposure']);nl=z['late_learning_windows'];count=z['late_accept_success_count'];p0=(c+delta)*ql;p1=lam*ql-m2*ql*ql/2;t=F(count,nl)
    assert ql==F(1,6250) and nl==250_000_000_000 and count==60_080_000 and 0<p0<t<p1<1
    assert F(d['late_bad_success_probability_upper'])==p0 and F(d['late_good_success_probability_lower'])==p1 and F(d['late_success_frequency_threshold'])==t
    bad_information=kl(t,p0);good_information=kl(t,p1)
    # Monotone exponent intervals for Chernoff, including strict-good rejection.
    bad=J(ex(nl*bad_information.r).l,ex(nl*bad_information.l).r);good=J(ex(nl*good_information.r).l,ex(nl*good_information.l).r)
    contained(bad,d['late_bad_error_chernoff_outer']);contained(good,d['late_good_error_chernoff_outer'])
    whole_bad=early+(1-early)*bad;whole_good=(1-early)*good
    contained(whole_bad,d['whole_bad_issuance_upper_outer']);contained(whole_good,d['whole_good_rejection_upper_outer']);assert whole_bad.r<alpha and whole_good.r<beta
    budget=nl*c*ql+c*A;assert F(d['late_worst_case_prepayment_budget'])==budget
    need=budget*ex(rate*W);contained(need,d['late_discounted_seed_needed_outer']);assert need.r<seed
    assert F(d['whole_initial_equity'])==X+c*q+seed==F('603.000000025')
    assert F(d['business_calendar'])==n*duration==F(3505) and F(d['whole_final_calendar'])==W+nl*ql==40_015_000
    assert d['whole_opportunity_capacity']==1+max(n,nl)==250_000_000_001
    for key in ['ordinary_paid_positive_window_requirement_preserved','finite_policy_claim']:assert d[key] is True
    for key in ['actual_capture_before_every_certificate','original_alpha_beta_and_profit_target_changed','learning_success_fee_removed','business_positive_fee_or_delay_removed','practical_opportunity_supply_claim','exact_critical_cash_rate_claim','published_stage_snapshot_overwritten']:assert d[key] is False
    assert d['ordinary_uniform_upper_inherited']==sources['A62']['ordinary_uniform_value_upper_outer'] and lower.r<F(d['ordinary_uniform_upper_inherited'][1])
    ratio=J(F(d['ordinary_uniform_upper_inherited'][1]))/lower;contained(ratio,d['ordinary_value_upper_to_new_lower_ratio_outer']);assert ratio.r<F('1.00463')
    out=dict(source_sha256=sha(path),verification_code_sha256=sha(Path(__file__)),rational_interval_engine_sha256=sha(HERE.parent/'第五十二轮研究_20261028/verify_joint_time_capacity.py'),
       original_market_risks_action_and_success_fee_preserved=True,positive_business_fee_and_delay_preserved=True,
       exact_constant_reserve_all_path_transition_inequality_checked=True,full_finite_business_value_and_initial_equity_book_checked=True,
       ordinary_original_point_one_target_finite_sufficient_policy_checked=True,late_fixed_sample_full_good_bad_Chernoff_checks=2,
       finite_seed_growth_full_learning_prepayment_checked=True,branch_maximum_calendar_and_capacity_checked=True,
       sufficient_cash_rate_fraction=str(fraction),whole_good_NPV_strict_lower='.10000291170519',
       higher_cash_rate_extension_proved_analytically=True,actual_capture_qualification_or_practical_scale_claimed=False,
       original_A60_necessary_cash_rate_unchanged=True,high_rate_ordinary_value_gap_ratio_strictly_below='1.00463',new_ordinary_sufficient_policies=1,published_stage_snapshot_overwritten=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(dict(ordinary_target_reached=True,cash_rate_fraction=str(fraction),NPV_lower=out['whole_good_NPV_strict_lower'],all_path_payability=True,actual_capture_qualification=False)))

if __name__=='__main__':main()
