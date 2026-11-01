"""Finite ordinary-qualification policy reaching the original 0.1 NPV task.
The independent verifier checks all enclosures and analytic inequality inputs.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json,hashlib
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(x):
    x=F(x);return Decimal(x.numerator)/Decimal(x.denominator)
def box(x):v=F(x);e=F(1,10**65);return [str(v-e),str(v+e)]
def kl(t,p):return t*(t/p).ln()+(1-t)*((1-t)/(1-p)).ln()
def main():
    parent=HERE.parent/'第六十二轮研究_20261029/uniform_value_ceiling.json';raw=json.loads(parent.read_text('utf8'));check=parent.parent/'artifact_check.json';assert sha(parent)==json.loads(check.read_text('utf8'))['source_sha256']
    business_parent=HERE.parent/'第六十一轮研究_20261029/positive_rate_point0875_certificate.json';business_check=business_parent.parent/'artifact_check.json';assert sha(business_parent)==json.loads(business_check.read_text('utf8'))['source_sha256']
    pa=dict(cost_rate='3/2',good_mean='301/200',moment_bound='25',delta='1/1000',gap='1/250',alpha='1/20',beta='1/20',valuation_rate='249/100000',
      learning_success_fee='1/100',action_exposure='1/25000',action_profit_threshold='1/50000000',project_NPV_target='1/10',
      business_fixed_access_fee='1/1000000000000',business_access_time='1/1000000000')
    control=dict(cash_rate_fraction='1999999/2000000',constant_business_reserve='603',business_exposure='7/10000000',business_windows=5_000_000_000,
      first_paid_learning_exposure='1/100000000',early_independent_certificate_probability='4999/100000',late_bank_seed='1/100000000',
      late_learning_absolute_start='15000',late_learning_exposure='1/6250',late_learning_windows=250_000_000_000,late_accept_success_count=60_080_000)
    with localcontext() as ctx:
        ctx.prec=160;c,lam,m2,rho=map(dec,[pa['cost_rate'],pa['good_mean'],pa['moment_bound'],pa['valuation_rate']]);alpha=dec(pa['alpha']);fraction=dec(control['cash_rate_fraction']);rate=rho*fraction
        a=dec(control['business_exposure']);delay=dec(pa['business_access_time']);d=a+delay;k=c*a+dec(pa['business_fixed_access_fee']);pg=lam*a-m2*a*a/2
        x=dec(control['constant_business_reserve']);N=control['business_windows'];v=(-rho*d).exp();yc=(x-k)*(rate*d).exp()-x
        value=(v*(yc+pg)/(1-v)-x)*(1-(-rho*d*N).exp())
        q=dec(control['first_paid_learning_exposure']);early=dec(control['early_independent_certificate_probability']);seed=dec(control['late_bank_seed']);W=dec(control['late_learning_absolute_start'])
        waitloss=x*(1-(-(rho-rate)*q).exp());lower=early*(-rho*q).exp()*value-c*q-seed-waitloss
        ql=dec(control['late_learning_exposure']);nl=control['late_learning_windows'];p0=(c+dec('1/1000'))*ql;p1=lam*ql-m2*ql*ql/2;t=dec(F(control['late_accept_success_count'],nl));tail_bad=(-nl*kl(t,p0)).exp();tail_good=(-nl*kl(t,p1)).exp()
        late_budget=nl*c*ql+c*dec(pa['action_exposure']);seed_need=late_budget*(-rate*W).exp()
        out=dict(sources={'A62':dict(path='第六十二轮研究_20261029/uniform_value_ceiling.json',sha256=sha(parent),check_sha256=sha(check)),
                         'A61':dict(path='第六十一轮研究_20261029/positive_rate_point0875_certificate.json',sha256=sha(business_parent),check_sha256=sha(business_check))},constructor_sha256=sha(Path(__file__)),parameters=pa,control=control,
          business_window_duration=str(F(d)),business_prepayment=str(F(k)),business_good_success_probability_lower=str(F(pg)),business_failure_state_refund_outer=box(yc),
          business_finite_NPV_lower_outer=box(value),first_window_reserve_discount_loss_outer=box(waitloss),whole_good_NPV_lower_outer=box(lower),
          late_bad_success_probability_upper=str(F(p0)),late_good_success_probability_lower=str(F(p1)),late_success_frequency_threshold=str(F(t)),
          late_bad_error_chernoff_outer=box(tail_bad),late_good_error_chernoff_outer=box(tail_good),whole_bad_issuance_upper_outer=box(early+(1-early)*tail_bad),whole_good_rejection_upper_outer=box((1-early)*tail_good),
          late_worst_case_prepayment_budget=str(F(late_budget)),late_discounted_seed_needed_outer=box(seed_need),
          whole_initial_equity=str(F(x+c*q+seed)),business_calendar=str(F(N*d)),
          whole_opportunity_capacity=1+max(N,nl),whole_final_calendar=str(F(W+nl*ql)),
          ordinary_paid_positive_window_requirement_preserved=True,actual_capture_before_every_certificate=False,
          original_alpha_beta_and_profit_target_changed=False,learning_success_fee_removed=False,business_positive_fee_or_delay_removed=False,
          finite_policy_claim=True,practical_opportunity_supply_claim=False,exact_critical_cash_rate_claim=False,
          ordinary_uniform_upper_inherited=raw['ordinary_uniform_value_upper_outer'],published_stage_snapshot_overwritten=False)
        out['ordinary_value_upper_to_new_lower_ratio_outer']=box(dec(raw['ordinary_uniform_value_upper_outer'][1])/lower)
        assert yc>0 and p0<t<p1 and early+(1-early)*tail_bad<alpha and (1-early)*tail_good<dec(pa['beta']) and seed_need<seed
        print(json.dumps(dict(cash_fraction=float(fraction),business_value=float(value),whole_NPV_lower=float(lower),bad_error_upper=float(early+(1-early)*tail_bad),late_seed_needed=float(seed_need),target_reached=lower>dec('1/10')),ensure_ascii=False))
        assert lower>dec('1/10')
    (HERE/'ordinary_target_reserve.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')

if __name__=='__main__':main()
