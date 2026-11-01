"""Independent, 110-digit rational interval reconstruction of A64.
No Decimal and no constructor import. Finite policy formulas and cash invariants
are checked; general Taylor/coupling/self-financing arguments are analytic.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

engine = HERE.parent/'第五十二轮研究_20261028/verify_joint_time_capacity.py'
old_check = json.loads((HERE.parent/'第六十三轮研究_20261029/artifact_check.json').read_text('utf8'))
assert sha(engine) == old_check['rational_interval_engine_sha256']
sys.path.insert(0, str(engine.parent))
from verify_joint_time_capacity import J, ex, log

def contained(x, stored):
    left, right = map(F, stored)
    assert left < x.l <= x.r < right, (float(x.l-left), float(right-x.r))
    assert right-left == F(2, 10**65)

def neg_exp_interval(x):
    assert x.l >= 0
    return J(ex(x.r).l, ex(x.l).r)

def power(p, n):
    assert 0 < p <= 1 and n >= 0
    return neg_exp_interval(-n*log(p))

def kl(t, p):
    return t*log(t/p)+(1-t)*log((1-t)/(1-p))

def main():
    path = HERE/'capture_target_approach.json'
    d = json.loads(path.read_text('utf8'))
    assert d['constructor_sha256'] == sha(HERE/'capture_target_approach.py')
    loaded = {}
    for key, row in d['sources'].items():
        source = HERE.parent/row['path']
        check = source.parent/'artifact_check.json'
        assert sha(source) == row['sha256'] == json.loads(check.read_text('utf8'))['source_sha256']
        assert sha(check) == row['check_sha256']
        loaded[key] = json.loads(source.read_text('utf8'))
    p, z = d['parameters'], d['control']
    assert p == loaded['A63']['parameters']
    c, lam, m2, rho, A, gA, fee, delay, s, alpha, beta = [F(p[k]) for k in
        ['cost_rate','good_mean','moment_bound','valuation_rate','action_exposure',
         'action_profit_threshold','business_fixed_access_fee','business_access_time',
         'learning_success_fee','alpha','beta']]
    assert c == F(3,2) and lam == F(301,200) and m2 == 25 and rho == F(249,100000)
    assert A == F(1,25000) and gA == F(1,50_000_000) and s == F(1,100)
    assert alpha == beta == F(1,20) and fee == F(1,10**12) and delay == F(1,10**9)
    mbad = (c+gA/A+m2*A/2+c*m2*A*A/6)/(1+c*A/2)
    assert F(d['strict_bad_mean_upper']) == mbad < c+F(p['delta'])
    # Rearrange the analytic integrated second-order Taylor bound exactly.
    assert (mbad-c)*A-m2*A*A/2+c*mbad*A*A/2-c*m2*A**3/6 == gA
    q = F(z['early_learning_exposure']); ne = z['early_learning_windows']; risk = F(z['early_bad_risk_allocation'])
    assert q == F(1,10**9) and ne == 34_173_771 and risk == F(49_999_999,10**9) < alpha
    p0 = mbad*q; p1 = lam*q-m2*q*q/2; T = ne*q; L = c*T
    assert F(d['early_bad_success_probability_upper']) == p0
    assert F(d['early_good_success_probability_lower']) == p1 > p0 > 0
    assert F(d['early_maximum_duration']) == T == F('0.034173771')
    assert F(d['early_worst_case_invoice_budget']) == L
    early_bad = 1-power(1-p0,ne); early_next = 1-power(1-p0,ne+1)
    early_good = 1-power(1-p1,ne)
    for value, field in [(early_bad,'early_bad_certificate_probability_outer'),
                         (early_next,'early_bad_probability_if_one_more_window_outer'),
                         (early_good,'early_good_certificate_probability_lower_outer')]:
        contained(value,d[field])
    assert early_bad.r < risk < early_next.l
    vq = ex(rho*q)
    surviving_power = neg_exp_interval(ne*(rho*q-log(1-p1)))
    S = (1-surviving_power)/(1-vq*(1-p1))
    hit = vq*p1*S; invoices = c*q*S
    contained(hit,d['early_discounted_hit_probability_lower_outer'])
    contained(invoices,d['early_discounted_invoice_cost_upper_outer'])
    a = F(z['business_exposure']); nb = z['business_windows']; X = F(z['constant_business_reserve'])
    assert a == F(73,10**8) <= A and nb == 10_000_000_000 and X == 603
    duration = a+delay; k = c*a+fee; pg = lam*a-m2*a*a/2
    unused = c*(lam*a*a/2-m2*a**3/6)
    for field, value in [('business_window_duration',duration),('business_prepayment',k),
            ('business_good_success_probability_lower',pg),
            ('business_expected_unused_exposure_refund_lower',unused)]:
        assert F(d[field]) == value > 0
    assert X >= c*A and X > k and 0 < pg < 1
    ql = F(z['late_learning_exposure']); nl = z['late_learning_windows']; count = z['late_accept_success_count']
    assert ql == F(1,6250) and nl == 500_000_000_000 and count == 120_160_000 > 0
    t = F(count,nl); lp0 = (c+F(p['delta']))*ql; lp1 = lam*ql-m2*ql*ql/2
    for field,value in [('late_success_frequency_threshold',t),('late_bad_success_probability_upper',lp0),('late_good_success_probability_lower',lp1)]:
        assert F(d[field]) == value
    assert 0 < lp0 < t < lp1 < 1
    bad = neg_exp_interval(nl*kl(t,lp0)); good = neg_exp_interval(nl*kl(t,lp1))
    contained(bad,d['late_bad_error_chernoff_outer']); contained(good,d['late_good_error_chernoff_outer'])
    whole_bad = early_bad+(1-early_bad)*bad; whole_good = (1-early_good)*good
    contained(whole_bad,d['whole_bad_certificate_probability_upper_outer'])
    contained(whole_good,d['whole_good_rejection_probability_upper_outer'])
    assert whole_bad.r < alpha and whole_good.r < beta
    seed = F(z['late_bank_seed']); W = F(z['late_learning_absolute_start'])
    budget = nl*c*ql+c*A
    assert seed == F(1,10**8) and W == 15000 > T
    assert F(d['late_worst_case_prepayment_budget']) == budget == F('120000000.00006')
    rates = list(map(F,z['cash_rate_fractions']))
    assert rates == [F(9999999,10**7),F(99999999,10**8),F(1)]
    bounds = [F('.09999105'),F('.09999378'),F('.09999408')]
    v = ex(rho*duration)
    upper = F(loaded['A62']['actual_capture_uniform_value_upper'])
    assert F(d['inherited_actual_capture_uniform_upper']) == upper
    for row,fraction,strict_lower in zip(d['rate_menu'],rates,bounds,strict=True):
        assert F(row['cash_rate_fraction']) == fraction
        rate = rho*fraction
        y = J(X-k)/ex(rate*duration)-X
        contained(y,row['failure_window_surplus_outer']); assert y.l > 0
        value = (v*(y+pg+unused)/(1-v)-X)*(1-ex(rho*duration*nb))
        contained(value,row['business_finite_NPV_outer'])
        spread = J((rho-rate)*(X+L)*T)/ex(rate*T)
        contained(spread,row['early_cash_spread_loss_upper_outer'])
        lower = hit*(value+1-s)-invoices-spread-seed
        contained(lower,row['whole_good_NPV_lower_outer'])
        assert lower.l > strict_lower > F('.09999') and lower.r < F(1,10)
        needed = budget*ex(rate*W)
        contained(needed,row['late_discounted_seed_needed_outer']); assert needed.r < seed
        ratio = J(upper)/lower
        contained(ratio,row['inherited_actual_capture_upper_to_policy_lower_ratio_outer'])
        assert ratio.r < F('1.00370')
    assert F(d['whole_initial_equity']) == X+L+seed == F('603.0512606665')
    assert F(d['business_calendar']) == nb*duration == 7310
    assert d['whole_opportunity_capacity'] == ne+max(nb,nl) == 500_034_173_771
    assert F(d['whole_final_calendar']) == max(T+nb*duration,W+nl*ql) == 80_015_000
    for name in ['actual_capture_before_every_certificate','original_market_risks_action_and_fees_preserved','business_unused_exposure_refund_already_permitted']:
        assert d[name] is True
    for name in ['original_target_reached_by_this_lower_bound','original_target_proved_impossible','within_learning_window_capture_time_observed','practical_scale_or_optimal_policy_claim','published_papers_overwritten']:
        assert d[name] is False
    assert F(d['original_target']) == F(p['project_NPV_target']) == F(1,10)
    out = dict(source_sha256=sha(path),verification_code_sha256=sha(Path(__file__)),
        rational_interval_engine_sha256=sha(engine),original_parameters_and_source_hashes_checked=True,
        strict_bad_mean_refinement_rearrangement_checked=True,
        completed_binary_window_first_hit_discounted_geometric_book_checked=True,
        early_maximal_integer_window_cap_checked=True,
        business_refund_and_all_path_constant_reserve_checked=True,
        all_good_bad_error_composition_checked=True,late_Chernoff_tails_checked=2,
        finite_all_path_learning_budgets_and_branch_resources_checked=True,
        whole_value_and_seed_finance_cash_rate_rows_checked=len(rates),
        all_certificates_have_actual_paid_capture=True,
        high_rate_actual_capture_value_strict_lower='.09999105',
        cash_rate_fraction_sufficient_for_this_lower=str(rates[0]),
        high_rate_value_upper_to_lower_ratio_strictly_below='1.00370',
        original_point_one_target_reached=False,original_point_one_target_excluded=False,
        all_policy_optimality_or_practical_scale_claimed=False,
        general_analytic_arguments_formally_verified=False,published_papers_overwritten=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(dict(rate_rows=len(rates),strict_value_lower=out['high_rate_actual_capture_value_strict_lower'],actual_capture=True,original_target_reached=False)))

if __name__ == '__main__':
    main()
