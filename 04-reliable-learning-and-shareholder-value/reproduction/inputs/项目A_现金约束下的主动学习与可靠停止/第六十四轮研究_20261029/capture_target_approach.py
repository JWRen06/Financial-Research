"""Finite capture-qualified policy near the unchanged 0.1 target.
Decimal proposes enclosures only; the separate verifier uses rational intervals.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext, ROUND_FLOOR
import hashlib
import json

HERE = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dec(q):
    q = F(q)
    return D(q.numerator)/D(q.denominator)

def box(q):
    q = F(q)
    e = F(1, 10**65)
    return [str(q-e), str(q+e)]

def kl(t, p):
    return t*(t/p).ln()+(1-t)*((1-t)/(1-p)).ln()

def main():
    sources = {}
    loaded = {}
    for key, name in [('A62', '第六十二轮研究_20261029/uniform_value_ceiling.json'),
                      ('A63', '第六十三轮研究_20261029/ordinary_target_reserve.json')]:
        path = HERE.parent/name
        check = path.parent/'artifact_check.json'
        assert sha(path) == json.loads(check.read_text('utf8'))['source_sha256']
        sources[key] = dict(path=name, sha256=sha(path), check_sha256=sha(check))
        loaded[key] = json.loads(path.read_text('utf8'))
    parameters = loaded['A63']['parameters']
    control = dict(constant_business_reserve='603', business_exposure='73/100000000',
        business_windows=10_000_000_000, early_learning_exposure='1/1000000000',
        early_bad_risk_allocation='49999999/1000000000',
        late_bank_seed='1/100000000', late_learning_absolute_start='15000',
        late_learning_exposure='1/6250', late_learning_windows=500_000_000_000,
        late_accept_success_count=120_160_000,
        cash_rate_fractions=['9999999/10000000', '99999999/100000000', '1'])
    cf, af, mf, gaf = [F(parameters[k]) for k in
                      ['cost_rate','action_exposure','moment_bound','action_profit_threshold']]
    mbad_fraction = (cf+gaf/af+mf*af/2+cf*mf*af*af/6)/(1+cf*af/2)
    with localcontext() as ctx:
        ctx.prec = 160
        c, lam, m2, rho, fee, delay, A, gA, s = [dec(parameters[k]) for k in
            ['cost_rate', 'good_mean', 'moment_bound', 'valuation_rate',
             'business_fixed_access_fee', 'business_access_time', 'action_exposure',
             'action_profit_threshold', 'learning_success_fee']]
        # F_nu(A) >= (m-c)A-M2*A^2/2+c*m*A^2/2-c*M2*A^3/6.
        mbad = dec(mbad_fraction)
        q = dec(control['early_learning_exposure'])
        p0 = mbad*q
        p1 = lam*q-m2*q*q/2
        risk = dec(control['early_bad_risk_allocation'])
        ne = int(((1-risk).ln()/(1-p0).ln()).to_integral_value(rounding=ROUND_FLOOR))
        control['early_learning_windows'] = ne
        T = ne*q
        early_bad = 1+(ne*(1-p0).ln()).exp()*(-1)
        early_bad_next = 1+((ne+1)*(1-p0).ln()).exp()*(-1)
        early_good = 1+(ne*(1-p1).ln()).exp()*(-1)
        vq = (-rho*q).exp()
        survival = vq*(1-p1)
        S = (1+(ne*survival.ln()).exp()*(-1))/(1-survival)
        disc_hit = vq*p1*S
        disc_cost = c*q*S
        a = dec(control['business_exposure'])
        d = a+delay
        k = c*a+fee
        pg = lam*a-m2*a*a/2
        operating_refund = c*(lam*a*a/2-m2*a*a*a/6)
        X = dec(control['constant_business_reserve'])
        nb = control['business_windows']
        v = (-rho*d).exp()
        seed = dec(control['late_bank_seed'])
        W = dec(control['late_learning_absolute_start'])
        ql = dec(control['late_learning_exposure'])
        nl = control['late_learning_windows']
        t = dec(F(control['late_accept_success_count'], nl))
        late_p0 = (c+dec(parameters['delta']))*ql
        late_p1 = lam*ql-m2*ql*ql/2
        bad_tail = (-nl*kl(t, late_p0)).exp()
        good_tail = (-nl*kl(t, late_p1)).exp()
        whole_bad = early_bad+(1-early_bad)*bad_tail
        whole_good = (1-early_good)*good_tail
        budget = nl*c*ql+c*A
        L = c*T
        rows = []
        for fraction in control['cash_rate_fractions']:
            rate = rho*dec(fraction)
            y = (X-k)*(rate*d).exp()-X
            value = (v*(y+pg+operating_refund)/(1-v)-X)*(1-(-rho*d*nb).exp())
            spread = (rho-rate)*(X+L)*(rate*T).exp()*T
            lower = disc_hit*(value+1-s)-disc_cost-spread-seed
            seed_need = budget*(-rate*W).exp()
            assert y > 0 and seed_need < seed
            rows.append(dict(cash_rate_fraction=fraction,
                failure_window_surplus_outer=box(y), business_finite_NPV_outer=box(value),
                early_cash_spread_loss_upper_outer=box(spread),
                whole_good_NPV_lower_outer=box(lower),
                late_discounted_seed_needed_outer=box(seed_need),
                inherited_actual_capture_upper_to_policy_lower_ratio_outer=box(
                    dec(loaded['A62']['actual_capture_uniform_value_upper'])/lower)))
        assert whole_bad < dec(parameters['alpha']) and whole_good < dec(parameters['beta'])
        assert early_bad < risk < early_bad_next
        assert F(rows[0]['whole_good_NPV_lower_outer'][0]) > F('.09999105')
        out = dict(sources=sources, constructor_sha256=sha(Path(__file__)), parameters=parameters,
            control=control, strict_bad_mean_upper=str(mbad_fraction),
            early_bad_success_probability_upper=str(mbad_fraction*F(control['early_learning_exposure'])),
            early_good_success_probability_lower=str(F(p1)), early_maximum_duration=str(F(T)),
            early_worst_case_invoice_budget=str(F(L)),
            early_bad_certificate_probability_outer=box(early_bad),
            early_bad_probability_if_one_more_window_outer=box(early_bad_next),
            early_good_certificate_probability_lower_outer=box(early_good),
            early_discounted_hit_probability_lower_outer=box(disc_hit),
            early_discounted_invoice_cost_upper_outer=box(disc_cost),
            business_window_duration=str(F(d)), business_prepayment=str(F(k)),
            business_good_success_probability_lower=str(F(pg)),
            business_expected_unused_exposure_refund_lower=str(F(operating_refund)),
            business_calendar=str(F(nb*d)), late_bad_success_probability_upper=str(F(late_p0)),
            late_good_success_probability_lower=str(F(late_p1)),
            late_success_frequency_threshold=str(F(t)), late_bad_error_chernoff_outer=box(bad_tail),
            late_good_error_chernoff_outer=box(good_tail),
            whole_bad_certificate_probability_upper_outer=box(whole_bad),
            whole_good_rejection_probability_upper_outer=box(whole_good),
            late_worst_case_prepayment_budget=str(F(budget)), whole_initial_equity=str(F(X+L+seed)),
            whole_opportunity_capacity=ne+max(nb,nl),
            whole_final_calendar=str(F(max(T+nb*d, W+nl*ql))),
            rate_menu=rows,
            inherited_actual_capture_uniform_upper=loaded['A62']['actual_capture_uniform_value_upper'],
            original_target='1/10', original_target_reached_by_this_lower_bound=False,
            original_target_proved_impossible=False, actual_capture_before_every_certificate=True,
            original_market_risks_action_and_fees_preserved=True,
            within_learning_window_capture_time_observed=False,
            business_unused_exposure_refund_already_permitted=True,
            practical_scale_or_optimal_policy_claim=False, published_papers_overwritten=False)
        (HERE/'capture_target_approach.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
        print(json.dumps(dict(early_windows=ne, early_duration=float(T), capital=float(X+L+seed),
            bad_error=float(whole_bad), good_rejection=float(whole_good),
            rate_values=[(r['cash_rate_fraction'],float(F(r['whole_good_NPV_lower_outer'][0]))) for r in rows],
            original_target_reached=False),ensure_ascii=False))

if __name__ == '__main__':
    main()

