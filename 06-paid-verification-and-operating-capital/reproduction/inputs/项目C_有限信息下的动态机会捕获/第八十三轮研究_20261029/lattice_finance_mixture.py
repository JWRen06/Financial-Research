"""Rational integer-score controllers and original full-market finance.
Decimal root proposals are independently checked by verify_lattice_finance_mixture.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext,ROUND_CEILING
import hashlib,json,sys
sys.set_int_max_str_digits(100_000)
HERE=Path(__file__).resolve().parent
VERT=[(F(0),F(0)),(F(1,89),F(0)),(F(0),F(1,11))]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(q):q=F(q);return D(q.numerator)/D(q.denominator)
def box(q):q=F(q);e=F(1,10**65);return [str(q-e),str(q+e)]
def upper(q,den):return F((F(q)*den).__ceil__(),den)

def envelope(H,L,first,last,power,quality,left,right):
    # Factors (1,first) for H and (1,last) for L, times a common prefactor.
    b=first-1;d=1+(last-1)*quality/F(11,100);e=-(last-1)*F(89,11)
    critical=-(H*b*d+L*e)/((H+L)*b*e);points=[left,right]
    if left<critical<right:points.append(critical)
    values=[power*(1+b*x)**H*(d+e*x)**L for x in points]
    return dict(quality=str(quality),H_domain=list(map(str,[left,right])),critical_H_error=str(critical),
        tested_H_errors=list(map(str,points)),exact_values=list(map(str,values)),maximum=str(max(values)))

def main():
    sources={};loaded={}
    for key,name in [('C79','第七十九轮研究_20261029/finance_typed_cohort.json'),('C82','第八十二轮研究_20261029/monotone_bias_endpoints.json')]:
        p=HERE.parent/name;check=p.parent/'artifact_check.json'
        assert sha(p)==json.loads(check.read_text('utf8'))['source_sha256']
        sources[key]=dict(path=name,sha256=sha(p),check_sha256=sha(check));loaded[key]=json.loads(p.read_text('utf8'))
    old=loaded['C79'];pilot=json.loads((HERE/'lattice_mixture_pilot.json').read_text('utf8'))
    H,L=22,2;cost=F(old['procurement_fee_envelope']['expected_cohort_net_invoice_upper_over_p'])
    alpha=F(old['alpha']);beta=F(old['beta']);tau=F(9,1000)
    weights=[F(878893,10**6),F(121107,10**6)];controllers=[];policies=[]
    with localcontext() as ctx:
        ctx.prec=160
        for chosen,weight in zip(pilot['selected'],weights):
            r=chosen['candidate'];assert r['quotas']==[H,L]
            m=r['positive_score_step'];kh=r['H_error_score_penalty'];kl=r['L_error_score_penalty'];z=F(str(r['evidence_base']));K=r['score_threshold']
            tail=F(format((-D(str(r['tail_delta']))).exp(),'.15f'))
            bad=envelope(H,L,z**(-kh),z**(-kl),z**m,F(41,1000),F(0),F(41,890))
            moment=envelope(H,L,tail**(-kh),tail**(-kl),tail**m,F(19,1000),tau,tau+F(1,89))
            maximum=F(moment['maximum']);ds=[m-kh*H*(x+tau)-kl*L*(y+tau) for x,y in VERT]
            assert F(bad['maximum'])<1 and maximum<1 and min(ds)>0
            assert z**(K-1)<1/alpha<=z**K
            controllers.append(dict(weight=str(weight),positive_score_step=m,H_error_penalty=kh,L_error_penalty=kl,
                evidence_base=str(z),score_threshold=K,score_at_stopping_upper=K+m-1,
                bad_product_envelope=bad,good_tail_base=str(tail),good_tail_moment_envelope=moment,
                good_score_drift_at_vertices=list(map(str,ds)),expected_cohort_at_vertices=list(map(str,[F(K+m-1)/d for d in ds]))))
        f=D('.001');c=D('.2');g=D('.1');b0l,b0h=map(F,old['free_reference_capital_outer'])
        def phi(B,x,y):
            t=(B-f)/c
            return D('.45')*(1-dec(x))*(1-(-2*t).exp())-D('9.5')*dec(y)*(1-(-t/100).exp())-f
        for previous in old['policies']:
            p=F(previous['label_price']);branches=[]
            for ctrl in controllers:
                rate=-dec(ctrl['good_tail_moment_envelope']['maximum']).ln();delta=-dec(ctrl['good_tail_base']).ln();K=ctrl['score_threshold']
                n=int(((dec(1/p).ln()+(K-1)*delta)/rate).to_integral_value(rounding=ROUND_CEILING))
                bill=64*n*p;means=[cost*p*F(x) for x in ctrl['expected_cohort_at_vertices']]
                branches.append(dict(cohort_cap=n,label_cap=24*n,maximum_invoice=str(bill),
                    cap_tail_exponent_outer=box(n*rate-(K-1)*delta),
                    vertex_expected_invoice_uppers=list(map(str,means))))
            def lip(B,bill):
                k=dec(bill)
                return D('.45')*((-2*(B-k-f)/c).exp()-(-2*(B-f)/c).exp())/k
            def values(B):
                return [(1-dec(p))*phi(B,x,y)-sum(dec(w)*(1+lip(B,br['maximum_invoice']))*dec(br['vertex_expected_invoice_uppers'][j])
                    for w,br in zip(weights,branches)) for j,(x,y) in enumerate(VERT)]
            def pure_values(B,i):
                br=branches[i]
                return [(1-dec(p))*phi(B,x,y)-(1+lip(B,br['maximum_invoice']))*dec(br['vertex_expected_invoice_uppers'][j]) for j,(x,y) in enumerate(VERT)]
            lo=dec(b0h)-D('.001');hi=dec(previous['whole_initial_capital'])+D('.001')
            assert min(values(lo))<g<min(values(hi))
            for _ in range(500):
                mid=(lo+hi)/2
                if min(values(mid))>g:hi=mid
                else:lo=mid
            proposal=F((lo+hi)/2);bl=proposal-F(1,10**80);bh=proposal+F(1,10**80);B=upper(bh,10**90)
            assert min(values(dec(B)))>g and B<F(previous['whole_initial_capital'])
            a69=F(previous['all_r_C69_necessary_lower']);a77=F(previous['r_at_least_one_C77_necessary_lower'])
            maxbill=max(F(b['maximum_invoice']) for b in branches)
            pure_roots=[]
            if F(previous['price_multiplier'])==1:
                for i in range(2):
                    l=dec(b0h)-D('.001');r=dec(previous['whole_initial_capital'])+D('.001')
                    for _ in range(500):
                        t=(l+r)/2
                        if min(pure_values(t,i))>g:r=t
                        else:l=t
                    q=F((l+r)/2);pure_roots.append([str(q-F(1,10**80)),str(q+F(1,10**80))])
            policies.append(dict(price_multiplier=previous['price_multiplier'],label_price=str(p),branches=branches,
                capital_root_outer=list(map(str,[bl,bh])),whole_initial_capital=str(B),
                mixed_full_good_vertex_profit_outer=[box(v) for v in values(dec(B))],
                separate_branch_vertex_profit_outer_at_mixed_capital=[[box(v) for v in pure_values(dec(B),i)] for i in range(2)],
                original_price_pure_certificate_capital_roots=pure_roots,
                maximum_label_cap=max(b['label_cap'] for b in branches),maximum_actual_invoice=str(maxbill),
                maximum_natural_stage_prepayment=str(41*p),all_path_exposure_interval=[str((B-maxbill-F(1,1000))/F(1,5)),str((B-F(1,1000))/F(1,5))],
                previous_C79_capital=previous['whole_initial_capital'],previous_C79_label_cap=previous['deterministic_label_cap'],previous_C79_maximum_invoice=previous['maximum_actual_total_invoice'],
                capital_reduction=str(F(previous['whole_initial_capital'])-B),
                all_r_necessary_lower=str(a69),r_at_least_one_necessary_lower=str(a77),
                all_r_excess_ratio=str((B-b0l)/(a69-b0h)),r_at_least_one_excess_ratio=str((B-b0l)/(a77-b0h))))
            print(json.dumps(dict(mu=previous['price_multiplier'],capital=float(B),labels=policies[-1]['maximum_label_cap'],invoice=float(maxbill),gap_r_ge1=float(F(policies[-1]['r_at_least_one_excess_ratio'])))),flush=True)
        out=dict(sources=sources,constructor_sha256=sha(Path(__file__)),proposal_source=dict(file='lattice_mixture_pilot.json',sha256=sha(HERE/'lattice_mixture_pilot.json')),
            alpha=str(alpha),beta=str(beta),certificate_quality=old['certificate_quality'],true_good_quality=old['true_good_quality'],service_accuracy=str(tau),
            cohort_quotas={'H':H,'L':L},procurement_fee_envelope=old['procurement_fee_envelope'],free_reference_capital_outer=old['free_reference_capital_outer'],
            controllers=controllers,policies=policies,integer_score_and_private_initial_coin=True,
            all_good_financial_triangle_covered_by_concave_mixture=True,
            market_fixed_before_private_coin_and_all_observations=True,
            new_full_market_mixture_policies=len(policies),generic_overshoot_or_randomization_method_novelty_claimed=False,
            finite_candidate_LP_or_all_policy_optimum_claimed=False,old_artifacts_or_published_papers_overwritten=False)
        (HERE/'lattice_finance_mixture.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')

if __name__=='__main__':main()
