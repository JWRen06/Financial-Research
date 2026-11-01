"""Independent rational full-market audit; no constructor import or Decimal."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,sys
sys.set_int_max_str_digits(100_000)
HERE=Path(__file__).resolve().parent
ENGINE=HERE.parent.parent/'项目A_现金约束下的主动学习与可靠停止/第五十二轮研究_20261028/verify_joint_time_capacity.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(ENGINE)=='af87115e344b638fea22c0a0cf65fd66bfedca99ccff6feb2df34cfb7b8b0c51'
sys.path.insert(0,str(ENGINE.parent))
from verify_joint_time_capacity import J,ex,log

def contained(x,bs):
    l,r=map(F,bs);assert l<x.l<=x.r<r,(float(x.l-l),float(r-x.r));assert r-l==F(2,10**65)

def bounded_log(x):
    # Round enormous exact moment products first, keeping a 100-digit outward
    # enclosure; monotonicity gives a rational interval for the genuine log.
    den=10**100;l=F((x*den).__floor__(),den);r=F((x*den).__ceil__(),den)
    assert 0<l<=x<=r
    return J(log(l).l,log(r).r)

def main():
    path=HERE/'lattice_finance_mixture.json';d=json.loads(path.read_text('utf8'))
    assert d['constructor_sha256']==sha(HERE/'lattice_finance_mixture.py')
    src={}
    for key,row in d['sources'].items():
        p=HERE.parent/row['path'];ck=p.parent/'artifact_check.json'
        assert sha(p)==row['sha256']==json.loads(ck.read_text('utf8'))['source_sha256'] and sha(ck)==row['check_sha256']
        src[key]=json.loads(p.read_text('utf8'))
    assert sha(HERE/d['proposal_source']['file'])==d['proposal_source']['sha256']
    old=src['C79'];assert all(d[k]==old[k] for k in ['alpha','beta','certificate_quality','true_good_quality'])
    alpha=F(d['alpha']);beta=F(d['beta']);tau=F(d['service_accuracy']);assert tau==F(9,1000)
    assert d['cohort_quotas']==old['cohort_quotas']=={'H':22,'L':2}
    H,L=22,2;vertices=[(F(0),F(0)),(F(1,89),F(0)),(F(0),F(1,11))]
    fee=d['procurement_fee_envelope'];assert fee==old['procurement_fee_envelope']
    den=fee['state_denominator'];states=[[F(n,den) for n in row] for row in fee['state_numerators']]
    interior=boundary=0
    for h in range(H+1):
        for l in range(L+1):
            if h==H:assert states[h][l]==20*(L-l);boundary+=1
            elif l==L:assert states[h][l]==2*(H-h);boundary+=1
            else:
                for s in [F(1,10),F(3,20)]:
                    assert states[h][l]>=(1-s)*(1+states[h+1][l])+s*(10+states[h][l+1]);interior+=1
    cost=F(fee['expected_cohort_net_invoice_upper_over_p']);assert cost==states[0][0]
    assert fee['maximum_cohort_net_invoice_over_p']==64==max(H+20*L,2*H+10*L)
    assert fee['maximum_natural_stage_prepayment_over_p']==H+10*L-1==41<64
    def envelope(first,last,mult,q,left,right,stored):
        b=first-1;slope=-(last-1)*F(89,11);intercept=1+(last-1)*q/F(11,100)
        critical=(H*b*intercept+L*slope)/(-(H+L)*b*slope)
        nodes=[left,right]
        if left<critical<right:nodes.append(critical)
        values=[mult*(1+b*x)**H*(intercept+slope*x)**L for x in nodes]
        assert list(map(F,stored['H_domain']))==[left,right] and F(stored['quality'])==q and F(stored['critical_H_error'])==critical
        assert list(map(F,stored['tested_H_errors']))==nodes and list(map(F,stored['exact_values']))==values and F(stored['maximum'])==max(values)
        assert all(1+b*x>0 and intercept+slope*x>0 for x in nodes)
        return max(values),len(nodes)
    weights=[];constants=[];bad_nodes=moment_nodes=0
    expected=[(11,15,23,F('1.04754781528508'),203,F(878893,10**6)),(3,4,6,F('1.17443039813211'),59,F(121107,10**6))]
    for ctrl,expected_row in zip(d['controllers'],expected,strict=True):
        m=ctrl['positive_score_step'];kh=ctrl['H_error_penalty'];kl=ctrl['L_error_penalty'];z=F(ctrl['evidence_base']);K=ctrl['score_threshold'];w=F(ctrl['weight'])
        assert (m,kh,kl,z,K,w)==expected_row and z>1 and K>0
        assert z**(K-1)<1/alpha<=z**K and ctrl['score_at_stopping_upper']==K+m-1
        bad,nodes=envelope(z**(-kh),z**(-kl),z**m,F(41,1000),F(0),F(41,890),ctrl['bad_product_envelope']);assert bad<1;bad_nodes+=nodes
        s=F(ctrl['good_tail_base']);assert 0<s<1
        moment,nodes=envelope(s**(-kh),s**(-kl),s**m,F(19,1000),tau,tau+F(1,89),ctrl['good_tail_moment_envelope']);assert moment<1;moment_nodes+=nodes
        drifts=[m-kh*H*(x+tau)-kl*L*(y+tau) for x,y in vertices]
        assert min(drifts)>0 and list(map(F,ctrl['good_score_drift_at_vertices']))==drifts
        expected_cohorts=[F(K+m-1)/x for x in drifts]
        assert list(map(F,ctrl['expected_cohort_at_vertices']))==expected_cohorts
        # Exact bounded-step overshoot on all integer states immediately below
        # K, with every feasible pair of H/L error counts in the next cohort.
        for previous in range(K-m,K):
            for eh in range(H+1):
                for el in range(L+1):
                    nxt=previous+m-kh*eh-kl*el
                    if nxt>=K:assert nxt<=K+m-1
        weights.append(w);constants.append((K,-bounded_log(moment),-bounded_log(s),expected_cohorts))
    assert sum(weights)==1 and all(w>0 for w in weights)
    b0l,b0h=map(F,d['free_reference_capital_outer']);assert d['free_reference_capital_outer']==old['free_reference_capital_outer']
    f=F(1,1000);c=F(1,5);g=F(1,10)
    def phi(B,x,y):
        t=(B-f)/c
        return F(9,20)*(1-x)*(1-ex(2*t))-F(19,2)*y*(1-ex(t/100))-f
    def chord(B,bill):return F(9,20)*(ex(2*(B-bill-f)/c)-ex(2*(B-f)/c))/bill
    max69=max77=F(0);policy_count=0;branch_checks=0
    for row in d['policies']:
        oldrow=next(r for r in old['policies'] if r['price_multiplier']==row['price_multiplier']);p=F(row['label_price']);assert p==F(oldrow['label_price'])<beta
        logp=log(1/p);bills=[];fees=[]
        for br,(K,rate,delta,cohorts) in zip(row['branches'],constants,strict=True):
            n=br['cohort_cap'];bill=F(br['maximum_invoice']);exponent=n*rate-(K-1)*delta
            contained(exponent,br['cap_tail_exponent_outer']);assert exponent.l>logp.r
            assert ((n-1)*rate-(K-1)*delta).r<logp.l
            assert bill==64*n*p and br['label_cap']==24*n
            means=[cost*p*v for v in cohorts];assert list(map(F,br['vertex_expected_invoice_uppers']))==means
            assert max(means)<bill
            bills.append(bill);fees.append(means);branch_checks+=1
        def profits(B):
            return [(1-p)*phi(B,x,y)-sum((weights[i]*(1+chord(B,bills[i]))*fees[i][j] for i in range(2)),J(0)) for j,(x,y) in enumerate(vertices)]
        def pure(B,i):return [(1-p)*phi(B,x,y)-(1+chord(B,bills[i]))*fees[i][j] for j,(x,y) in enumerate(vertices)]
        bl,bh=map(F,row['capital_root_outer']);assert bh-bl==F(2,10**80) and min(x.r for x in profits(bl))<g<min(x.l for x in profits(bh))
        B=F(row['whole_initial_capital']);assert B==F((bh*10**90).__ceil__(),10**90)<F(oldrow['whole_initial_capital'])
        for value,stored in zip(profits(B),row['mixed_full_good_vertex_profit_outer']):
            contained(value,stored);assert value.l>g
        for i in range(2):
            for value,stored in zip(pure(B,i),row['separate_branch_vertex_profit_outer_at_mixed_capital'][i]):contained(value,stored)
        if F(row['price_multiplier'])==1:
            assert len(row['original_price_pure_certificate_capital_roots'])==2
            for i,bracket in enumerate(row['original_price_pure_certificate_capital_roots']):
                l,r=map(F,bracket);assert B<l<r and r-l==F(2,10**80) and min(v.r for v in pure(l,i))<g<min(v.l for v in pure(r,i))
            assert pure(B,0)[2].r<g<pure(B,0)[1].l and pure(B,1)[1].r<g<pure(B,1)[2].l
        else:assert not row['original_price_pure_certificate_capital_roots']
        assert F(row['maximum_actual_invoice'])==max(bills)<F(oldrow['maximum_actual_total_invoice'])==F(row['previous_C79_maximum_invoice'])
        assert row['maximum_label_cap']==max(br['label_cap'] for br in row['branches'])<oldrow['deterministic_label_cap']==row['previous_C79_label_cap']
        assert F(row['maximum_natural_stage_prepayment'])==41*p
        low,high=map(F,row['all_path_exposure_interval']);assert low==(B-max(bills)-f)/c and high==(B-f)/c and 0<low<high<F(3,20)
        assert all(phi(B-max(bills),x,y).l>0 for x,y in vertices)
        assert F(row['previous_C79_capital'])==F(oldrow['whole_initial_capital']) and F(row['capital_reduction'])==F(oldrow['whole_initial_capital'])-B>0
        n69=F(row['all_r_necessary_lower']);n77=F(row['r_at_least_one_necessary_lower']);assert n69==F(oldrow['all_r_C69_necessary_lower'])<n77==F(oldrow['r_at_least_one_C77_necessary_lower'])<B
        r69=(B-b0l)/(n69-b0h);r77=(B-b0l)/(n77-b0h)
        assert r69==F(row['all_r_excess_ratio']) and r77==F(row['r_at_least_one_excess_ratio'])
        max69=max(max69,r69);max77=max(max77,r77);policy_count+=1
    assert policy_count==d['new_full_market_mixture_policies']==7 and branch_checks==14
    first=d['policies'][0];assert F(first['whole_initial_capital'])<F('.026826085820') and first['maximum_label_cap']==13920 and F(first['maximum_actual_invoice'])==F('.0007424')
    assert max69<F('3.73') and max77<F('1.196')
    for key in ['integer_score_and_private_initial_coin','all_good_financial_triangle_covered_by_concave_mixture','market_fixed_before_private_coin_and_all_observations']:assert d[key] is True
    for key in ['generic_overshoot_or_randomization_method_novelty_claimed','finite_candidate_LP_or_all_policy_optimum_claimed','old_artifacts_or_published_papers_overwritten']:assert d[key] is False
    out=dict(source_sha256=sha(path),verification_code_sha256=sha(Path(__file__)),rational_interval_engine_sha256=sha(ENGINE),
        original_alpha_beta_quality_profit_and_service_contract_preserved=True,
        integer_score_controllers_checked=2,bad_boundary_endpoint_stationary_values_checked=bad_nodes,
        good_tail_endpoint_stationary_values_checked=moment_nodes,positive_exact_drift_vertices_checked=6,
        exact_threshold_and_integer_overshoot_checked=True,procurement_fee_affine_inequalities_checked=interior,procurement_fee_boundary_states_checked=boundary,
        finite_cap_tail_certificates_checked=branch_checks,whole_original_full_market_mixture_policies_checked=policy_count,
        full_good_finance_covered_by_concavity_not_grid=True,all_path_prepayments_and_complete_cash_profit_checked=True,
        fixed_market_private_coin_mixture_semantics_proved_analytically=True,
        original_price_two_pure_bound_roots_and_opposite_financial_vertex_deficits_checked=True,
        all_seven_capital_label_and_max_invoice_bounds_improve_C79=True,
        original_price_sufficient_capital_strict_upper='.026826085820',original_price_maximum_labels=13920,original_price_maximum_invoice='.0007424',
        all_r_excess_ratio_strict_upper='3.73',r_at_least_one_excess_ratio_strict_upper='1.196',
        max_all_r_excess_ratio_display=float(max69),max_r_at_least_one_excess_ratio_display=float(max77),
        all_resources_expected_invoice_or_optimal_policy_improvement_claimed=False,general_analytic_proofs_formally_verified=False,
        published_papers_overwritten=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(dict(policies=policy_count,branch_tail_checks=branch_checks,capital_upper=out['original_price_sufficient_capital_strict_upper'],labels=13920,all_r_gap=float(max69),conditional_gap=float(max77))))

if __name__=='__main__':main()
