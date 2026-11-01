"""Independent C66 log/exp, root exclusion and same-quality comparisons."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'项目A_现金约束下的主动学习与可靠停止'/'第五十二轮研究_20261028'))
from verify_joint_time_capacity import J,ex,log
def kl(p,q):return p*log(p/q)+(1-p)*log((1-p)/(1-q))
def near(x,bb,tol):
    a,b=map(F,bb);assert x.l<=b and a<=x.r and max(abs(x.l-a),abs(x.r-b))<tol
def main():
    raw=(HERE/'boundary_world_cash_necessity.json').read_bytes();d=json.loads(raw);src={}
    for k,s in d['sources'].items():
        pp=HERE.parent/s['path'];rr=pp.read_bytes();cr=(pp.parent/'artifact_check.json').read_bytes();assert hashlib.sha256(rr).hexdigest()==s['sha256']==json.loads(cr)['source_sha256'] and hashlib.sha256(cr).hexdigest()==s['check_sha256'];src[k]=json.loads(rr)
    al=F(d['alpha']);be=F(d['beta']);assert al==F(src['C63']['alpha']) and be==F(src['C63']['beta']);b0l,b0h=map(F,d['free_reference_capital_outer']);assert d['free_reference_capital_outer']==src['C63']['free_reference_capital_outer']
    A=F(198,445);f=F(1,1000);c=F(1,5);g=F(1,10);d0=2*(A-f-g)/c;assert d0==F(d['supporting_tangent_slope_at_reference']) and g-d0*b0h>0
    def value(B):return A*(1-ex(2*(B-f)/c))-f
    def slope(B):return 2*A/c*ex(2*(B-f)/c)
    constants={}
    for z in d['quality_constants']:
        q=F(z['certificate_quality']);D=kl(F(1,89),q/F(89,100));Du=F(z['forward_KL_upper_rational']);assert Du>D.r and F(1,100)<q<=F(2,25)
        near(D,z['per_H_forward_KL_outer'],F(1,10**75));coef=(1+d0)/d0*log(1/al)/D;near(coef,z['asymptotic_capital_excess_per_p_necessary_lower_outer'],F(1,10**75));constants[q]=Du
    n=oldchecked=0;qualities={};pricegroups=set()
    for z in d['policies']:
        old=src['C63']['policies'][z['C63_policy_index']];q=F(z['certificate_quality']);mu=F(z['price_multiplier']);p=F(z['label_price']);Du=constants[q]
        assert q==F(old['certificate_quality']) and mu==F(old['price_multiplier']) and p==F(old['label_price']) and F(z['boundary_bad_H_mean'])==q/F(89,100)
        bl,bh=map(F,z['safe_necessary_capital_root_outer']);assert bh-bl==F(2,10**80) and b0h<bl<bh<F(d['local_capital_domain_upper'])==F(28,1000)
        def residual(B):
            vv=value(B);der=slope(B);assert vv.l>g and (vv-der*B).l>0
            aL=max(1-be,g/vv.r);aH=max(1-be,g/vv.l);assert al<aL<=aH<1
            kL=kl(aL,al).l;kH=kl(aH,al).r
            return J(vv.l-g-(1+der.r)*p/Du*kH,vv.r-g-(1+der.l)*p/Du*kL),[aL,aH]
        fl,aa=residual(bl);fr,bb=residual(bh);assert fl.r<0<fr.l
        tolerance=min(-fl.r,fr.l)/100;near(fl,z['residual_signs_outer'][0],tolerance);near(fr,z['residual_signs_outer'][1],tolerance)
        # Source acceptance boxes use the constructor's exponential bounds;
        # independent acceptance bounds must overlap them at high precision.
        near(J(*aa),z['acceptance_probability_lower_outer_at_root_lower'],F(1,10**75));near(J(*bb),z['acceptance_probability_lower_outer_at_root_upper'],F(1,10**75))
        assert F(z['complete_capital_necessary_lower'])==bl and F(z['capital_excess_necessary_lower'])==bl-b0h
        upper=F(old['whole_initial_capital']);assert upper==F(z['same_quality_C63_capital_upper'])>bh and upper/bl==F(z['feasible_to_necessary_total_capital_ratio']) and (upper-b0l)/(bl-b0h)==F(z['feasible_to_necessary_excess_ratio'])
        prior=next((s for s in src['C62']['policies'] if F(s['price_multiplier'])==mu and s['computation_label_cap']==256),None)
        if prior:
            gap=F(prior['capital_excess_necessary_lower']);assert gap==F(z['C62_necessary_excess_lower'])<bl-b0h and z['stronger_than_C62_when_available'] is True;oldchecked+=1
        else:assert z['C62_necessary_excess_lower'] is None and z['stronger_than_C62_when_available'] is None
        qualities[q]=qualities.get(q,0)+1;pricegroups.add(mu);n+=1
    assert n==21 and oldchecked==6 and len(pricegroups)==7 and len(qualities)==3 and all(v==7 for v in qualities.values())
    assert d['new_statistical_policies']==0 and d['new_necessary_capital_root_brackets']==21 and d['original_services_and_full_financial_target_kept'] and d['fixed_quality_not_silently_relaxed']
    for key in ['boundary_world_assumed_subject_to_bad_class_error_constraint','coefficient_claimed_optimal_or_achieved','original_results_overwritten']:assert d[key] is False
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_rational_math_digits=110,new_fixed_quality_necessary_capital_roots_checked=n,
        certificate_quality_groups_checked=3,inherited_C63_policy_comparisons_checked=21,inherited_C62_necessary_bounds_improved=oldchecked,
        asymptotic_cash_shadow_price_coefficients_checked=3,positive_tangent_intercept_checked=True,
        boundary_limit_uses_only_strict_bad_worlds_proved_analytically=True,all_adaptive_service_controls_and_full_good_profit_preserved=True,
        original_all_market_alpha_beta_kept=True,new_sufficient_policies_claimed=0,
        finite_optimum_or_achieved_asymptotic_coefficient_claimed=False,certificate_quality_relaxed_in_comparison=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out))
if __name__=='__main__':main()
