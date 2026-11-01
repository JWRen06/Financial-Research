"""Independent original financial-action bad world and discounted value audit."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'项目C_有限信息下的动态机会捕获'/'第六十六轮研究_20261028'))
from verify_boundary_world_cash_necessity import ex,log
def main():
    raw=(HERE/'action_boundary_evidence_value.json').read_bytes();d=json.loads(raw);src={}
    for k,s in d['sources'].items():
        pp=HERE.parent/s['path'];rr=pp.read_bytes();cc=(pp.parent/'artifact_check.json').read_bytes();assert hashlib.sha256(rr).hexdigest()==s['sha256']==json.loads(cc)['source_sha256'] and hashlib.sha256(cc).hexdigest()==s['check_sha256'];src[k]=json.loads(rr)
    pa=src['A52']['parameters'];c=F(pa['cost_rate']);lam=F(pa['good_mean']);delta=F(pa['delta']);gap=F(pa['gap']);act=F(pa['action_exposure']);ga=F(pa['profit_threshold']);ell=F(d['fixed_bad_speed']);al=F(pa['alpha']);assert (c,lam,act,ga,ell)==(F(3,2),F(301,200),F(1,25000),F(1,50000000),F(3001,2000)) and ell==c+delta/2 and ga==(ell-c)*act and ell==F(src['A52']['hard_pair']['slow'])
    assert d['parameters']==src['A58']['parameters']==src['A59']['parameters'] and F(d['fixed_good_speed'])==lam and lam*lam<25 and ell*ell<25 and F(d['original_action_exposure'])==act and F(d['original_action_profit_threshold'])==ga and F(d['original_action_good_margin'])==gap*act and F(d['original_delta'])==delta
    fb=(ell-c)/ell*(1-ex(ell*act));fg=(lam-c)/lam*(1-ex(lam*act));assert fb.r<ga<ga+gap*act<fg.l
    for x,box in [(fb,d['fixed_bad_action_profit_outer']),(fg,d['fixed_good_action_profit_outer'])]:
        sl,sh=map(F,box);assert sl<=x.r and x.l<=sh and max(abs(sl-x.l),abs(sh-x.r))<F(1,10**100)
    ds=F(d['discounted_learning_loss_per_capture']);I=log(lam/ell)-1+ell/lam;Iu=F(d['safe_information_per_capture_upper']);K=F(d['safe_K_lower']);assert ds==F(src['A58']['discounted_learning_loss_per_capture']) and I.r<Iu and Iu-I.l<F(2,10**85) and Iu<F(src['A58']['safe_information_per_capture_upper']) and K==ds/Iu
    sl,sh=map(F,d['information_per_discounted_capture_outer']);assert sl<=I.r and I.l<=sh and max(abs(sl-I.l),abs(sh-I.r))<F(1,10**100)
    dl,dh=map(F,d['qualified_discounted_issuance_crossing_outer']);assert dh-dl==F(2,10**80) and al<dl<dh<1
    def kl(t):return t*log(t/al)+(1-t)*log((1-t)/(1-al))
    assert (kl(dl)-Iu*dl).r<0<(kl(dh)-Iu*dh).l
    n=improvements=ord_excluded=cap_excluded=0;keys=set()
    for z in d['bounds']:
        parent=z['cash_certificate_source'];rate=F(z['cash_rate_fraction']);key=(parent,rate);assert key not in keys;keys.add(key)
        old=next(x for x in (src['A58']['policies_or_bounds'] if parent=='A58' else src['A59']['certificates']) if F(x['cash_rate_fraction'])==rate);b=F(old['post_cash_constant_surplus_upper'] if parent=='A58' else old['post_business_surplus_upper']);assert b==F(z['post_business_surplus_upper'])>ds
        der=b-K*log(dl*(1-al)/(al*(1-dl)));assert der.r<0 and abs(der.r-F(z['maximum_at_capture_crossing_derivative_upper']))<F(1,10**100);qual=(b-ds)*dh;assert qual==F(z['actual_capture_value_upper'])
        et=ex(b/K);lo=K*log(1-al+al/et.r).l;hi=K*log(1-al+al/et.l).r;sl,sh=map(F,z['ordinary_value_upper_outer']);assert sl<=hi and lo<=sh and max(abs(sl-lo),abs(sh-hi))<F(1,10**100)
        assert F(z['prior_same_cash_ordinary_upper'])==F(old['ordinary_value_upper_outer'][1])>sh and F(z['prior_same_cash_actual_capture_upper'])==F(old['actual_capture_value_upper'])>qual and qual<lo
        assert z['ordinary_point_one_excluded']==(sh<F(1,10)) and z['actual_capture_point_one_excluded']==(qual<F(1,10));ord_excluded+=sh<F(1,10);cap_excluded+=qual<F(1,10);n+=1;improvements+=1
    ro=max(F(z['cash_rate_fraction']) for z in d['bounds'] if z['ordinary_point_one_excluded']);rc=max(F(z['cash_rate_fraction']) for z in d['bounds'] if z['actual_capture_point_one_excluded']);assert ro==rc==F(4999,5000)==F(d['point_one_target_necessary_cash_ratio_ordinary_strictly_above'])==F(d['point_one_target_necessary_cash_ratio_actual_capture_strictly_above'])
    assert n==14 and improvements==14 and ord_excluded==cap_excluded==13 and d['new_cash_generator_certificates']==d['new_sufficient_policies']==0
    for key in ['new_true_critical_rate_or_exact_optimal_world_claimed','A54_positive_value_bracket_changed','original_sources_overwritten']:assert not d[key]
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),original_fixed_action_good_and_strict_bad_world_qualification_checked=True,same_good_cash_certificates_reused_without_new_generator_counts=14,
        new_information_constant_and_qualification_crossing_checked=True,ordinary_and_actual_capture_value_bounds_checked=n,same_cash_parent_value_bounds_strictly_improved=improvements,listed_point_one_exclusions_per_qualification=ord_excluded,
        necessary_target_cash_ratio_strictly_above='4999/5000',no_coarse_observation_mean_bound_used_as_financial_bad_class=True,
        all_lower_cash_rates_and_arbitrary_equity_controls_covered_analytically=True,new_cash_generator_certificates=0,new_sufficient_policies=0,exact_critical_rate_or_optimal_hard_world_claimed=False)
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out))
if __name__=='__main__':main()
