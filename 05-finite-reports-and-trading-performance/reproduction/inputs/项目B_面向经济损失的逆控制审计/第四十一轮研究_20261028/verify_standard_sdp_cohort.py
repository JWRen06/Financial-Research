"""B41 independent per-piece caches and complete/partial source cover audit.
No proposal/solver implementation imported. Partial snapshots remain frozen.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第三十八轮研究_20261027'))
from verify_independent_direction_pilot import arithmetic,Q

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(path,d):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(path)

def verify_piece(d,sourcepart,band,baseline,cal,pieces,strips,reused_b40):
    assert d['beta_interval']==band['beta_interval'] and d['allowed_pole_indices']==band['allowed_pole_indices'] and d['root_gap_upper']==band['root_gap_upper']
    assert d['direction']==sourcepart['runs'][0]['direction'];a,b=map(F,d['direction']);assert a*a+b*b==1 and b>=0
    lo,hi=map(F,band['beta_interval']);e0=arithmetic.root(arithmetic.b0,F(1));delta=F(d['root_gap_upper']);whole=Q(arithmetic.root(hi,F(1)).l,arithmetic.root(lo,F(1)).r);assert (whole-e0).r<=delta
    axes=[z for z in band['facets'] if z['origin'] in ['nonnegative_perpendicular','inherited_B20_B27_axis']];facets=list(sourcepart['final_facets']);assert len(axes)==4 and len(d['methods'])==2
    bcells=pcells=rawcells=qualified=proposals=solved=near=errors=sdpbetter=0;statuses={};budgets={}
    shifts=list(map(str,[F(0),F(1,10**7),F(1,10**6),F(1,10**5),F(1,10000),F(3,10000),F(1,1000),F(3,1000),F(1,100),F(3,100),F(1,10)]))
    for z,method in zip(d['methods'],['LP','SDP'],strict=True):
        assert z['method']==method and z['reused_from_B40']==d['reused_B40_pair']
        if 'proposal' in z:
            r=z['proposal'];assert r['method']==method and r['scalar_support_samples_included'] is True and r['solver']=='Clarabel 0.11.1';proposals+=1
            statuses[method]=r['status'];solved+=r['status']=='Solved';near+=r['status']=='AlmostSolved'
            nbase=9+8*r['pole_samples'];assert r['variables']==nbase+(978 if method=='SDP' else 0)
            assert r['PSD_orders']==([19,16]*3 if method=='SDP' else []) and r['equalities']==(258 if method=='SDP' else 0)
        else:assert 'construction_error_type' in z;statuses[method]=z['construction_error_type']
        if 'construction_error_type' in z:errors+=1
        if d['reused_B40_pair']:
            old=next(x for x in reused_b40['methods'] if x['method']==method);assert z['certified']==old['certified'] and z['qualification_trials']==old['qualification_trials']
            if z['certified']:assert z['proof']==old['proof']
        else:assert z['qualification_finished'] is True
        if not z['certified']:
            assert 'proof' not in z and all(not t['passed'] for t in z['qualification_trials']);continue
        r=z['proposal'];pr=z['proof'];tr=z['qualification_trials'];assert [x['shift'] for x in tr]==shifts[:len(tr)] and tr[-1]['passed'] and all(not x['passed'] for x in tr[:-1])
        shift=F(pr['coefficient_shift']);cc=list(map(F,pr['coefficients']));assert shift==F(tr[-1]['shift']) and cc==[F(x)+shift for x in r['candidate_coefficients']] and min(cc[0],cc[4])>=0
        last=lo;assert len(pr['complete_beta_kernel_cover'])==16
        for j,row in enumerate(pr['complete_beta_kernel_cover']):
            ll,hh=map(F,row['beta_interval']);assert ll==last==lo+(hi-lo)*F(j,16) and hh==lo+(hi-lo)*F(j+1,16);last=hh
            es=[Q(arithmetic.root(hh,F(k)).l,arithmetic.root(ll,F(k)).r) for k in [1,2]];point=F(0);kp=row['kernel_proof']
            for cell in kp['cover']:
                l,r=map(F,cell['persistence']);assert l==point<r;point=r;polys,initial=arithmetic.polynomials(es,e0,cc,a,b,l,r,kp['method']);assert min(x.l for x in initial)>=0 and all(min(x.l for x in poly)>=0 for poly in polys);pcells+=1
            assert point==1;bcells+=1
        assert last==hi;M=F(pr['raw_budget_upper']);rawcells+=arithmetic.raw_check(cc,M,pr['raw_calibration_proof'],cal,pieces,strips,band['allowed_pole_indices'])
        facet=dict(origin='B41_'+method,halfspace=[str(a),str(b),str(M)]);facets.append(facet);_,radius=arithmetic.polygon(axes+[facet]);assert delta*delta*radius==F(z['axis_plus_one_facet_loss_upper'])
        old=sourcepart['runs'][0].get('proof');assert z['old_B39_same_direction_budget']==(None if old is None else old['raw_budget_upper']) and z['improved_old_same_direction_budget']==bool(old and M<F(old['raw_budget_upper']))
        budgets[method]=M;qualified+=1
    _,radius=arithmetic.polygon(facets);new=delta*delta*radius;assert d['final_facets']==facets and F(d['final_squared_radius'])==radius and F(d['new_local_loss_upper'])==new
    hybrid=min(baseline,new);assert F(d['prior_local_loss_upper'])==baseline and F(d['hybrid_local_loss_upper'])==hybrid and d['qualified_methods']==qualified
    return dict(qualified_methods=qualified,proposal_records=proposals,exact_solver_solved_flags=solved,near_solved_flags=near,construction_errors=errors,statuses=statuses,
        beta_subintervals=bcells,persistence_cells=pcells,raw_calibration_cells=rawcells,hybrid_local_loss_upper=str(hybrid),improved=hybrid<baseline,
        both_methods_qualified=('LP' in budgets and 'SDP' in budgets),SDP_smaller_certified_budget=('LP' in budgets and 'SDP' in budgets and budgets['SDP']<budgets['LP']),
        LP_smaller_certified_budget=('LP' in budgets and 'SDP' in budgets and budgets['LP']<budgets['SDP']))

def main():
    start=time.monotonic();partial='--partial' in sys.argv;path=HERE/('standard_sdp_cohort_checkpoint.json' if partial else 'standard_sdp_cohort.json');raw=path.read_bytes();d=json.loads(raw);data={}
    assert digest(HERE/'standard_sdp_cohort.py')==d['construction_code_hashes']['constructor'] and digest(HERE.parent/'第四十轮研究_20261028/sdp_support_pilot.py')==d['construction_code_hashes']['B40_model']
    assert digest(HERE.parent/'第四十轮研究_20261028/solver_environment.json')==d['solver_environment_sha256']
    for k,s in d['sources'].items():
        p=HERE.parent/s['path'];rr=p.read_bytes();assert hashlib.sha256(rr).hexdigest()==s['sha256'];data[k]=json.loads(rr)
        if 'check_sha256' in s:
            ck=(p.parent/'artifact_check.json').read_bytes();assert hashlib.sha256(ck).hexdigest()==s['check_sha256'] and json.loads(ck)['source_sha256']==s['sha256']
    b39=data['B39'];b40=data['B40'];prior={(z['source_parent_index'],z['source_child_index']):F(z['hybrid_local_loss_upper']) for z in b39['completed_subbands']};prior[b40['parent'],b40['child']]=F(b40['selected_local_hybrid_bound'])
    jobs=sorted(b39['completed_subbands'],key=lambda z:prior[z['source_parent_index'],z['source_child_index']],reverse=True);assert len(jobs)==d['total_subbands']==137
    bands={(p['source_beta_certificate_index'],j):c for p in data['B28']['beta_records'] for j,c in enumerate(p.get('subbands',[]))}
    cal=data['cal'];strips=data['excl']['strips'];cy=list(map(F,cal['calibration_centers']));ce=2*F(cal['calibration_raw_radius']);l0,l1,l2,l3=[v-ce for v in cy];u0,u1,u2,u3=[v+ce for v in cy]
    rlo=(l3-u2)/(u2-l1);rhi=(u3-l2)/(l2-u1);r0=(l3-u2)/(u2-u1);r1=(l3-l2)/(l2-u1);pieces=[([rlo,r0],[l3-u2,l3-u2,l0-u2],[F(0),F(1)]),([r0,r1],[l3-u1,l0-u1,l0-u1],[F(1),F(1)]),([r1,rhi],[l2-u1,l0-u1,F(0)],[F(1),F(0)])]
    # Pin inherited piece hashes before using their already checked bounds.
    sourcepieces={};b39dir=HERE.parent/'第三十九轮研究_20261027';oldverify=digest(b39dir/'verify_adaptive_direction_cover.py')
    for e in b39['completed_subbands']:
        pp=b39dir/e['path'];assert digest(pp)==e['sha256'];v=json.loads((b39dir/'分区核验'/(pp.stem+'.json')).read_text(encoding='utf-8'));assert v['source_sha256']==e['sha256'] and v['verification_code_sha256']==oldverify and v['hybrid_local_loss_upper']==e['hybrid_local_loss_upper'];sourcepieces[e['source_parent_index'],e['source_child_index']]=pp
    codehash=digest(Path(__file__));deps={n:digest(HERE.parent/p) for n,p in [('B28','第二十八轮研究_20261026/verify_support_polygon.py'),('B27','第二十七轮研究_20261026/verify_joint_restricted_performance.py'),('B38','第三十八轮研究_20261027/verify_independent_direction_pilot.py')]}
    cache=HERE/'分区核验';cache.mkdir(exist_ok=True);results=[];replacements={}
    for e,job in zip(d['completed_subbands'],jobs):
        key=(e['source_parent_index'],e['source_child_index']);assert key==(job['source_parent_index'],job['source_child_index']);piecepath=HERE/e['path'];h=digest(piecepath);assert h==e['sha256'];cached=cache/piecepath.name
        if cached.exists():
            z=json.loads(cached.read_text(encoding='utf-8'));assert z['source_sha256']==h and z['verification_code_sha256']==codehash and z['verification_dependencies']==deps
        else:
            obj=json.loads(piecepath.read_text(encoding='utf-8'));assert obj['construction_code_hashes']==d['construction_code_hashes'] and obj['source_piece_sha256']==job['sha256']
            z=verify_piece(obj,json.loads(sourcepieces[key].read_text(encoding='utf-8')),bands[key],prior[key],cal,pieces,strips,b40);z.update(source_sha256=h,verification_code_sha256=codehash,verification_dependencies=deps);save(cached,z)
        assert z['qualified_methods']==e['qualified_methods'] and z['hybrid_local_loss_upper']==e['hybrid_local_loss_upper'] and z['improved']==e['improved'];results.append(z);replacements[key]=F(z['hybrid_local_loss_upper'])
        print(json.dumps(dict(independently_checked=len(results),snapshot_completed=len(d['completed_subbands']),seconds=round(time.monotonic()-start,1))),flush=True)
    uppers=[]
    for parent in data['B28']['beta_records']:
        if not parent['allowed_pole_indices']:continue
        ix=parent['source_beta_certificate_index']
        value=max(replacements.get((ix,j),prior.get((ix,j),F(c['absolute_cost_upper']))) for j,c in enumerate(parent['subbands']) if c['allowed_pole_indices']) if parent.get('subbands') else F(parent['absolute_cost_upper'])
        uppers.append((value,ix))
    upper,attainer=max(uppers);assert upper==F(d['new_hybrid_global_upper'])<=F(d['prior_B40_global_upper']) and attainer==d['new_hybrid_global_attaining_parent']
    lowerpath=HERE.parent/'第三十四轮研究_20261027/standard_sqp_comparison.json';lr=lowerpath.read_bytes();lc=json.loads((lowerpath.parent/'artifact_check.json').read_text(encoding='utf-8'));assert hashlib.sha256(lr).hexdigest()==lc['source_sha256'];lowerdata=json.loads(lr);lower=F(lc['best_loss_lower']);assert upper>lower
    worldchecks=[]
    for world in lowerdata['certified_worlds']:
        beta=F(world['parameters']['beta']);pole=F(world['parameters']['report_pole']);loss=F(world['original_objective_gap_outer'][1]);matches=[]
        for parent in data['B28']['beta_records']:
            pl,pr=map(F,parent['beta_interval'])
            if not pl<=beta<=pr or not parent['allowed_pole_indices']:continue
            ix=parent['source_beta_certificate_index'];children=parent.get('subbands')
            for j,child in enumerate(children if children else [parent]):
                l,r=map(F,child['beta_interval'])
                allowed=[k for k in child['allowed_pole_indices'] if F(strips[k]['pole_interval'][0])<=pole<=F(strips[k]['pole_interval'][1])]
                if not l<=beta<=r or not allowed:continue
                bound=replacements.get((ix,j),prior.get((ix,j),F(child['absolute_cost_upper']))) if children else F(parent['absolute_cost_upper']);assert bound>=loss
                matches.append(dict(parent=ix,child=j if children else None,local_bound=str(bound),pole_strips=allowed))
        assert matches and upper>=loss;worldchecks.append(dict(source_index=world['index'],matching_local_certificates=matches))
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_interval_digits=110,completed_subbands_checked=len(results),total_subbands=137,
        qualified_method_candidates=sum(z['qualified_methods'] for z in results),solver_proposal_records=sum(z['proposal_records'] for z in results),exact_solver_solved_flags=sum(z['exact_solver_solved_flags'] for z in results),near_solved_flags=sum(z['near_solved_flags'] for z in results),construction_errors=sum(z['construction_errors'] for z in results),
        beta_subintervals=sum(z['beta_subintervals'] for z in results),persistence_cells=sum(z['persistence_cells'] for z in results),raw_calibration_cells=sum(z['raw_calibration_cells'] for z in results),improved_subbands=sum(z['improved'] for z in results),
        both_methods_qualified_subbands=sum(z['both_methods_qualified'] for z in results),SDP_smaller_certified_budget_subbands=sum(z['SDP_smaller_certified_budget'] for z in results),LP_smaller_certified_budget_subbands=sum(z['LP_smaller_certified_budget'] for z in results),
        inherited_B39_pieces_rehashed=137,original_parent_records_preserved=len(data['B28']['beta_records']),known_legal_world_local_checks=worldchecks,new_compatible_worlds_fitted=0,
        new_global_hybrid_upper=str(upper),known_compatible_loss_lower=str(lower),global_upper_to_known_lower_ratio=str(upper/lower),
        first_direction_inherited_not_autonomous=True,solver_failures_not_counted_as_physical_counterexamples=True,
        original_twelve_report_error_union_and_joint_calibration_preserved=True,full_51_variable_SDP_or_exact_optimum_claimed=False,equal_resource_superiority_claimed=False)
    if partial:
        snaps=HERE/'快照';snaps.mkdir(exist_ok=True);sp=snaps/f'cohort_{len(results):03d}.json'
        if sp.exists():assert sp.read_bytes()==raw
        else:sp.write_bytes(raw)
        out.update(partial_only=True,construction_snapshot_path=str(sp.relative_to(HERE)).replace('\\','/'),all_137_new_subband_method_comparisons_complete=False,global_cover_uses_inherited_safe_fallbacks=True)
        save(HERE/'partial_audit.json',out)
    else:
        assert d['construction_complete'] and len(results)==137 and d['reused_B40_method_pairs']==1 and d['newly_invoked_solver_runs']==272
        out.update(partial_only=False,all_137_new_subband_method_comparisons_complete=True,one_B40_method_pair_reused=True,new_method_attempts=272);save(HERE/'artifact_check.json',out)
    print(json.dumps(dict(partial=partial,checked=len(results),qualified=out['qualified_method_candidates'],upper=float(upper),ratio=float(upper/lower),seconds=round(time.monotonic()-start,1))))
if __name__=='__main__':main()
