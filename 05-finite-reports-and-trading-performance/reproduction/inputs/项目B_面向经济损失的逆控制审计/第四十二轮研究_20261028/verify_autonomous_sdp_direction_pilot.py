"""Independent SDP-only vertex selection and continuous certificate audit.
No constructor, proposal model, or solver imported.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第三十八轮研究_20261027'))
from verify_independent_direction_pilot import direction,arithmetic,Q
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');q.replace(p)
def piece(obj,old,band,baseline,cal,strips,pieces):
    assert obj['beta_interval']==band['beta_interval'] and obj['allowed_pole_indices']==band['allowed_pole_indices'] and obj['root_gap_upper']==band['root_gap_upper']
    axes=[z for z in band['facets'] if z['origin'] in ['nonnegative_perpendicular','inherited_B20_B27_axis']];assert axes==obj['axis_facets'] and len(axes)==4
    facets=axes.copy();used=set();lo,hi=map(F,band['beta_interval']);e0=arithmetic.root(arithmetic.b0,F(1));delta=F(band['root_gap_upper']);whole=Q(arithmetic.root(hi,F(1)).l,arithmetic.root(lo,F(1)).r);assert (whole-e0).r<=delta
    shifts=list(map(F,['0','1/10000000','1/1000000','1/100000','1/10000','3/10000','1/1000','3/1000','1/100','3/100','1/10']));pcells=bcells=rcells=qualified=newqualified=newruns=0;statuses={}
    assert 1<=len(obj['runs'])<=4
    for i,r in enumerate(obj['runs']):
        assert r['iteration']==i and r['facets_before']==facets and r['qualification_finished'] and r['reused_first_B41_SDP']==(i==0)
        vertices,_=arithmetic.polygon(facets);vertices.sort(key=lambda z:sum(v*v for v in z),reverse=True);expected=None
        for vertex in vertices:
            if vertex[1]<=0:continue
            a,b=direction(vertex)
            if (a,b) not in used:expected=vertex,a,b;break
        assert expected is not None;vertex,a,b=expected;assert list(map(F,r['selected_vertex']))==list(vertex) and list(map(F,r['direction']))==[a,b] and a*a+b*b==1 and b>=0;used.add((a,b))
        if 'proposal' in r:
            p=r['proposal'];assert p['method']=='SDP' and p['solver']=='Clarabel 0.11.1' and p['scalar_support_samples_included'];assert p['variables']==9+8*p['pole_samples']+978 and p['equalities']==258 and p['PSD_orders']==[19,16]*3;statuses[p['status']]=statuses.get(p['status'],0)+1
        else:assert 'construction_error_type' in r
        if i==0:
            s=next(z for z in old['methods'] if z['method']=='SDP');assert old['direction']==r['direction']
            for k in ['proposal','certified','qualification_trials','proof']:
                assert r.get(k)==s.get(k)
        else:newruns+=1
        if not r['certified']:assert 'proof' not in r and all(not z['passed'] for z in r['qualification_trials'])
        else:
            pr=r['proof'];M=F(pr['raw_budget_upper']);qualified+=1
            if i>0:
                tr=r['qualification_trials'];assert [F(z['shift']) for z in tr]==shifts[:len(tr)] and tr[-1]['passed'] and all(not z['passed'] for z in tr[:-1]);shift=F(pr['coefficient_shift']);assert shift==F(tr[-1]['shift']);cc=list(map(F,pr['coefficients']));assert cc==[F(z)+shift for z in r['proposal']['candidate_coefficients']] and min(cc[0],cc[4])>=0
                last=lo;assert len(pr['complete_beta_kernel_cover'])==16
                for j,row in enumerate(pr['complete_beta_kernel_cover']):
                    ll,hh=map(F,row['beta_interval']);assert ll==last==lo+(hi-lo)*F(j,16) and hh==lo+(hi-lo)*F(j+1,16);last=hh;es=[Q(arithmetic.root(hh,F(k)).l,arithmetic.root(ll,F(k)).r) for k in [1,2]];point=F(0);kp=row['kernel_proof']
                    for cell in kp['cover']:
                        l,h=map(F,cell['persistence']);assert l==point<h;point=h;polys,initial=arithmetic.polynomials(es,e0,cc,a,b,l,h,kp['method']);assert min(z.l for z in initial)>=0 and all(min(z.l for z in poly)>=0 for poly in polys);pcells+=1
                    assert point==1;bcells+=1
                assert last==hi;rcells+=arithmetic.raw_check(cc,M,pr['raw_calibration_proof'],cal,pieces,strips,band['allowed_pole_indices']);newqualified+=1
            facets.append(dict(origin='autonomous_SDP_vertex',iteration=i,halfspace=[str(a),str(b),str(M)]))
        _,radius=arithmetic.polygon(facets);assert radius==F(r['normalized_squared_radius_after'])
    _,radius=arithmetic.polygon(facets);loss=delta*delta*radius;hybrid=min(loss,baseline);assert obj['final_facets']==facets and radius==F(obj['final_squared_radius']) and loss==F(obj['new_local_loss_upper']) and baseline==F(obj['prior_B41_local_upper']) and hybrid==F(obj['hybrid_local_loss_upper'])
    return dict(SDP_direction_runs=len(obj['runs']),reused_first_SDP_runs=1,new_SDP_method_attempts=newruns,qualified_facets=qualified,newly_qualified_facets=newqualified,solver_status_counts=statuses,beta_subintervals=bcells,persistence_cells=pcells,raw_calibration_cells=rcells,hybrid_local_loss_upper=str(hybrid),improved=hybrid<baseline)
def main():
    start=time.monotonic();partial='--partial' in sys.argv;manifest=HERE/('autonomous_sdp_direction_checkpoint.json' if partial else 'autonomous_sdp_direction_pilot.json');raw=manifest.read_bytes();d=json.loads(raw)
    codes=d['construction_code_hashes'];assert sha(HERE/'autonomous_sdp_direction_pilot.py')==codes['constructor'] and sha(HERE.parent/'第四十一轮研究_20261028/standard_sdp_cohort.py')==codes['B41_helper'] and sha(HERE.parent/'第四十轮研究_20261028/sdp_support_pilot.py')==codes['B40_model']
    s=d['source_B41'];p=HERE.parent/s['path'];assert sha(p)==s['sha256'];source=json.loads(p.read_text(encoding='utf-8'));cp=p.parent/'artifact_check.json';assert sha(cp)==s['check_sha256'];oldcheck=json.loads(cp.read_text(encoding='utf-8'));assert oldcheck['source_sha256']==s['sha256'] and oldcheck['all_137_new_subband_method_comparisons_complete'] and not oldcheck['partial_only']
    data={}
    for k,s in d['source_inputs'].items():
        p=HERE.parent/s['path'];assert sha(p)==s['sha256'] and s==source['sources'][k];data[k]=json.loads(p.read_text(encoding='utf-8'))
    jobs=sorted(source['completed_subbands'],key=lambda z:F(z['hybrid_local_loss_upper']),reverse=True)[:12];assert d['total_pilot_subbands']==12;oldpieces={};vh=sha(HERE.parent/'第四十一轮研究_20261028/verify_standard_sdp_cohort.py')
    # Rehash the complete inherited cover before using any B41 fallback.
    for e in source['completed_subbands']:
        p=HERE.parent/'第四十一轮研究_20261028'/e['path'];assert sha(p)==e['sha256'];c=json.loads((p.parent.parent/'分区核验'/p.name).read_text(encoding='utf-8'));assert c['source_sha256']==e['sha256'] and c['verification_code_sha256']==vh and c['hybrid_local_loss_upper']==e['hybrid_local_loss_upper'];oldpieces[e['source_parent_index'],e['source_child_index']]=p
    bands={(p['source_beta_certificate_index'],j):z for p in data['B28']['beta_records'] for j,z in enumerate(p.get('subbands',[]))};cal=data['cal'];strips=data['excl']['strips'];cy=list(map(F,cal['calibration_centers']));ce=2*F(cal['calibration_raw_radius']);l0,l1,l2,l3=[z-ce for z in cy];u0,u1,u2,u3=[z+ce for z in cy];rlo=(l3-u2)/(u2-l1);rhi=(u3-l2)/(l2-u1);r0=(l3-u2)/(u2-u1);r1=(l3-l2)/(l2-u1);pieces=[([rlo,r0],[l3-u2,l3-u2,l0-u2],[F(0),F(1)]),([r0,r1],[l3-u1,l0-u1,l0-u1],[F(1),F(1)]),([r1,rhi],[l2-u1,l0-u1,F(0)],[F(1),F(0)])]
    codehash=sha(Path(__file__));deps={k:sha(HERE.parent/q) for k,q in [('B38','第三十八轮研究_20261027/verify_independent_direction_pilot.py'),('B28','第二十八轮研究_20261026/verify_support_polygon.py'),('B27','第二十七轮研究_20261026/verify_joint_restricted_performance.py')]};cache=HERE/'分区核验';cache.mkdir(exist_ok=True);results=[];replacements={(e['source_parent_index'],e['source_child_index']):F(e['hybrid_local_loss_upper']) for e in source['completed_subbands']}
    for e,job in zip(d['completed_subbands'],jobs):
        key=e['source_parent_index'],e['source_child_index'];assert key==(job['source_parent_index'],job['source_child_index']);pp=HERE/e['path'];h=sha(pp);assert h==e['sha256'];cached=cache/pp.name
        if cached.exists():z=json.loads(cached.read_text(encoding='utf-8'));assert z['source_sha256']==h and z['verification_code_sha256']==codehash and z['verification_dependencies']==deps
        else:
            obj=json.loads(pp.read_text(encoding='utf-8'));assert obj['construction_code_hashes']==codes and obj['source_B41_piece_sha256']==job['sha256'];z=piece(obj,json.loads(oldpieces[key].read_text(encoding='utf-8')),bands[key],F(job['hybrid_local_loss_upper']),cal,strips,pieces);z.update(source_sha256=h,verification_code_sha256=codehash,verification_dependencies=deps);save(cached,z)
        assert z['new_SDP_method_attempts']==e['new_SDP_method_attempts'] and z['qualified_facets']==e['qualified_facets'] and z['hybrid_local_loss_upper']==e['hybrid_local_loss_upper'] and z['improved']==e['improved'];replacements[key]=F(z['hybrid_local_loss_upper']);results.append(z);print(json.dumps(dict(independently_checked=len(results),snapshot_subbands=len(d['completed_subbands']),seconds=round(time.monotonic()-start,1))),flush=True)
    bounds=[]
    for p in data['B28']['beta_records']:
        if not p['allowed_pole_indices']:continue
        ix=p['source_beta_certificate_index'];val=max(replacements.get((ix,j),F(z['absolute_cost_upper'])) for j,z in enumerate(p['subbands']) if z['allowed_pole_indices']) if p.get('subbands') else F(p['absolute_cost_upper']);bounds.append((val,ix))
    upper,attainer=max(bounds);assert upper==F(d['new_global_hybrid_upper'])<=F(oldcheck['new_global_hybrid_upper']) and attainer==d['attaining_parent'];lowerdir=HERE.parent/'第三十四轮研究_20261027';lower=json.loads((lowerdir/'artifact_check.json').read_text(encoding='utf-8'));assert sha(lowerdir/'standard_sqp_comparison.json')==lower['source_sha256'];worlds=json.loads((lowerdir/'standard_sqp_comparison.json').read_text(encoding='utf-8'))['certified_worlds'];worldchecks=[]
    for w in worlds:
        beta=F(w['parameters']['beta']);pole=F(w['parameters']['report_pole']);loss=F(w['original_objective_gap_outer'][1]);matches=[]
        for p in data['B28']['beta_records']:
            ll,hh=map(F,p['beta_interval'])
            if not ll<=beta<=hh or not p['allowed_pole_indices']:continue
            ix=p['source_beta_certificate_index'];children=p.get('subbands')
            for j,z in enumerate(children if children else [p]):
                l,h=map(F,z['beta_interval']);allowed=[k for k in z['allowed_pole_indices'] if F(strips[k]['pole_interval'][0])<=pole<=F(strips[k]['pole_interval'][1])]
                if not l<=beta<=h or not allowed:continue
                bound=replacements.get((ix,j),F(z['absolute_cost_upper'])) if children else F(p['absolute_cost_upper']);assert loss<=bound;matches.append([ix,j if children else None])
        assert matches and loss<=upper;worldchecks.append(dict(index=w['index'],matching_local_certificates=matches))
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),partial_only=partial,independent_interval_digits=110,completed_subbands_checked=len(results),total_pilot_subbands=12,
        new_SDP_method_attempts=sum(z['new_SDP_method_attempts'] for z in results),reused_first_SDP_runs=len(results),newly_qualified_facets=sum(z['newly_qualified_facets'] for z in results),improved_subbands=sum(z['improved'] for z in results),beta_subintervals=sum(z['beta_subintervals'] for z in results),persistence_cells=sum(z['persistence_cells'] for z in results),raw_calibration_cells=sum(z['raw_calibration_cells'] for z in results),
        inherited_B41_pieces_rehashed=137,original_parent_records_preserved=len(data['B28']['beta_records']),known_legal_world_local_checks=worldchecks,new_global_hybrid_upper=str(upper),prior_global_upper=oldcheck['new_global_hybrid_upper'],global_upper_to_known_lower_ratio=str(upper/F(lower['best_loss_lower'])),
        SDP_only_accepted_facets_used_for_direction_selection=True,old_mixed_LP_facets_supplied_to_direction_selection=False,equal_resource_or_full_51_variable_optimum_claimed=False)
    if partial:
        snaps=HERE/'快照';snaps.mkdir(exist_ok=True);p=snaps/f'pilot_{len(results):02d}_{out["source_sha256"][:12]}.json'
        if p.exists():assert p.read_bytes()==raw
        else:p.write_bytes(raw)
        out['construction_snapshot_path']=str(p.relative_to(HERE)).replace('\\','/');save(HERE/'partial_audit.json',out)
    else:assert d['construction_complete'] and len(results)==12;save(HERE/'artifact_check.json',out)
    print(json.dumps(dict(partial=partial,checked=len(results),new_SDP=out['new_SDP_method_attempts'],upper=float(upper),ratio=float(upper/F(lower['best_loss_lower'])))))
if __name__=='__main__':main()
