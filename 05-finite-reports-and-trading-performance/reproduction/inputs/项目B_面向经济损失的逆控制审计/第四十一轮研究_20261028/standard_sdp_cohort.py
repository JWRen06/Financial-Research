"""B41: controlled LP/SDP coefficient synthesis for all 137 B39 first angles.
One already-completed B40 pair is reused and counted explicitly. Every other
pair is newly solved, continuously qualified, and checkpointed per method.
No new result is integrated before independent certificate reconstruction.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time,os,traceback
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第四十轮研究_20261028'))
import sdp_support_pilot as model
base=model.base;support=model.support

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(path,d):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(path)

def qualify(run,band,cal,pieces,step,strips,direction):
    a,b=direction;lo,hi=map(F,band['beta_interval']);e0=support.n.p.gb.root(support.reference,support.rho,F(1));boxes=[(lo+(hi-lo)*F(j,16),lo+(hi-lo)*F(j+1,16)) for j in range(16)]
    es=[support.n.roots(x) for x in boxes];seed=list(map(F,run['candidate_coefficients']));trials=[]
    shifts=[F(0),F(1,10**7),F(1,10**6),F(1,10**5),F(1,10000),F(3,10000),F(1,1000),F(3,1000),F(1,100),F(3,100),F(1,10)]
    for shift in shifts:
        cc=[x+shift for x in seed];cover=[]
        for (ll,hh),roots in zip(boxes,es):
            kp=support.cover(roots,e0,cc,a,b,max_depth=16,max_cells=2500)
            if kp is None:break
            cover.append(dict(beta_interval=[str(ll),str(hh)],kernel_proof=kp))
        passed=len(cover)==16;trials.append(dict(shift=str(shift),passed=passed,first_unproved_beta_index=None if passed else len(cover)))
        if passed:
            M,raw=base.budget(cc,cal,pieces,step,strips,band['allowed_pole_indices'])
            return dict(certified=True,qualification_trials=trials,proof=dict(coefficients=list(map(str,cc)),coefficient_shift=str(shift),complete_beta_kernel_cover=cover,raw_budget_upper=str(M),raw_calibration_proof=raw))
    return dict(certified=False,qualification_trials=trials)

def global_cover(data,manifest,replacements):
    old={(z['source_parent_index'],z['source_child_index']):F(z['hybrid_local_loss_upper']) for z in manifest['completed_subbands']}
    b40=data['B40'];old[b40['parent'],b40['child']]=F(b40['selected_local_hybrid_bound']);bounds=[]
    for parent in data['B28']['beta_records']:
        if not parent['allowed_pole_indices']:continue
        ix=parent['source_beta_certificate_index']
        if parent.get('subbands'):
            value=max(replacements.get((ix,j),old.get((ix,j),F(c['absolute_cost_upper']))) for j,c in enumerate(parent['subbands']) if c['allowed_pole_indices'])
        else:value=F(parent['absolute_cost_upper'])
        bounds.append((value,ix))
    return max(bounds)

def main():
    start=time.monotonic();support.n.p.iv.DEN=10**100;data={};sources={}
    for key,name in [('B39','第三十九轮研究_20261027/adaptive_direction_cover.json'),('B40','第四十轮研究_20261028/qualified_sdp_support.json'),('B28','第二十八轮研究_20261026/polygon_performance_certificate.json'),('cal','第七轮研究_20261023/physical_step_absolute_audit.json'),('excl','第二十七轮研究_20261026/joint_compatibility_exclusion.json')]:
        p=HERE.parent/name;rr=p.read_bytes();data[key]=json.loads(rr);sources[key]=dict(path=name,sha256=hashlib.sha256(rr).hexdigest());ck=p.parent/'artifact_check.json'
        if key in ['B39','B40','B28']:
            cc=ck.read_bytes();assert json.loads(cc)['source_sha256']==sources[key]['sha256'];sources[key]['check_sha256']=hashlib.sha256(cc).hexdigest()
    codes=dict(constructor=digest(Path(__file__)),B40_model=digest(HERE.parent/'第四十轮研究_20261028/sdp_support_pilot.py'))
    environment=HERE.parent/'第四十轮研究_20261028/solver_environment.json';envhash=digest(environment)
    assert model.clarabel.__version__=='0.11.1'
    source=data['B39'];b40=data['B40'];b40pilot=json.loads((HERE.parent/'第四十轮研究_20261028/sdp_support_pilot.json').read_text(encoding='utf-8'))
    assert digest(HERE.parent/'第四十轮研究_20261028/sdp_support_pilot.json')==b40['pilot_sha256'] and codes['B40_model']==b40['pilot_construction_code_sha256']
    cal=data['cal'];strips=data['excl']['strips'];_,_,step,pieces=base.calmod.frontier(cal)
    bands={(p['source_beta_certificate_index'],j):c for p in data['B28']['beta_records'] for j,c in enumerate(p.get('subbands',[]))}
    baseline={(z['source_parent_index'],z['source_child_index']):F(z['hybrid_local_loss_upper']) for z in source['completed_subbands']};baseline[b40['parent'],b40['child']]=F(b40['selected_local_hybrid_bound'])
    jobs=sorted(source['completed_subbands'],key=lambda z:baseline[z['source_parent_index'],z['source_child_index']],reverse=True);assert len(jobs)==137
    outpath=HERE/'standard_sdp_cohort_checkpoint.json';chunk=HERE/'分区证据';work=HERE/'进行中';chunk.mkdir(exist_ok=True);work.mkdir(exist_ok=True)
    if outpath.exists():
        d=json.loads(outpath.read_text(encoding='utf-8'));assert d['sources']==sources and d['construction_code_hashes']==codes and d['solver_environment_sha256']==envhash
        done=d['completed_subbands']
        for z in done:assert digest(HERE/z['path'])==z['sha256']
    else:done=[]
    lowercheck=json.loads((HERE.parent/'第三十四轮研究_20261027/artifact_check.json').read_text(encoding='utf-8'));knownlower=F(lowercheck['best_loss_lower'])
    def checkpoint(complete=False):
        replacements={(z['source_parent_index'],z['source_child_index']):F(z['hybrid_local_loss_upper']) for z in done};upper,attainer=global_cover(data,source,replacements)
        assert knownlower<=upper<=F(b40['new_global_hybrid_upper'])
        d=dict(construction_complete=complete,sources=sources,construction_code_hashes=codes,solver_environment_sha256=envhash,total_subbands=137,
            first_directions_inherited_from_B39=True,completed_subbands=done,new_hybrid_global_upper=str(upper),new_hybrid_global_attaining_parent=attainer,prior_B40_global_upper=b40['new_global_hybrid_upper'],
            newly_invoked_solver_runs=sum(z['new_solver_runs'] for z in done),reused_B40_method_pairs=sum(z['reused_B40_pair'] for z in done),
            continuously_constructed_facets=sum(z['qualified_methods'] for z in done),improved_subbands=sum(z['improved'] for z in done),
            scope='All 137 inherited B39 first directions, same original reports/error/calibration, shared B40 LP/SDP coefficient basis and sampling safeguards. One B40 pair reused explicitly. Other methods newly solved and continuously qualified; failures preserve old bounds. Construction only until independent whole-cover reconstruction. Not autonomous SDP directions, full 51-variable solve, sharp optimum or equal-resource superiority.')
        save(outpath,d)
        if complete:save(HERE/'standard_sdp_cohort.json',d)
    checkpoint()
    for ordinal,entry in enumerate(jobs):
        key=(entry['source_parent_index'],entry['source_child_index']);parent,child=key
        if ordinal<len(done):assert (done[ordinal]['source_parent_index'],done[ordinal]['source_child_index'])==key;continue
        sourcepath=HERE.parent/'第三十九轮研究_20261027'/entry['path'];assert digest(sourcepath)==entry['sha256'];part=json.loads(sourcepath.read_text(encoding='utf-8'));band=bands[key];direction=list(map(F,part['runs'][0]['direction']))
        runtime=dict(process_id=os.getpid(),current_job_index=ordinal,source_parent_index=parent,source_child_index=child,phase='starting',completed_subbands=len(done),wall_seconds=time.monotonic()-start)
        save(HERE/'runtime.json',runtime);jobfile=work/f'{parent}_{child:02d}.json'
        if jobfile.exists():
            obj=json.loads(jobfile.read_text(encoding='utf-8'));assert obj['construction_code_hashes']==codes and obj['source_piece_sha256']==entry['sha256']
        else:
            obj=dict(source_parent_index=parent,source_child_index=child,source_piece_sha256=entry['sha256'],construction_code_hashes=codes,beta_interval=band['beta_interval'],allowed_pole_indices=band['allowed_pole_indices'],root_gap_upper=band['root_gap_upper'],direction=list(map(str,direction)),prior_local_loss_upper=str(baseline[key]),methods=[])
        reused=key==(b40['parent'],b40['child'])
        if reused and not obj['methods']:
            for z in b40['methods']:
                row=dict(method=z['method'],proposal=b40pilot['runs'][z['pilot_run_index']],reused_from_B40=True,certified=z['certified'],qualification_trials=z['qualification_trials'])
                if z['certified']:row['proof']=z['proof']
                obj['methods'].append(row)
            save(jobfile,obj)
        if not reused:
            for method in ['LP','SDP']:
                z=next((q for q in obj['methods'] if q['method']==method),None)
                if z is not None and z.get('qualification_finished'):continue
                if z is None:z=dict(method=method,reused_from_B40=False);obj['methods'].append(z)
                runtime['phase']='proposal_'+method;save(HERE/'runtime.json',runtime)
                try:
                    if 'proposal' not in z:z['proposal']=model.program(band,cal,pieces,step,strips,direction,method,True);save(jobfile,obj)
                    runtime['phase']='continuous_qualification_'+method;save(HERE/'runtime.json',runtime)
                    if 'candidate_coefficients' in z['proposal']:z.update(qualify(z['proposal'],band,cal,pieces,step,strips,direction))
                    else:z.update(certified=False,qualification_trials=[],no_finite_solver_candidate=True)
                except BaseException as exc:
                    if isinstance(exc,(KeyboardInterrupt,SystemExit)):raise
                    z.update(certified=False,qualification_trials=z.get('qualification_trials',[]),construction_error_type=type(exc).__name__,construction_error=str(exc),traceback=traceback.format_exc())
                z['qualification_finished']=True;save(jobfile,obj)
        assert len(obj['methods'])==2
        facets=list(part['final_facets']);axes=[x for x in band['facets'] if x['origin'] in ['nonnegative_perpendicular','inherited_B20_B27_axis']];delta=F(band['root_gap_upper']);qualified=0
        for z in obj['methods']:
            if not z['certified']:continue
            M=F(z['proof']['raw_budget_upper']);facet=dict(origin='B41_'+z['method'],halfspace=[*obj['direction'],str(M)]);facets.append(facet);qualified+=1
            radius=base.norm_upper(axes+[facet]);z['axis_plus_one_facet_loss_upper']=str(delta*delta*radius)
            oldproof=part['runs'][0].get('proof');z['old_B39_same_direction_budget']=None if oldproof is None else oldproof['raw_budget_upper'];z['improved_old_same_direction_budget']=bool(oldproof and M<F(oldproof['raw_budget_upper']))
        radius=base.norm_upper(facets);new=delta*delta*radius;hybrid=min(baseline[key],new)
        obj.update(final_facets=facets,final_squared_radius=str(radius),new_local_loss_upper=str(new),hybrid_local_loss_upper=str(hybrid),qualified_methods=qualified,reused_B40_pair=reused)
        finalpath=chunk/jobfile.name;save(finalpath,obj)
        done.append(dict(source_parent_index=parent,source_child_index=child,path=str(finalpath.relative_to(HERE)).replace('\\','/'),sha256=digest(finalpath),new_solver_runs=0 if reused else 2,reused_B40_pair=reused,qualified_methods=qualified,hybrid_local_loss_upper=str(hybrid),prior_local_loss_upper=str(baseline[key]),improved=hybrid<baseline[key]))
        checkpoint();print(json.dumps(dict(completed=len(done),total=137,parent=parent,child=child,statuses=[z['proposal']['status'] if 'proposal' in z else z.get('construction_error_type') for z in obj['methods']],qualified=qualified,improved=hybrid<baseline[key],seconds=round(time.monotonic()-start,1))),flush=True)
    checkpoint(True);save(HERE/'runtime.json',dict(process_id=os.getpid(),phase='construction_complete',completed_subbands=len(done),wall_seconds=time.monotonic()-start));print('B41 all 137 construction subbands complete; independent audit still required',flush=True)
if __name__=='__main__':main()
