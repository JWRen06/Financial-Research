"""B42 planned 12-band pilot: start at axis facets, choose new SDP angles.
Requires B41 all137 independently complete. First identical-axis SDP run
is reused explicitly; three more angles selected using only SDP facets.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time,os,traceback
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第四十一轮研究_20261028'))
import standard_sdp_cohort as cohort
model=cohort.model;base=cohort.base;support=cohort.support
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8');q.replace(p)
def pick(facets,used):
    for x,y in sorted(base.vertices(facets),key=lambda z:sum(v*v for v in z),reverse=True):
        if y<=0:continue
        norm=float(x*x+y*y)**.5;t=F(round(float(y)/(norm+float(x))*10**6),10**6);a=(1-t*t)/(1+t*t);b=2*t/(1+t*t)
        if (a,b) not in used:return x,y,a,b
    return None
def main():
    start=time.monotonic();support.n.p.iv.DEN=10**100;b41dir=HERE.parent/'第四十一轮研究_20261028';sourcepath=b41dir/'standard_sdp_cohort.json';checkpath=b41dir/'artifact_check.json';raw=sourcepath.read_bytes();source=json.loads(raw);ck=json.loads(checkpath.read_text(encoding='utf-8'))
    assert ck['source_sha256']==hashlib.sha256(raw).hexdigest() and ck['all_137_new_subband_method_comparisons_complete'] is True and not ck['partial_only']
    source_ref=dict(path=str(sourcepath.relative_to(HERE.parent)).replace('\\','/'),sha256=sha(sourcepath),check_sha256=sha(checkpath));data={}
    for k in ['B28','cal','excl']:
        s=source['sources'][k];p=HERE.parent/s['path'];assert sha(p)==s['sha256'];data[k]=json.loads(p.read_text(encoding='utf-8'))
    codes=dict(constructor=sha(Path(__file__)),B41_helper=sha(b41dir/'standard_sdp_cohort.py'),B40_model=sha(HERE.parent/'第四十轮研究_20261028/sdp_support_pilot.py'))
    assert codes['B41_helper']==source['construction_code_hashes']['constructor'] and codes['B40_model']==source['construction_code_hashes']['B40_model']
    jobs=sorted(source['completed_subbands'],key=lambda z:F(z['hybrid_local_loss_upper']),reverse=True)[:12];assert len(jobs)==12
    bands={(p['source_beta_certificate_index'],j):z for p in data['B28']['beta_records'] for j,z in enumerate(p.get('subbands',[]))};cal=data['cal'];strips=data['excl']['strips'];_,_,step,pieces=base.calmod.frontier(cal)
    manifest=HERE/'autonomous_sdp_direction_checkpoint.json';chunk=HERE/'分区证据';work=HERE/'进行中';chunk.mkdir(exist_ok=True);work.mkdir(exist_ok=True)
    if manifest.exists():
        old=json.loads(manifest.read_text(encoding='utf-8'));assert old['source_B41']==source_ref and old['construction_code_hashes']==codes;done=old['completed_subbands']
        for e in done:assert sha(HERE/e['path'])==e['sha256']
    else:done=[]
    def checkpoint(complete=False):
        replacement={(e['source_parent_index'],e['source_child_index']):F(e['hybrid_local_loss_upper']) for e in source['completed_subbands']}
        replacement.update({(e['source_parent_index'],e['source_child_index']):F(e['hybrid_local_loss_upper']) for e in done});values=[]
        for parent in data['B28']['beta_records']:
            if not parent['allowed_pole_indices']:continue
            ix=parent['source_beta_certificate_index'];val=max(replacement.get((ix,j),F(z['absolute_cost_upper'])) for j,z in enumerate(parent['subbands']) if z['allowed_pole_indices']) if parent.get('subbands') else F(parent['absolute_cost_upper']);values.append((val,ix))
        upper,attainer=max(values);assert upper<=F(ck['new_global_hybrid_upper'])
        d=dict(construction_complete=complete,source_B41=source_ref,construction_code_hashes=codes,source_inputs={k:source['sources'][k] for k in ['B28','cal','excl']},total_pilot_subbands=12,directions_per_subband_cap=4,first_identical_axis_SDP_run_reused=True,
            completed_subbands=done,new_global_hybrid_upper=str(upper),attaining_parent=attainer,prior_global_upper=ck['new_global_hybrid_upper'],
            scope='Pilot highest12 B41 bands. Start from four inherited axis facets only, reuse identical first SDP direction explicitly, select remaining directions using ONLY accepted SDP facets. No B39 LP mixed coefficients/facets supplied to direction selection. B41 whole cover retained as safe global fallback. Not equally resourced or full51-variable optimum.')
        save(manifest,d)
        if complete:save(HERE/'autonomous_sdp_direction_pilot.json',d)
    checkpoint()
    for ordinal,entry in enumerate(jobs):
        key=entry['source_parent_index'],entry['source_child_index']
        if ordinal<len(done):assert key==(done[ordinal]['source_parent_index'],done[ordinal]['source_child_index']);continue
        pp=b41dir/entry['path'];assert sha(pp)==entry['sha256'];old=json.loads(pp.read_text(encoding='utf-8'));band=bands[key];axes=[z for z in band['facets'] if z['origin'] in ['nonnegative_perpendicular','inherited_B20_B27_axis']];assert len(axes)==4
        job=work/f'{key[0]:03d}_{key[1]:02d}.json'
        if job.exists():obj=json.loads(job.read_text(encoding='utf-8'));assert obj['construction_code_hashes']==codes and obj['source_B41_piece_sha256']==entry['sha256']
        else:obj=dict(source_parent_index=key[0],source_child_index=key[1],source_B41_piece_sha256=entry['sha256'],construction_code_hashes=codes,beta_interval=band['beta_interval'],allowed_pole_indices=band['allowed_pole_indices'],root_gap_upper=band['root_gap_upper'],axis_facets=axes,runs=[])
        facets=axes.copy();used=set()
        for i in range(4):
            selected=pick(facets,used)
            if selected is None:break
            x,y,a,b=selected;used.add((a,b));run=obj['runs'][i] if i<len(obj['runs']) else None
            if run is not None:assert run['direction']==list(map(str,[a,b])) and run['facets_before']==facets
            else:
                run=dict(iteration=i,selected_vertex=list(map(str,[x,y])),direction=list(map(str,[a,b])),facets_before=facets.copy(),reused_first_B41_SDP=(i==0));obj['runs'].append(run)
            save(HERE/'runtime.json',dict(process_id=os.getpid(),completed_subbands=len(done),current_parent=key[0],current_child=key[1],iteration=i,phase='proposal_or_qualification'))
            if not run.get('qualification_finished'):
                if i==0:
                    assert old['direction']==run['direction'];reused=next(z for z in old['methods'] if z['method']=='SDP');run.update({k:reused[k] for k in ['proposal','certified','qualification_trials'] if k in reused})
                    if reused['certified']:run['proof']=reused['proof']
                else:
                    try:
                        if 'proposal' not in run:run['proposal']=model.program(band,cal,pieces,step,strips,[a,b],'SDP',True);save(job,obj)
                        if 'candidate_coefficients' in run['proposal']:run.update(cohort.qualify(run['proposal'],band,cal,pieces,step,strips,[a,b]))
                        else:run.update(certified=False,qualification_trials=[],no_finite_solver_candidate=True)
                    except BaseException as exc:
                        if isinstance(exc,(KeyboardInterrupt,SystemExit)):raise
                        run.update(certified=False,qualification_trials=[],construction_error_type=type(exc).__name__,construction_error=str(exc),traceback=traceback.format_exc())
                run['qualification_finished']=True;save(job,obj)
            if run['certified']:facets.append(dict(origin='autonomous_SDP_vertex',iteration=i,halfspace=run['direction']+[run['proof']['raw_budget_upper']]))
            run['normalized_squared_radius_after']=str(base.norm_upper(facets));save(job,obj)
        radius=base.norm_upper(facets);loss=F(band['root_gap_upper'])**2*radius;hybrid=min(F(entry['hybrid_local_loss_upper']),loss);obj.update(final_facets=facets,final_squared_radius=str(radius),new_local_loss_upper=str(loss),prior_B41_local_upper=entry['hybrid_local_loss_upper'],hybrid_local_loss_upper=str(hybrid));final=chunk/job.name;save(final,obj)
        done.append(dict(source_parent_index=key[0],source_child_index=key[1],path=str(final.relative_to(HERE)).replace('\\','/'),sha256=sha(final),reused_first_SDP_runs=1,new_SDP_method_attempts=max(0,len(obj['runs'])-1),qualified_facets=sum(r['certified'] for r in obj['runs']),hybrid_local_loss_upper=str(hybrid),improved=hybrid<F(entry['hybrid_local_loss_upper'])));checkpoint()
        print(json.dumps(dict(completed=len(done),total=12,parent=key[0],child=key[1],qualified=sum(r['certified'] for r in obj['runs']),improved=done[-1]['improved'],seconds=round(time.monotonic()-start,1))),flush=True)
    checkpoint(True);save(HERE/'runtime.json',dict(process_id=os.getpid(),phase='construction_complete',completed_subbands=len(done)));print('B42 pilot constructed; independent audit still required',flush=True)
if __name__=='__main__':main()
