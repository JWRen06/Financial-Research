"""Independent direction provenance, continuous kernels and raw budgets.
Uses only B28 independent verification arithmetic, never synthesis code.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第二十八轮研究_20261026'))
import verify_support_polygon as arithmetic
Q=arithmetic.Q
def direction(vertex):
    x,y=vertex
    with localcontext() as ctx:
        ctx.prec=180
        def dd(z):return Decimal(z.numerator)/Decimal(z.denominator)
        norm=dd(x*x+y*y).sqrt();t=F(int((dd(y)/(norm+dd(x))*10**6).to_integral_value(rounding=ROUND_HALF_EVEN)),10**6)
    return (1-t*t)/(1+t*t),2*t/(1+t*t)
def main():
    start=time.monotonic();raw=(HERE/'independent_direction_pilot.json').read_bytes();d=json.loads(raw);data={}
    for key,ss in d['sources'].items():
        p=HERE.parent/ss['path'];rr=p.read_bytes();assert hashlib.sha256(rr).hexdigest()==ss['sha256'];data[key]=json.loads(rr)
        if 'check_sha256' in ss:
            ck=(p.parent/'artifact_check.json').read_bytes();assert hashlib.sha256(ck).hexdigest()==ss['check_sha256'] and json.loads(ck)['source_sha256']==ss['sha256']
    assert d['source_parent_index']==data['B36']['updated_global_attaining_parent'];parent=next(z for z in data['B28']['beta_records'] if z['source_beta_certificate_index']==d['source_parent_index'])
    band=parent['subbands'][d['source_child_index']];assert d['beta_interval']==band['beta_interval'] and d['allowed_pole_indices']==band['allowed_pole_indices']
    axes=[z for z in band['facets'] if z['origin'] in ['nonnegative_perpendicular','inherited_B20_B27_axis']];assert d['initial_axis_facets']==axes and len(axes)==4
    facets=axes.copy();_,radius=arithmetic.polygon(facets);assert radius==F(d['axis_normalized_squared_radius']);lo,hi=map(F,d['beta_interval']);e0=arithmetic.root(arithmetic.b0,F(1))
    cal=data['calibration'];strips=data['exclusion']['strips'];cy=list(map(F,cal['calibration_centers']));ce=2*F(cal['calibration_raw_radius']);l0,l1,l2,l3=[v-ce for v in cy];u0,u1,u2,u3=[v+ce for v in cy]
    rlo=(l3-u2)/(u2-l1);rhi=(u3-l2)/(l2-u1);r0=(l3-u2)/(u2-u1);r1=(l3-l2)/(l2-u1)
    pieces=[([rlo,r0],[l3-u2,l3-u2,l0-u2],[F(0),F(1)]),([r0,r1],[l3-u1,l0-u1,l0-u1],[F(1),F(1)]),([r1,rhi],[l2-u1,l0-u1,F(0)],[F(1),F(0)])]
    used=set();pcells=bcells=rawcells=nproved=success=0
    for index,run in enumerate(d['runs']):
        assert run['iteration']==index and run['current_polygon_facets_before']==facets
        vertices,_=arithmetic.polygon(facets);vertices.sort(key=lambda z:sum(v*v for v in z),reverse=True);expected=None
        for vertex in vertices:
            if vertex[1]<=0:continue
            aa,bb=direction(vertex)
            if (aa,bb) not in used:expected=vertex,aa,bb;break
        assert expected is not None;vertex,a,b=expected;assert list(map(F,run['selected_vertex']))==list(vertex) and list(map(F,run['direction']))==[a,b] and a*a+b*b==1 and b>=0;used.add((a,b))
        if run['solver_success']:assert run['solver_status']==0;success+=1
        if 'proof' in run:
            z=run['proof'];cc=list(map(F,z['coefficients']));shift=F(z['coefficient_shift']);seed=list(map(F,run['unshifted_coefficients']));assert cc==[v+shift for v in seed] and min(cc[0],cc[4])>=0
            assert len(z['complete_beta_kernel_cover'])==z['beta_divisions'];last=lo
            for j,row in enumerate(z['complete_beta_kernel_cover']):
                ll,hh=map(F,row['beta_interval']);assert ll==last==lo+(hi-lo)*F(j,z['beta_divisions']) and hh==lo+(hi-lo)*F(j+1,z['beta_divisions']);last=hh
                es=[Q(arithmetic.root(hh,F(k)).l,arithmetic.root(ll,F(k)).r) for k in [1,2]];point=F(0);proof=row['kernel_proof']
                for cell in proof['cover']:
                    left,right=map(F,cell['persistence']);assert left==point<right;point=right
                    tests,initial=arithmetic.polynomials(es,e0,cc,a,b,left,right,proof['method']);assert min(v.l for v in initial)>=0 and all(min(v.l for v in poly)>=0 for poly in tests);pcells+=1
                assert point==1;bcells+=1
            assert last==hi
            rawcells+=arithmetic.raw_check(cc,F(z['raw_budget_upper']),z['raw_calibration_proof'],cal,pieces,strips,d['allowed_pole_indices'])
            facets.append(dict(origin='independent_vertex_direction',id=index,halfspace=[str(a),str(b),z['raw_budget_upper']]));nproved+=1
        _,radius=arithmetic.polygon(facets);assert radius==F(run['safe_normalized_squared_radius_after'])
    vertices,radius=arithmetic.polygon(facets);assert facets==d['final_facets'] and [list(map(str,z)) for z in vertices]==d['final_normalized_vertices'] and radius==F(d['final_normalized_squared_radius'])
    delta=F(d['root_gap_upper']);es=[Q(arithmetic.root(hi,F(k)).l,arithmetic.root(lo,F(k)).r) for k in [1,2]];assert (es[0]-e0).r<=delta and d['root_gap_upper']==band['root_gap_upper']
    assert delta*delta*F(d['axis_normalized_squared_radius'])==F(d['local_axis_only_loss_upper']) and delta*delta*radius==F(d['local_direction_discovery_loss_upper'])
    source=next(z for z in data['B35']['refined_subbands'] if z['source_parent_index']==d['source_parent_index'] and z['source_child_index']==d['source_child_index']);assert source['new_absolute_cost_upper']==d['inherited_B35_local_upper']
    out=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_interval_digits=110,
        standard_LP_runs=len(d['runs']),standard_LP_successes=success,new_vertex_selected_directions=len(used),
        continuously_proved_new_facets=nproved,complete_beta_subintervals_checked=bcells,
        persistence_cells_checked=pcells,original_raw_calibration_quadratic_cells_checked=rawcells,
        original_axis_seed_facets_only=4,old_mixed_directions_or_coefficients_used_for_selection=False,
        beta_partition_and_pole_qualification_inherited=True,full_original_global_bound_updated=False,
        independent_SDP_or_full_problem_solution_claimed=False,algorithm_superiority_claimed=False,
        local_axis_upper=d['local_axis_only_loss_upper'],local_discovered_direction_upper=d['local_direction_discovery_loss_upper'],
        inherited_B35_local_upper=d['inherited_B35_local_upper'])
    (HERE/'artifact_check.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('_upper')}));print(json.dumps(dict(local_upper=float(delta*delta*radius),relative_to_old=float(delta*delta*radius/F(d['inherited_B35_local_upper'])),seconds=round(time.monotonic()-start,1))))
if __name__=='__main__':main()
