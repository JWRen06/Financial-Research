"""Independent 110-digit reconstruction of new support halfspaces.
Does not import B28 proposal, constructor, cover or vertex routines.
Previously checked B20/B27/B26 sources are reused only under pinned hashes.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
from itertools import combinations
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'第二十七轮研究_20261026'))
import verify_joint_restricted_performance as arithmetic
Q=arithmetic.Q;power=arithmetic.power;root=arithmetic.root
b0=F(9999,10000);rho=F(1,15)
def elevate(p,n):
    m=len(p)-1
    return [sum((p[j]*F(comb(m,j)*comb(n-m,i-j),comb(n,i)) for j in range(max(0,i-n+m),min(m,i)+1)),Q(0)) for i in range(n+1)]
def add(p,q):
    n=max(len(p),len(q))-1
    return [a+b for a,b in zip(elevate(p,n),elevate(q,n))]
def scale(p,z):return [x*z for x in p]
def multiply(p,q):
    m,n=len(p)-1,len(q)-1
    return [sum((p[j]*q[i-j]*F(comb(m,j)*comb(n,i-j),comb(m+n,i)) for j in range(max(0,i-n),min(m,i)+1)),Q(0)) for i in range(m+n+1)]
def product(*args):
    p=[Q(1)]
    for q in args:p=multiply(p,q)
    return p
def sqrtbox(v):
    assert v.l>0;den=10**105
    return Q(F(isqrt(v.l.numerator*den*den//v.l.denominator),den),F(isqrt(v.r.numerator*den*den//v.r.denominator)+1,den))
def polynomials(es,e0,c,a,b,left,right,method):
    # Reconstruct shapes directly as sums of affine powers of persistence;
    # the constructor instead translates monomial coefficients first.
    s=[Q(left),Q(right)];u=[Q(1-left),Q(1-right)]
    ds=[add([Q(rho)],scale(u,k*(1-e))) for k,e in zip([1,2],es)]
    d0=add([Q(rho)],scale(u,1-e0));ee=add([Q(1)],scale(s,-b0*es[0]));p=[Q(0)];initial=[]
    for g,e in enumerate(es):
        shape=[Q(0)]
        for n in range(4):
            for j in range(n+1):
                sp=[Q(1)]
                for _ in range(j):sp=multiply(sp,s)
                shape=add(shape,scale(sp,c[4*g+n]*power(e,n-j)))
        p=add(p,scale(multiply(shape,ds[1-g]),rho*(1-e)))
        initial.append(sum((c[4*g+n]*power(e,n+1) for n in range(4)),Q(0)))
    e=es[0];cc=(1-e)*(1-e0)/rho;rr=1-e+cc;v=1-b0*e*e;t=v+b0*rr*rr
    bn=v*(e+rr)+rr*cc;nn=add([t],scale(s,-b0*bn))
    x=add(scale(product(d0,ee,p),sqrtbox(2*e0*v*t)),scale(multiply(ds[1],nn),-a*rho*rho))
    if b==0:return [x],initial
    w=e*(2-e*(1+b0))+rr*rr
    uu=add(scale(add([Q(1-b0)],scale(u,2*b0)),(1-b0*e)*(1-b0*e)),scale(multiply(u,u),b0*b0*w))
    survival=[Q(1-b0*left*left),Q(1-b0*left*right),Q(1-b0*right*right)]
    if method=='inherited_axis_factorisation':
        assert (a,b)==(0,1)
        sq=add(scale(product(d0,d0,ee,ee,survival,p,p),2*e0*t),scale(product(s,s,uu,ds[1],ds[1]),-rho**4*(1+cc)*(1+cc)*b0))
        return [p,sq],initial
    sq=add(product(survival,x,x),scale(product(s,s,uu,ds[1],ds[1]),-b*b*rho**4*(1+cc)*(1+cc)*b0*v))
    return [x,sq],initial
def raw_check(c,m,rows,cal,pieces,strips,allowed):
    step=F(cal['known_physical_step']);ys=[F(y) for row in cal['main_centers'] for y in row];eps=2*F(cal['main_raw_radius'])
    assert m>=0 and len(c)==len(ys)==8
    assert set(r['pole_index'] for r in rows)==set(allowed)
    count=0
    for j in allowed:
        previous=F(strips[j]['pole_interval'][0]);end=F(strips[j]['pole_interval'][1])
        for row in [r for r in rows if r['pole_index']==j]:
            lo,hi=map(F,row['pole_interval']);assert lo==previous<hi;previous=hi
            interval,p,den=pieces[row['gain_piece']];assert interval[0]<=lo<hi<=interval[1]
            signs=row['raw_signs'];assert len(signs)==8;aa=bb=F(0)
            for n,(ci,sg,y) in enumerate(zip(c,signs,ys)):
                slope=-c[n+1] if n%4<3 else F(0)
                assert sg in [-1,1] and sg*(ci+slope*lo)>=0 and sg*(ci+slope*hi)>=0
                aa+=ci*(y+eps*sg);bb+=slope*(y+eps*sg)
            assert [aa,bb]==list(map(F,row['raw_numerator']))
            d0,d1=den;prod=[d0*aa,d0*bb+d1*aa,d1*bb];poly=[m*x-step*y for x,y in zip(p,prod)]
            minimum=arithmetic.quadratic(poly,lo,hi)
            assert poly==list(map(F,row['quadratic'])) and minimum==F(row['minimum'])>=0
            at=F(row['minimizer']);assert lo<=at<=hi and sum(poly[k]*at**k for k in range(3))==minimum
            assert arithmetic.quadratic(p,lo,hi)>0 and min(d0+d1*lo,d0+d1*hi)>0;count+=1
        assert previous==end
    return count
def polygon(facets):
    # Independently enumerate all feasible boundary intersections, including
    # degenerate duplicates. Old axes give compactness and the origin feasibility.
    hs=[tuple(map(F,z['halfspace'])) for z in facets];vv=set()
    for (a,b,m),(c,d,n) in combinations(hs,2):
        det=a*d-b*c
        if det==0:continue
        pt=((m*d-b*n)/det,(a*n-m*c)/det)
        if all(h*pt[0]+k*pt[1]<=v for h,k,v in hs):vv.add(pt)
    assert vv and all(v>=0 for _,_,v in hs)
    return sorted(vv),max(x*x+y*y for x,y in vv)
def main():
    started=time.monotonic();raw=(HERE/'polygon_performance_certificate.json').read_bytes();out=json.loads(raw);assert out['complete'] is True
    paths=dict(B20='第二十轮研究_20261025/global_directional_certificate.json',B20_check='第二十轮研究_20261025/artifact_check.json',B27='第二十七轮研究_20261026/joint_restricted_performance.json',B27_check='第二十七轮研究_20261026/artifact_check.json',exclusion='第二十七轮研究_20261026/joint_compatibility_exclusion.json',calibration='第七轮研究_20261023/physical_step_absolute_audit.json',lower='第二十六轮研究_20261026/continuous_five_node_witness.json',lower_check='第二十六轮研究_20261026/artifact_check.json')
    data={}
    for key,path in paths.items():
        rr=(HERE.parent/path).read_bytes();assert hashlib.sha256(rr).hexdigest()==out['source_hashes'][key];data[key]=json.loads(rr)
    for key in ['B20','B27','lower']:assert out['source_hashes'][key]==data[key+'_check']['source_sha256']
    assert data['B27_check']['exclusion_sha256']==out['source_hashes']['exclusion']
    for file,expected in out['construction_code_hashes'].items():assert hashlib.sha256((HERE/file).read_bytes()).hexdigest()==expected
    old=data['B20'];prior=data['B27'];strips=data['exclusion']['strips'];cal=data['calibration'];pieces=[]
    # Rebuild calibration frontier directly from the original twelve report data.
    cy=list(map(F,cal['calibration_centers']));ce=2*F(cal['calibration_raw_radius']);l0,l1,l2,l3=[v-ce for v in cy];u0,u1,u2,u3=[v+ce for v in cy]
    rlo=(l3-u2)/(u2-l1);rhi=(u3-l2)/(l2-u1);r0=(l3-u2)/(u2-u1);r1=(l3-l2)/(l2-u1)
    pieces=[([rlo,r0],[l3-u2,l3-u2,l0-u2],[F(0),F(1)]),([r0,r1],[l3-u1,l0-u1,l0-u1],[F(1),F(1)]),([r1,rhi],[l2-u1,l0-u1,F(0)],[F(1),F(0)])]
    e0=root(b0,F(1));totalcells=rawcells=facetcount=children_count=polygons=0;uppers=[];previous=F(old['beta_outer'][0])
    assert len(out['beta_records'])==len(prior['beta_records'])==430
    for idx,(row,parent) in enumerate(zip(out['beta_records'],prior['beta_records'])):
        assert row['source_beta_certificate_index']==parent['source_beta_certificate_index']==idx
        assert row['beta_interval']==parent['beta_interval']==old['certificates'][idx]['beta_interval']
        lo,hi=map(F,row['beta_interval']);assert lo==previous<hi;previous=hi
        allowed=[j for j,s in enumerate(strips) if hi>F(s['excluded_beta_interval'][1])];assert allowed==row['allowed_pole_indices']==parent['allowed_pole_indices']
        if row['status']=='inherited_joint_exclusion':assert not allowed;continue
        if row['status']=='inherited_performance':assert F(row['absolute_cost_upper'])==F(parent['absolute_cost_upper']);uppers.append(F(row['absolute_cost_upper']));continue
        assert row['status']=='adaptive_support_cover';point=lo;childbounds=[]
        for child in row['subbands']:
            cl,ch=map(F,child['beta_interval']);assert cl==point<ch<=hi;point=ch;children_count+=1
            ca=[j for j,s in enumerate(strips) if ch>F(s['excluded_beta_interval'][1])];assert ca==child['allowed_pole_indices'] and set(ca)<=set(allowed)
            if child['status']=='inherited_joint_exclusion':assert not ca;continue
            assert ca
            es=[Q(root(ch,F(k)).l,root(cl,F(k)).r) for k in [1,2]];delta=F(child['root_gap_upper']);gap=es[0]-e0
            assert 0<gap.l<=gap.r<=delta<=F(old['certificates'][idx]['root_gap_outer'][1])
            axes=[dict(origin='nonnegative_perpendicular',halfspace=['0','-1','0'])]
            for kind,a,b in [('plus',1,0),('minus',-1,0),('perpendicular',0,1)]:
                key=parent['components'][kind];mm=F(prior['registry'][key]['raw_upper']);assert mm>=0
                axes.append(dict(origin='inherited_B20_B27_axis',kind=kind,source_budget_key=key,halfspace=list(map(str,[a,b,mm]))))
            if child['status']=='inherited_axes_subband':vv,radius=polygon(axes)
            else:
                assert child['status']=='new_support_polygon' and child['facets'][:4]==axes;polygons+=1
                assert len(child['facets'])==4+len(child['new_facets'])
                for j,facet in enumerate(child['new_facets']):
                    assert facet['id']==j;c=list(map(F,facet['normalized_coefficients']));assert len(c)==8 and min(c[0],c[4])>=0
                    a,b=map(F,facet['direction']);assert a*a+b*b==1 and b>=0;proof=facet['kernel_proof'];last=F(0)
                    if proof['method']=='inherited_axis_factorisation':assert (a,b) in [(1,0),(-1,0),(0,1)]
                    else:assert proof['method']=='mixed_support_factorisation'
                    for cell in proof['cover']:
                        left,right=map(F,cell['persistence']);assert left==last<right;last=right
                        tests,initial=polynomials(es,e0,c,a,b,left,right,proof['method'])
                        assert all(v.l>=0 for v in initial) and all(min(v.l for v in p)>=0 for p in tests),(idx,cl,ch,j,left,right)
                        totalcells+=1
                    assert last==1
                    mm=F(facet['raw_budget_upper']);rawcells+=raw_check(c,mm,facet['raw_calibration_proof'],cal,pieces,strips,ca)
                    assert child['facets'][4+j]==dict(origin='new_support_facet',id=j,halfspace=list(map(str,[a,b,mm])))
                    facetcount+=1
                vv,radius=polygon(child['facets']);assert vv==[tuple(map(F,z)) for z in child['normalized_vertices']]
                assert radius==F(child['normalized_squared_radius'])
            value=delta*delta*radius;assert value==F(child['absolute_cost_upper']);childbounds.append(value)
        assert point==hi and max(childbounds)==F(row['absolute_cost_upper']);uppers.append(max(childbounds))
        if idx%10==0:print(json.dumps(dict(parent_index=idx,new_facets=facetcount,whole_support_cells=totalcells,seconds=round(time.monotonic()-started,1))),flush=True)
    assert previous==F(old['beta_outer'][1]);upper=max(uppers);lower=F(data['lower']['genuine_compatible_loss_lower']);nominal=F(prior['nominal_compatible_lower'])
    assert lower<upper<=F(21,10**9)<F(prior['absolute_cost_upper']) and upper==F(out['absolute_cost_upper'])
    assert lower==F(out['genuine_compatible_lower']) and nominal==F(out['nominal_compatible_lower'])
    assert upper/lower==F(out['union_gap_ratio']) and upper/nominal==F(out['nominal_gap_ratio'])
    assert F(out['old_B27_upper'])==F(prior['absolute_cost_upper']) and 1-upper/F(prior['absolute_cost_upper'])==F(out['relative_upper_reduction'])
    result=dict(source_sha256=hashlib.sha256(raw).hexdigest(),independent_precision_digits=110,original_parent_bands_accounted_for=430,
        independently_checked_subbands=children_count,new_support_polygons=polygons,new_support_facets=facetcount,
        independent_persistence_cells=totalcells,independent_raw_quadratic_cells=rawcells,
        inherited_target_exclusion_and_lower_proofs_reused_by_hash=True,full_unknown_beta_and_persistence_coverage_verified=True,
        same_twelve_reports_actual_error_union_positive_target_and_initial_states=True,grid_interpolation_used_as_proof=False,
        sharpness_or_joint_attainability_of_polygon_vertices_claimed=False,
        absolute_cost_upper=str(upper),genuine_compatible_lower=str(lower),union_gap_ratio=str(upper/lower),nominal_gap_ratio=str(upper/nominal))
    (HERE/'artifact_check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(complete=True,upper=float(upper),gap=float(upper/lower),new_facets=facetcount,persistence_cells=totalcells,raw_quadratic_cells=rawcells,seconds=round(time.monotonic()-started,1))))
if __name__=='__main__':main()
