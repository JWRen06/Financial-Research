"""Focused numerical study for the rewritten opportunity paper. All data are synthetic."""
from pathlib import Path
import sys,importlib.util,json,math,time,hashlib,subprocess,os,platform
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]
if os.environ.get('OPPORTUNITY_EXPERIMENTS'):EXPS=Path(os.environ['OPPORTUNITY_EXPERIMENTS'])
elif (HERE/'local_rerun/experiments').exists():EXPS=HERE/'local_rerun/experiments'
else:EXPS=HERE/'experiments'
for p in [BASE/'工作区/04_第二第三篇共同研究/原始工作目录/experiments/_local_dependencies',BASE/'工作区/04_第二第三篇共同研究/原始工作目录/experiments/_plot_dependencies']:
    if p.exists():sys.path.insert(0,str(p))
import numpy as np, scipy
from scipy.optimize import linprog,minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
spec=importlib.util.spec_from_file_location('cap',EXPS/'application_contract_bridge.py')
cap=importlib.util.module_from_spec(spec);spec.loader.exec_module(cap)
OUT=Path(os.environ.get('OPPORTUNITY_OUTPUT',str(HERE/'focused_results')));OUT.mkdir(exist_ok=True,parents=True)
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1.3,'savefig.dpi':300,'pdf.fonttype':42,'ps.fonttype':42})
CN=FontProperties(fname='C:/Windows/Fonts/simsun.ttc',size=10)
start=time.perf_counter()
K=100;c=.8;rho=1.4;m=8;z0=.9
true_v=-np.expm1(-np.array([rho/K,z0]));true_w=np.array([K-m,m])
def revenue(q):return np.sum(true_w*(-np.expm1(np.asarray(q)[...,None]*np.log1p(-true_v))),axis=-1)
p=np.array([revenue(q) for q in [1,2,3]])
def robust(p,S,c,Q):
    options=[(0.,0.)]
    for lo,hi,f in [(0.,1.,cap.peak_h),(1.,2.,cap.peak_g),(2.,3.,cap.peak_h),(3.,Q,cap.peak_g)]:
        if hi<lo:continue
        options.append(f(p,S,c,lo,hi) if f==cap.peak_h else f(p,c,lo,hi))
    return max(options,key=lambda x:x[1])
rows=[]
for Q in [3.,5.,10.,20.,40.,80.,160.]:
    q,v=robust(p,K,c,Q);qu,u=cap.capped_upper(p,K,c,Q)
    rows.append({'Q':Q,'robust_action':q,'guaranteed_profit':v,'upper_profit':u,'true_profit_at_robust_action':float(revenue(q)-c*q),'optimistic_action':qu})

# Arbitrary continuous experimental intensities: finite-grid LP, independently refined.
# These LP results are numerical approximations, never used as rigorous guarantees.
arbitrary=[]
for design in [[1.,2.,3.],[.75,1.8,3.]]:
    readings=np.array([revenue(q) for q in design])
    for resolution in [501,2001,8001]:
        nodes=np.unique(np.r_[np.linspace(0.,1.,resolution),true_v,cap.capped_atoms(p,K)[0],cap.parameters(p)[4]])
        kernel=lambda q:-np.expm1(q*np.log1p(-np.minimum(nodes,1-1e-16)))
        A=np.array([kernel(t) for t in design]);A[:,nodes==1]=1
        def lp(q,sign=1.):
            response=kernel(q);response[nodes==1]=1.
            sol=linprog(sign*response,A_eq=A,b_eq=readings,A_ub=np.ones((1,len(nodes))),b_ub=[K],bounds=(0,None),method='highs')
            if not sol.success:raise RuntimeError(sol.message)
            return sol.fun*sign
        opt=minimize_scalar(lambda q:-(lp(q)-c*q),bounds=(3,12),method='bounded',options={'xatol':2e-6})
        q=float(opt.x);value=-float(opt.fun)
        upper20=lp(20.,-1)-c*20.
        arbitrary.append({'test_intensities':design,'grid_nodes':len(nodes),'approximate_robust_action':q,'approximate_guarantee':value,'approximate_upper_at_q20':upper20,'discretization_only_not_uniform_certificate':True})
        if design==[1.,2.,3.]:
            exact=robust(p,K,c,20.)
            assert abs(value-exact[1])<.002

native=json.loads((EXPS/'application_contract_bridge.json').read_text('utf-8'))
groups=native['native_action_groups'];paths=native['native_action_paths']
assert native['native_paths']==7200 and all(x['joint_successes']==x['paths'] for x in groups)
summary=[]
for S in [100,400,1600]:
    for Q in [20,80]:
        data=[x for x in paths if x['cap']==S and x['Q']==Q]
        summary.append({'K':S,'Q':Q,'paths':len(data),'median_records':float(np.median([x['records'] for x in data])),'median_action':float(np.median([x['action'] for x in data])),'max_actual_loss':max(x['actual_best_robust']-x['actual_guarantee'] for x in data),'joint_successes':sum(x['joint_success'] for x in data)})
q=np.linspace(0,80,1000);low,high=cap.capped_values(p,K,q);low[0]=high[0]=0
qstar,vstar=robust(p,K,c,80)
fig,ax=plt.subplots(1,2,figsize=(9.0,3.3))
ax[0].fill_between(q,low-c*q,high-c*q,color='.88',label='试验结果允许的利润范围')
ax[0].plot(q,revenue(q)-c*q,color='black',ls='--',label='生成数据的真实曲线')
ax[0].plot(q,low-c*q,color='#254e70',label='最坏情形利润')
ax[0].scatter([qstar],[vstar],color='#254e70',s=24,zorder=5)
ax[0].set_ylim(-10,20);ax[0].set_xlim(0,80);ax[0].set_xlabel('连续投入强度 q',fontproperties=CN);ax[0].set_ylabel('每窗口期望利润',fontproperties=CN)
ax[0].legend(prop=CN,frameon=False,loc='lower left',fontsize=8)
shown=[x for x in rows if x['Q']<=80]
ax[1].plot([x['Q'] for x in shown],[x['upper_profit'] for x in shown],'o-',color='black',label='最高可能利润')
ax[1].plot([x['Q'] for x in shown],[x['guaranteed_profit'] for x in shown],'s-',color='#254e70',label='最优保守保证')
ax[1].set_xlim(0,82);ax[1].set_xticks([0,20,40,60,80]);ax[1].set_xlabel('允许的最大投入 Q',fontproperties=CN);ax[1].set_ylabel('每窗口期望利润',fontproperties=CN);ax[1].legend(prop=CN,frameon=False)
for a in ax:a.grid(axis='y',alpha=.2,linewidth=.4)
fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(OUT/('resource_decision.'+ext),bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(9.,3.3))
for kind,label,style in [('fixed_cap','固定 K=100','o-'),('cap_proportional_to_Q','K 随 Q 扩大（K=10Q）','s--')]:
    rr=[x for x in native['finite_scale_necessary_records'] if x['regime']==kind]
    ax[0].plot([x['Q'] for x in rr],[x['necessary_continuous_pilot_records'] for x in rr],style,label=label,color='black' if kind=='fixed_cap' else '#254e70',markersize=4)
ax[0].set(xscale='log',yscale='log');ax[0].set_xlabel('允许的最大投入 Q',fontproperties=CN);ax[0].set_ylabel('解析检验下界（记录条数）',fontproperties=CN);ax[0].legend(prop=CN,frameon=False)
for Q,col,shift in [(20,'#254e70',-.08),(80,'.5',.08)]:
    rr=[x for x in summary if x['Q']==Q]
    ax[1].plot(np.arange(3)+shift,[x['median_records'] for x in rr],'o',color=col,label=f'Q={Q}')
ax[1].set_xticks(range(3),['100','400','1600']);ax[1].set_ylim(0,470000);ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[1].set_xlabel('已知机会数 K',fontproperties=CN);ax[1].set_ylabel('行动程序的记录中位数',fontproperties=CN);ax[1].legend(frameon=False)
for a in ax:a.grid(axis='y',alpha=.2,linewidth=.4)
fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(OUT/('measurement_and_stopping.'+ext),bbox_inches='tight')
plt.close(fig)
report={'synthetic':True,'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'baseline':{'K':K,'fast_count':m,'slow_speed':rho/K,'fast_speed':z0,'cost':c,'readings':p.tolist(),'robust_action':qstar,'guaranteed_profit':vstar,'minimum_mass':cap.parameters(p)[-1]},'decision_table':rows,'arbitrary_continuous_design':arbitrary,'native_action_summary':summary,'native_total_paths':len(paths),'native_total_success':sum(x['joint_success'] for x in paths),'joint_group_success_lower_bound':min(x['simultaneous_one_sided_95_lower'] for x in groups),'seconds':time.perf_counter()-start,'data_source':'fresh local_rerun/experiments/application_contract_bridge.json'}
(OUT/'study_results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),'utf-8')
print(json.dumps(report,ensure_ascii=True),flush=True)
