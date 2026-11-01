from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font=FontProperties(fname='C:/Windows/Fonts/simsun.ttc')
plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False})
d=json.loads((HERE/'results/case_results.json').read_text('utf-8'))
fig,ax=plt.subplots(1,2,figsize=(7.1,2.8),layout='constrained')
for key,color,marker,label in [('first','#285777','o','首次反馈'),('buy','#a15e2f','s','购买重复消息')]:
    f=[r for r in d['frontier'] if r['permission']==key and r['H']<=4096 and r['B_certificate'] is not None]
    ax[0].plot([r['H'] for r in f],[r['B_certificate'] for r in f],color=color,marker=marker,ms=4,label=label,lw=1.3)
    f=[r for r in d['validation'] if r['permission']==key]
    ax[1].plot([r['theta'] for r in f],[r['windows_mean'] for r in f],color=color,marker=marker,ms=4,label=label,lw=1.3)
for a in ax:
    a.grid(axis='y',color='.9',lw=.6);a.legend(prop=font,frameon=False,fontsize=8)
ax[0].set_xlabel('完成期限 H（窗口数）',fontproperties=font);ax[0].set_ylabel('充分准备金',fontproperties=font)
ax[0].set_xticks([1024,2048,3072,4096]);ax[0].tick_params(axis='x',labelsize=8)
ax[1].set_xlabel('真实参数 θ',fontproperties=font);ax[1].set_ylabel('平均支付窗口数（对数刻度）',fontproperties=font)
ax[1].set_yscale('log',base=2);ax[1].set_ylim(64,5200);ax[1].set_yticks([64,256,1024,4096],labels=['64','256','1024','4096'])
ax[1].axhline(4096,color='.5',ls='--',lw=.8);ax[1].text(-.95,4400,'期限 H=4096',fontproperties=font,fontsize=8)
images=HERE/'figures';images.mkdir(exist_ok=True)
fig.savefig(images/'financing_case.pdf');fig.savefig(HERE/'results/financing_case.png',dpi=180);plt.close(fig)
