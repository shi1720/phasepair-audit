from pathlib import Path
import json,sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from run_study import construct
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':240,
 'axes.grid':True,'grid.alpha':.16,'legend.frameon':False})
BLUE='#245a81';RED='#b44f3b';TEAL='#21897e';GRAY='#a6afb6'
records={r['id']:r for r in json.loads((ROOT/'data/profiles.json').read_text())}
df=pd.read_csv(ROOT/'results/primary.csv');summary=json.loads((ROOT/'results/primary_summary.json').read_text())
selected=df[df.spacing_mm==3].sort_values('worst',ascending=False).iloc[0]
r=records[selected.id];p=construct(r)

def save(fig,name):
 fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');fig.savefig(OUT/(name+'.png'),bbox_inches='tight');plt.close(fig)

fig=plt.figure(figsize=(7.1,4.8),layout='constrained');gs=fig.add_gridspec(2,2)
ax=fig.add_subplot(gs[0,0],projection='3d')
obj=json.loads((ROOT/'data/example_contours.json').read_text())['target']
for i,pts in enumerate(obj['LPS_contour_points']):
 a=np.array(pts);ax.add_collection3d(Poly3DCollection([a],alpha=.28,facecolor=BLUE,edgecolor=BLUE,lw=.6))
a=np.concatenate(obj['LPS_contour_points']); lo=a.min(0);hi=a.max(0)
ax.set_xlim(lo[0]-1,hi[0]+1);ax.set_ylim(lo[1]-1,hi[1]+1);ax.set_zlim(lo[2]-1,hi[2]+1)
ax.set_box_aspect(np.maximum(hi-lo,1));ax.view_init(23,-50)
ax.set_xlabel('x (mm)',labelpad=-4);ax.set_ylabel('y (mm)',labelpad=-4);ax.set_zlabel('z (mm)',labelpad=-3)
ax.tick_params(labelsize=6,pad=-2)
ax.set_title('a  Public target contours: case 141',loc='left',pad=8)
ax=fig.add_subplot(gs[0,1]);ax.plot(p.z,p.area,color=BLUE,lw=2);ax.scatter(r['z'],r['area'],s=18,color=BLUE,zorder=3)
ax.set(xlabel='Relative axial position (mm)',ylabel=r'Cross-sectional area (mm$^2$)')
ax.set_title(f'b  Area reconstruction ({p.volume:.2f} mm³)',loc='left')
ax=fig.add_subplot(gs[1,0]);u=np.linspace(0,1,600);v=p.sample_volume(3,3*u)/p.volume
ax.plot(u,v,color=BLUE,lw=2);ax.axhline(1,color=GRAY,ls='--',lw=1)
ax.set(xlabel='Grid origin / spacing',ylabel='Sampled / integrated volume',ylim=(.75,1.4));ax.set_title('c  One unbiased volume distribution',loc='left')
ax=fig.add_subplot(gs[1,1]);u=(np.arange(400)+.5)/400;v=p.sample_volume(3,3*u)
event=(v[:,None]>=1.2*v[None,:])
ax.imshow(event,extent=(0,1,0,1),origin='lower',aspect='auto',cmap=ListedColormap(['#edf1f3',RED]),interpolation='nearest')
ax.plot([0,1],[0,1],color=TEAL,lw=1.5,label='Identical origins')
ax.set(xlabel='Baseline origin / spacing',ylabel='Follow-up origin / spacing');ax.set_title('d  Apparent ≥20% increase (red)',loc='left');ax.grid(False)
save(fig,'fig1_mechanism')

fig,axs=plt.subplots(1,2,figsize=(7.1,2.7),layout='constrained')
ax=axs[0];x=np.arange(4);w=.31
for i,(key,color,label) in enumerate([('independent',BLUE,'Independent origins'),('worst',RED,'Worst admissible dependence')]):
 y=np.array([s[key]*100 for s in summary]);cis=np.array([s[key+'_ci'] for s in summary])*100
 ax.bar(x+(i-.5)*w,y,w,color=color,label=label,zorder=3)
 ax.errorbar(x+(i-.5)*w,y,yerr=np.stack([y-cis[:,0],cis[:,1]-y]),fmt='none',color=color,capsize=2,lw=1)
ax.set_xticks(x,['1','1.5','2','3']);ax.set(xlabel='Slice spacing (mm)',ylabel='Mean false-growth probability (%)',ylim=(0,2.85))
ax.legend(fontsize=7,loc='upper left');ax.set_title('a  Primary cohort (n = 217)',loc='left')
ax=axs[1];d=df[df.spacing_mm==3];s=d.susceptible
ax.scatter(d.loc[~s,'volume_mm3'],100*(d.loc[~s,'maximum_ratio']/d.loc[~s,'minimum_ratio']-1),s=16,c=GRAY,alpha=.8)
ax.scatter(d.loc[s,'volume_mm3'],100*(d.loc[s,'maximum_ratio']/d.loc[s,'minimum_ratio']-1),s=23,c=RED,label='13 susceptible objects')
ax.axhline(20,color=RED,ls='--',lw=1);ax.set_xscale('log');ax.set(xlabel=r'Integrated volume (mm$^3$)',ylabel='Largest apparent increase (%)',ylim=(-1,68))
ax.legend(fontsize=7,loc='upper right');ax.set_title('b  Origin envelope at 3 mm',loc='left');save(fig,'fig2_cohort')

fig,axs=plt.subplots(1,2,figsize=(7.1,2.65),layout='constrained')
s=pd.read_csv(ROOT/'results/relative_shift.csv');ax=axs[0]
for h,color in [(2,BLUE),(3,RED)]:
 a=s[s.spacing_mm==h].groupby('relative_shift').false_growth.mean()
 ax.plot(a.index,a*100,color=color,label=f'{h} mm',lw=2)
 ax.axhline(df[df.spacing_mm==h].independent.mean()*100,color=color,ls=':',lw=1)
ax.set(xlabel='Fixed relative origin shift / spacing',ylabel='Mean false-growth probability (%)',ylim=(-.04,1.15))
ax.legend(fontsize=8);ax.set_title('a  Fixed shifts preserve uniform marginals',loc='left')
ax=axs[1];curve=pd.read_csv(ROOT/'results/exploratory_resolution_curve.csv')
ax.plot(curve.spacing_mm,curve.independent*100,color=BLUE,lw=2,label='Independent')
ax.plot(curve.spacing_mm,curve.worst*100,color=RED,lw=2,label='Worst dependence')
for h in [1,1.5,2,3]:ax.axvline(h,color=GRAY,lw=.6,ls=':')
ax.set(xlabel='Slice spacing (mm)',ylabel='Mean false-growth probability (%)')
ax.legend(fontsize=8,loc='upper left');ax.set_title('b  Exploratory spacing sweep',loc='left');save(fig,'fig3_dependence')

fig,axs=plt.subplots(1,2,figsize=(7.1,2.75),layout='constrained');ax=axs[0]
labels=['Linear','Constant slab','Half-length caps','Repaired polygons'];tags=['primary','slab','half_cap','repaired'];vals=[]
for tag in tags:vals.append(json.loads((ROOT/f'results/{tag}_summary.json').read_text())[-1])
x=np.arange(4);w=.31
for i,(key,color,label) in enumerate([('independent',BLUE,'Independent'),('worst',RED,'Worst dependence')]):
 ax.bar(x+(i-.5)*w,[v[key]*100 for v in vals],w,color=color,label=label,zorder=3)
ax.set_xticks(x,labels,rotation=18,ha='right');ax.set(ylabel='Mean false-growth probability (%)',ylim=(0,3.8))
ax.legend(fontsize=7,loc='upper left');ax.set_title('a  Reconstruction sensitivity at 3 mm',loc='left')
ax=axs[1];d=pd.read_csv(ROOT/'results/known_changes_and_thresholds.csv');d=d[d.threshold.isna()]
for h,color in [(1,TEAL),(2,BLUE),(3,RED)]:
 a=d[d.spacing_mm==h].groupby('true_ratio').independent_growth.mean()
 ax.plot(a.index,a*100,'o-',ms=4,color=color,label=f'{h} mm')
ax.axvline(1.2,color=GRAY,ls='--',lw=1);ax.set(xlabel='Known synthetic volume ratio',ylabel='Probability of a growth call (%)',ylim=(-3,103))
ax.legend(fontsize=8);ax.set_title('b  Controlled cross-sectional area scaling',loc='left');save(fig,'fig4_sensitivity')
print('Four figures generated from stored results.')
