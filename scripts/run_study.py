from pathlib import Path
import json,time,platform,sys
import numpy as np,pandas as pd
from phasepair import Profile,independent_growth,worst_growth
ROOT=Path(__file__).resolve().parents[1]
SPACINGS=[1.,1.5,2.,3.]
SEED=1720406

def construct(r,kind='linear'):
    z=np.array(r['z']);a=np.array(r['area']);d0=z[1]-z[0];d1=z[-1]-z[-2]
    if kind=='slab':return Profile(np.r_[z[0]-d0/2,(z[:-1]+z[1:])/2,z[-1]+d1/2],a,'slab')
    cap=.5 if kind=='half_cap' else 1.
    return Profile(np.r_[z[0]-cap*d0,z,z[-1]+cap*d1],np.r_[0,a,0])

def shifted_growth(p,h,delta,g=1.,q=1.2):
    # Exact one-dimensional phase length for a fixed relative origin shift.
    r=np.mod(p.z/h,1);knots=np.sort(np.r_[0,1,r,np.mod(r-delta,1)])
    knots=knots[np.r_[True,np.diff(knots)>1e-9]]
    if p.kind=='slab':
        u=(knots[:-1]+knots[1:])/2
        val=g*p.sample_volume(h,h*np.mod(u+delta,1))-q*p.sample_volume(h,h*u)
        return float(np.sum(np.diff(knots)*(val>=-1e-9)))
    u=knots
    vals=g*p.sample_volume(h,h*np.mod(u+delta,1))-q*p.sample_volume(h,h*u)
    a,b=vals[:-1],vals[1:];w=np.diff(knots)
    ans=0.
    for x,y,weight in zip(a,b,w):
        if min(x,y)>=-1e-9:frac=1.
        elif max(x,y)<-1e-9:frac=0.
        else: frac=max(x,y)/abs(y-x)
        ans+=weight*np.clip(frac,0,1)
    return float(ans)

def run(records,kind='linear',tag='primary',full=False):
    rows=[];known=[];shift=[];start=time.perf_counter()
    for j,r in enumerate(records):
        p=construct(r,kind)
        for h in SPACINGS:
            law=p.phase_law(h)
            row=dict(id=r['id'],representation=kind,volume_mm3=p.volume,nplanes=len(r['z']),
                     native_min_mm=r['native_min'],spacing_mm=h,phase_mean_ratio=law.mean,
                     minimum_ratio=float(law.lo.min()),maximum_ratio=float(law.hi.max()),
                     components=len(law.w),non_detection_possible=law.lo.min()<=0)
            if row['non_detection_possible']:
                row.update(independent=np.nan,worst=np.nan,susceptible=True)
            else:
                row.update(independent=independent_growth(law,law),worst=worst_growth(law,law),
                           susceptible=law.hi.max() >=1.2*law.lo.min()+1e-9)
            rows.append(row)
            if full and law.lo.min()>0:
                for g in [1.,1.1,1.2,1.3,1.5]:
                    known.append(dict(id=r['id'],spacing_mm=h,true_ratio=g,
                         independent_growth=independent_growth(law,law.scaled(g)),
                         worst_growth=worst_growth(law,law.scaled(g))))
                for q in [1.1,1.15,1.25,1.3]:
                    known.append(dict(id=r['id'],spacing_mm=h,true_ratio=1.,threshold=q,
                         independent_growth=independent_growth(law,law,q),
                         worst_growth=worst_growth(law,law,q)))
                for delta in np.linspace(0,1,41):
                    shift.append(dict(id=r['id'],spacing_mm=h,relative_shift=float(delta),
                                     false_growth=shifted_growth(p,h,float(delta))))
        if (j+1)%50==0:print(tag,j+1,flush=True)
    df=pd.DataFrame(rows);df.to_csv(ROOT/f'results/{tag}.csv',index=False)
    if full:
        pd.DataFrame(known).to_csv(ROOT/'results/known_changes_and_thresholds.csv',index=False)
        pd.DataFrame(shift).to_csv(ROOT/'results/relative_shift.csv',index=False)
    summary=[];rng=np.random.default_rng(SEED)
    for h,d in df.groupby('spacing_mm'):
        vals=d[['independent','worst','susceptible']].astype(float).to_numpy()
        means=vals.mean(0);ix=rng.integers(0,len(d),size=(10000,len(d)))
        boots=vals[ix].mean(1);ci=np.quantile(boots,[.025,.975],axis=0)
        result=dict(spacing_mm=h,n=len(d),mean_phase_bias=float(np.max(abs(d.phase_mean_ratio-1))),
                    non_detection_cases=int(d.non_detection_possible.sum()),
                    max_phase_ratio=float((d.maximum_ratio/d.minimum_ratio).max()))
        for i,key in enumerate(['independent','worst','susceptible']):
            result[key]=float(means[i]);result[key+'_ci']=[float(v) for v in ci[:,i]]
        result['susceptible_count']=int(d.susceptible.sum());summary.append(result)
    (ROOT/f'results/{tag}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(tag,'seconds',round(time.perf_counter()-start,3),json.dumps(summary),flush=True)
    return df

if __name__=='__main__':
    records=json.loads((ROOT/'data/profiles.json').read_text())
    (ROOT/'results').mkdir(exist_ok=True)
    run(records,full=True)
    run(records,'slab','slab')
    run(records,'half_cap','half_cap')
    repaired=json.loads((ROOT/'data/profiles_repaired.json').read_text())
    run(repaired,'linear','repaired')
    vols=np.array([construct(r).volume for r in records])
    info=dict(n=len(records),volume_quantiles_mm3=np.quantile(vols,[0,.25,.5,.75,1]).tolist(),
       native_spacing_counts=pd.Series([round(r['native_min'],3) for r in records]).value_counts().to_dict(),
       python=sys.version,platform=platform.platform(),seed=SEED,bootstrap_resamples=10000)
    (ROOT/'results/environment.json').write_text(json.dumps(info,indent=2)+'\n')
