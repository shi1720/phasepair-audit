"""Independent checks of cohort computations and finite transport optima."""
from pathlib import Path
import json,time
import numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import qmc
from phasepair import PhaseLaw,independent_growth,worst_growth
from run_study import construct,shifted_growth
ROOT=Path(__file__).resolve().parents[1]
records={r['id']:r for r in json.loads((ROOT/'data/profiles.json').read_text())}
df=pd.read_csv(ROOT/'results/primary.csv')
checks=[];start=time.perf_counter()
# Deterministic stratified phase integration, independent of the exact law engine.
for row in df.itertuples():
    p=construct(records[row.id]);u=(np.arange(4096)+.5)/4096
    samples=p.sample_volume(row.spacing_mm,row.spacing_mm*u)/p.volume
    ref=float(np.mean([np.mean(samples>=1.2*x) for x in samples]))
    checks.append(dict(id=row.id,spacing_mm=row.spacing_mm,exact=row.independent,
                       deterministic_sweep=ref,absolute_difference=abs(ref-row.independent)))
pd.DataFrame(checks).to_csv(ROOT/'results/verification_grid.csv',index=False)
# 1024-phase assignment is a separately solved finite coupling problem.
rows=[]
for row in df[df.spacing_mm==3].sort_values('worst',ascending=False).head(15).itertuples():
    p=construct(records[row.id]);n=1024;u=(np.arange(n)+.5)/n
    v=p.sample_volume(3,3*u)/p.volume
    cost=-(v[None,:]>=1.2*v[:,None]).astype(float)
    a,b=linear_sum_assignment(cost)
    finite=float(-cost[a,b].mean())
    rows.append(dict(id=row.id,sharp_continuous=row.worst,assignment_1024=finite,
                     absolute_difference=abs(finite-row.worst)))
pd.DataFrame(rows).to_csv(ROOT/'results/verification_transport.csv',index=False)
# Monte Carlo on the highest-risk object (two independent origin draws).
r=records[df[df.spacing_mm==3].sort_values('worst',ascending=False).iloc[0].id]
p=construct(r);rng=np.random.default_rng(1720406);v=p.sample_volume(3,rng.uniform(0,3,(2,1000000)))
mc=float(np.mean(v[1]>=1.2*v[0]));law=p.phase_law(3)
# Slab integration and linear interpolation must agree on uniform source spacings.
slab_discrepancies=[abs(construct(r).volume-construct(r,'slab').volume)/construct(r).volume for r in records.values()]
summary=dict(cohort_grid_checks=len(checks),max_grid_error=max(r['absolute_difference'] for r in checks),
 mean_grid_error=float(np.mean([r['absolute_difference'] for r in checks])),
 transport_checks=len(rows),max_transport_error=max(r['absolute_difference'] for r in rows),
 monte_carlo_case=r['id'],monte_carlo_n=1000000,monte_carlo=mc,
 monte_carlo_exact=independent_growth(law,law),
 max_linear_slab_relative_volume_difference=max(slab_discrepancies),runtime_seconds=time.perf_counter()-start)
(ROOT/'results/verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
