"""Exploratory spacing curve, declared separately from primary four spacings."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from phasepair import independent_growth,worst_growth
from run_study import construct
ROOT=Path(__file__).resolve().parents[1]
profiles=[construct(r) for r in json.loads((ROOT/'data/profiles.json').read_text())]
rows=[]
for h in np.round(np.arange(.75,4.001,.025),4):
    ps=[];ws=[];ss=[]
    for p in profiles:
        law=p.phase_law(h)
        ps.append(independent_growth(law,law));ws.append(worst_growth(law,law))
        ss.append(law.hi.max()>1.2*law.lo.min()+1e-9)
    rows.append(dict(spacing_mm=float(h),independent=float(np.mean(ps)),worst=float(np.mean(ws)),
                     susceptible_count=int(sum(ss))))
pd.DataFrame(rows).to_csv(ROOT/'results/exploratory_resolution_curve.csv',index=False)
print('resolution curve complete',len(rows))
