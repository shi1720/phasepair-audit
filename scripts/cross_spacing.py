from pathlib import Path
import json,numpy as np,pandas as pd
from phasepair import independent_growth,worst_growth
from run_study import construct
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for r in json.loads((ROOT/'data/profiles.json').read_text()):
 p=construct(r);laws={h:p.phase_law(h) for h in [1.,1.5,2.,3.]}
 for h0,x in laws.items():
  for h1,y in laws.items():
   rows.append(dict(id=r['id'],baseline_spacing_mm=h0,followup_spacing_mm=h1,
    independent=independent_growth(x,y),worst=worst_growth(x,y)))
pd.DataFrame(rows).to_csv(ROOT/'results/cross_spacing.csv',index=False)
print(pd.DataFrame(rows).groupby(['baseline_spacing_mm','followup_spacing_mm'])[['independent','worst']].mean())
