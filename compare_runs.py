"""Compare frozen and repeated PhasePair results; no third-party packages required."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
keys=['primary.csv','primary_summary.json','slab.csv','slab_summary.json',
      'half_cap.csv','half_cap_summary.json','repaired.csv','repaired_summary.json',
      'known_changes_and_thresholds.csv','relative_shift.csv','cross_spacing.csv',
      'exploratory_resolution_curve.csv','verification_grid.csv','verification_transport.csv']
report={name:(ROOT/'results'/name).read_bytes()==(ROOT/'replication/results'/name).read_bytes() for name in keys}
a=json.loads((ROOT/'results/verification_summary.json').read_text())
b=json.loads((ROOT/'replication/results/verification_summary.json').read_text())
for x in (a,b):x.pop('runtime_seconds',None)
report['verification_summary_non_timing_identical']=a==b
print(json.dumps(report,indent=2))
if not all(report.values()):raise SystemExit('At least one repeated result differs')
