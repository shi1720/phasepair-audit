# PhasePair numerical study protocol
Frozen before inspecting outcome rates, 17 October 2026.

Question: how much does dependence between two slice-grid origins change the probability of a spurious >=20% volumetric increase, while both single-scan volume distributions remain identical?

Design: contour-derived, unchanged digital objects; not a longitudinal clinical study. Public Vestibular-Schwannoma-SEG March 2021 contour archive, T1-coordinate target contours, one target per patient. Recognized target names: TV/tv, AN/an, TvT1, Rt AN, tumour. Exclude files without an unambiguous target, nonaxial contours (z range >0.001 mm), or invalid planar polygons. Preserve exclusions and reasons. Multiple coplanar polygons use even-odd fill. No selection by outcomes.

Primary representation: exact planar polygon areas; piecewise-linear area between native planes, zero at one additional native interval before/after support. Primary spacings: 1, 1.5, 2, 3 mm. Reference volume is the integral of this explicit area profile, not biological ground truth. Primary event: V_followup >=1.2 V_baseline with positive baseline; non-detection recorded separately.

Primary endpoints: per-object independent-origin false-growth probability; supremum over all origin couplings preserving uniform marginals; fraction of objects with any threshold crossing. Equal-origin null pairs as negative control. Compute piecewise-linear phase response exactly, rather than using a finite random phase sweep. Report equally weighted object means and patient bootstrap 95% intervals (10,000 resamples, seed 1720406); intervals describe resampling this cohort only.

Secondary: thresholds 10/15/25/30%; constant-slab reconstruction, half-interval terminal caps; optional invalid-polygon repair analyzed separately; known synthetic area scalings 1.1, 1.2, 1.3, 1.5; fixed relative phase shifts; cross-spacing grid pairs. No learned model, no train/test performance claims. Existing native-volume integration is the invariant comparator.

Verification: analytic rectangle/triangle phantoms; exact distribution mean equals profile integral; event integration against Monte Carlo; coupling formula against finite optimal assignment/transport; no-crossing and vanishing-volume cases; cohort completeness and numeric precision.

Scope: this isolates numerical sampling. It does not simulate MRI intensity acquisition, voxel point-spread functions, segmentation errors or clinical outcomes. Clinical false-positive rates cannot be estimated here. Makarov coupling bounds, Cavalieri unbiasedness, polygon integration and volume-preserving integration are established; contribution sought is the paired endpoint audit, exact implementation and VS contour benchmark.
