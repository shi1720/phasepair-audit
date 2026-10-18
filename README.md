# PhasePair reproducibility artifact

**Shivam Gupta — Independent Researcher, Delhi, India**

This artifact audits paired slice-origin sensitivity of explicit contour-derived volume measurements. It contains code and numerical data, not a clinical prediction model or a journal manuscript. Existing Cavalieri sampling and fixed-marginal coupling results are credited in the research manuscript; this implementation specializes those methods to the stated audit.

## Frozen numerical findings

The primary cohort comprises 217 geometrically admissible public targets from 242 candidates. Repaired-polygon sensitivity retains 238. Under the primary area reconstruction, no target can cross a 20% numerical-growth threshold at 1 or 1.5 mm. At 3 mm, 13/217 are susceptible: the mean independent-origin probability is 0.442588%, compared with a mean object-specific worst-dependence bound of 1.422326%. These are conditional numerical results for unchanged digital objects, **not clinical false-positive rates**. Reconstruction materially affects the results.

## Contents

- `src/`, `tests/`, `scripts/`: frozen computation and verification code.
- `data/`: contour-derived areas, full inclusion/exclusion audits and one example geometry.
- `results/`: original complete numerical outputs, including secondary and exploratory analyses.
- `replication/results/`: separately executed results from 18 October 2026.
- `replication/logs/`, `replication/Provenance.json`: repeat logs, source hashes and extraction checks.
- `figures/`: original programmatic vector figures and previews.
- `compare_runs.py`: standalone comparison of original and repeated numerical results.
- `MANIFEST.json`: SHA-256 of every distributed file except the manifest itself.

The original and repeat are executions of the same deterministic analysis. They are not independent patient cohorts and are never pooled as clinical evidence. Repeated output files remain separate; runtime metadata are not required to match.

## Reproduce

Tested with Python 3.12.14. Use the pinned `requirements-lock.txt`:

```sh
python -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m pip install -e .
.venv/bin/python -m pytest -q
.venv/bin/python compare_runs.py
.venv/bin/python scripts/run_study.py
.venv/bin/python scripts/cross_spacing.py
.venv/bin/python scripts/resolution_curve.py
.venv/bin/python scripts/verify_cohort.py
.venv/bin/python scripts/make_figures.py
```

Run in a copy to preserve the frozen `results/` files. The calculations are offline using bundled profiles. `scripts/extract_profiles.py` can repeat extraction from the original contour ZIP; it checks the frozen SHA-256 and downloads about 16.7 MB if the ZIP is absent. It does not retrieve MRI intensities. Network availability at the original data provider is required for that optional download. The retained `data/provenance.json` describes the original acquisition; the repeat date is recorded separately.

## Data attribution and licensing

Shapey, J., Kujawa, A., Dorent, R., Wang, G., Bisdas, S., Dimitriadis, A., Grishchuck, D., Paddick, I., Kitchen, N., Bradford, R., Saeed, S., Ourselin, S., and Vercauteren, T. (2021). *Segmentation of Vestibular Schwannoma from Magnetic Resonance Imaging: An Open Annotated Dataset and Baseline Algorithm*, version 2 [dataset]. The Cancer Imaging Archive. https://doi.org/10.7937/TCIA.9YTJ-5Q73

Related source article: Shapey et al., *Scientific Data* 8, 286 (2021), https://doi.org/10.1038/s41597-021-01064-w.

Source and derived contour geometry retain **CC BY 4.0**, https://creativecommons.org/licenses/by/4.0/. The derived areas, profile translations and inclusion decisions are computational adaptations; they are not new clinical annotations. See `data/LICENSE.txt`. Original software is **MIT licensed**; see `LICENSE`. Source investigators retain credit for acquisition, delineation and the previously documented contour-conversion issues. The original ZIP hash is `318261d120722f45be2bdf01f7970a2a5cc4a220f2d9efc478ff13544a1eb5aa`.

## Numerical conventions

The implementation uses floating point. It analytically integrates the finite event functions subject to the documented tolerances; it does not provide interval-certified real arithmetic. Phase residues within 1e-9 mm are merged, intervals shorter than 1e-11 in normalized volume collapse to atoms, and atom comparisons use 1e-10 tolerance. Susceptibility requires a 1e-9 excess margin. The unrestricted dependence bound is a stress test, not a measured scanner coupling.

## Research record

The original protocol and computation were retained in a private author repository. This public artifact contains no earlier manuscript or conference submission text. A separate VS 2027 Seoul abstract-only poster proposal (FP-0006) is pending review; that does not establish journal publication or acceptance. No patient identities, private keys or registration credentials are included.
