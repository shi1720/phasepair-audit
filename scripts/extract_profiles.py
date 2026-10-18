"""Download the cited contour archive and extract reproducible, audited area profiles."""
from pathlib import Path
import io,json,zipfile,hashlib,sys
import requests,numpy as np,pandas as pd
from shapely.geometry import Polygon
from shapely import make_valid
ROOT=Path(__file__).resolve().parents[1]
URL='https://www.cancerimagingarchive.net/wp-content/uploads/Vestibular-Schwannoma-SEG-contours-Mar-2021.zip'

def extract(archive,repair=False):
    zf=zipfile.ZipFile(archive); profiles=[]; audit=[]
    for name in sorted(zf.namelist()):
        if not name.endswith('_t1/contours.json'):continue
        objects=json.loads(zf.read(name));pid=name.split('/')[-2]
        targets=[r for r in objects if r['structure_name'].lower() in ('tv','an','tvt1','rt an','tumour')]
        reason=[]; planes={};invalid=0; maxz=0
        if len(targets)!=1:reason.append('no_unambiguous_target')
        else:
            for contour in targets[0]['LPS_contour_points']:
                a=np.array(contour);maxz=max(maxz,float(np.ptp(a[:,2])))
                poly=Polygon(a[:,:2]); invalid+=int(not poly.is_valid)
                if repair:poly=make_valid(poly)
                z=round(float(a[:,2].mean()),5)
                planes.setdefault(z,[]).append(poly)
            if maxz>0.001:reason.append('nonaxial')
            if invalid and not repair:reason.append('invalid_polygon')
        if not reason:
            zs=sorted(planes);area=[]
            for z in zs:
                ps=planes[z];geom=ps[0]
                for p in ps[1:]:geom=geom.symmetric_difference(p)
                area.append(geom.area)
            zs=np.array(zs); zs-=zs[0]
            profiles.append(dict(id=pid,z=zs.tolist(),area=area,native_min=float(np.diff(zs).min()),
                                 native_max=float(np.diff(zs).max()),invalid_repaired=invalid))
        audit.append(dict(id=pid,included=not reason,reason=';'.join(reason),invalid_polygons=invalid,
                          max_within_contour_z_range_mm=maxz,structure=targets[0]['structure_name'] if len(targets)==1 else ''))
    return profiles,audit

if __name__=='__main__':
    raw=ROOT/'data/raw'; raw.mkdir(exist_ok=True,parents=True)
    path=raw/'contours.zip'
    if not path.exists():
        r=requests.get(URL,timeout=180);r.raise_for_status();path.write_bytes(r.content)
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    expected='318261d120722f45be2bdf01f7970a2a5cc4a220f2d9efc478ff13544a1eb5aa'
    if digest != expected:
        raise ValueError('Source archive hash differs from the frozen March 2021 export')
    for repair in (False,True):
        p,a=extract(path,repair)
        suffix='_repaired' if repair else ''
        (ROOT/f'data/profiles{suffix}.json').write_text(json.dumps(p,indent=2)+'\n')
        pd.DataFrame(a).to_csv(ROOT/f'data/cohort_audit{suffix}.csv',index=False)
        print(suffix or 'primary','included',len(p),'of',len(a))
    (ROOT/'data/provenance.json').write_text(json.dumps(dict(url=URL,sha256=digest,bytes=path.stat().st_size,
       dataset_doi='10.7937/TCIA.9YTJ-5Q73',paper_doi='10.1038/s41597-021-01064-w',
       retrieval_date='2026-10-17',license='CC BY 4.0',source='March 2021 public JSON contour export'),indent=2)+'\n')
