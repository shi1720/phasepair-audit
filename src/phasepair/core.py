from dataclasses import dataclass
import numpy as np

@dataclass
class PhaseLaw:
    """Mixture of uniform intervals [lo,hi] and atoms, with probabilities w."""
    lo: np.ndarray
    hi: np.ndarray
    w: np.ndarray

    def __post_init__(self):
        self.lo, self.hi, self.w = map(lambda a: np.asarray(a, float), (self.lo,self.hi,self.w))
        if not (self.lo.shape == self.hi.shape == self.w.shape): raise ValueError('shape mismatch')
        if np.any(self.hi < self.lo) or np.any(self.w < 0) or not np.isclose(self.w.sum(),1):
            raise ValueError('invalid probability law')

    @property
    def mean(self): return float(np.dot(self.w,(self.lo+self.hi)/2))

    def cdf(self, x, side='right'):
        x=np.asarray(x,float)
        a,b=self.lo,self.hi
        span=b-a; atoms=span <= 1e-11
        frac=np.clip((x[...,None]-a)/np.where(atoms,1,span),0,1)
        # Numerical tolerance is used only to identify coincident atom locations.
        delta=x[...,None]-a
        frac[...,atoms]=(delta[...,atoms] > 1e-10) if side=='left' else (delta[...,atoms] >= -1e-10)
        return np.sum(frac*self.w,axis=-1)

    def scaled(self, scale):
        if scale <= 0: raise ValueError('positive scale required')
        return PhaseLaw(self.lo*scale,self.hi*scale,self.w.copy())

@dataclass
class Profile:
    """Area in mm^2 at physical z knots in mm; zero outside support."""
    z: np.ndarray
    area: np.ndarray
    kind: str = 'linear'

    def __post_init__(self):
        self.z=np.asarray(self.z,float); self.area=np.asarray(self.area,float)
        if self.kind not in ('linear','slab'): raise ValueError('unknown profile')
        expected=len(self.z) if self.kind=='linear' else len(self.z)-1
        if len(self.area)!=expected or np.any(np.diff(self.z)<=0) or np.any(self.area<0):
            raise ValueError('invalid area profile')
        if self.kind=='linear' and (self.area[0]!=0 or self.area[-1]!=0):
            raise ValueError('linear support endpoints must be zero')

    @property
    def volume(self):
        if self.kind=='slab': return float(np.dot(np.diff(self.z),self.area))
        return float(np.sum(np.diff(self.z)*(self.area[:-1]+self.area[1:])/2))

    def evaluate(self,z):
        z=np.asarray(z,float)
        if self.kind=='linear': return np.interp(z,self.z,self.area,left=0,right=0)
        i=np.searchsorted(self.z,z,side='right')-1
        return np.where((i>=0)&(i<len(self.area)),self.area[np.clip(i,0,len(self.area)-1)],0)

    def sample_volume(self,h,phase):
        if h<=0: raise ValueError('spacing must be positive')
        phases=np.asarray(phase,float)
        k=np.arange(int(np.floor(self.z[0]/h))-1,int(np.ceil(self.z[-1]/h))+2)
        return h*np.sum(self.evaluate(phases[...,None]+h*k),axis=-1)

    def phase_law(self,h):
        if h<=0: raise ValueError('spacing must be positive')
        # Every slope/discontinuity event occurs at a source knot modulo h.
        r=np.mod(self.z,h)
        r[np.isclose(r,h,rtol=0,atol=1e-9)]=0
        r[np.isclose(r,0,rtol=0,atol=1e-9)]=0
        raw=np.sort(np.r_[0,r,h]); keep=np.r_[True,np.diff(raw)>1e-9]
        knots=raw[keep]
        if knots[-1] < h-1e-9: knots=np.r_[knots,h]
        weights=np.diff(knots)/h
        if self.kind=='linear':
            y=self.sample_volume(h,knots)/self.volume
            lo=np.minimum(y[:-1],y[1:]);hi=np.maximum(y[:-1],y[1:])
            same=hi-lo<1e-11;lo[same]=hi[same]=(lo[same]+hi[same])/2
        else:
            y=self.sample_volume(h,(knots[:-1]+knots[1:])/2)/self.volume;lo=y.copy();hi=y.copy()
        return PhaseLaw(lo,hi,weights)


def independent_growth(x,y,q=1.2):
    """P(Y>=qX) for independent positive-baseline laws; exact linear integration."""
    if q<=0 or np.min(x.lo)<=0: raise ValueError('positive baseline and threshold required')
    ans=0.
    for a,b,wx in zip(x.lo,x.hi,x.w):
        for c,d,wy in zip(y.lo,y.hi,y.w):
            if b-a<=1e-11:
                if d-c<=1e-11: p=float(c>=q*a-1e-10)
                else:p=float(np.clip((d-q*a)/(d-c),0,1))
            elif d-c<=1e-11:
                p=float(np.clip((c/q-a)/(b-a),0,1))
            else:
                cuts=np.unique(np.r_[a,b,np.clip([c/q,d/q],a,b)])
                v=np.clip((d-q*cuts)/(d-c),0,1)
                p=float(np.sum(np.diff(cuts)*(v[1:]+v[:-1])/2)/(b-a))
            ans+=wx*wy*p
    return float(np.clip(ans,0,1))


def worst_growth(x,y,q=1.2):
    """Sharp supremum of P(Y>=qX), over all couplings with the supplied marginals.

    Ordered transport / Makarov bound: 1-sup_t (F_Y(qt-)-F_X(t-)).
    Both one-sided limits are checked at every piecewise-linear CDF breakpoint.
    """
    if q<=0 or np.min(x.lo)<=0: raise ValueError('positive baseline and threshold required')
    knots=np.unique(np.r_[x.lo,x.hi,y.lo/q,y.hi/q])
    deficiency=max(0.,float(np.max(y.cdf(q*knots,'left')-x.cdf(knots,'left'))),
                   float(np.max(y.cdf(q*knots,'right')-x.cdf(knots,'right'))))
    return float(np.clip(1-deficiency,0,1))
