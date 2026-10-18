import numpy as np
import pytest
from scipy.optimize import linprog
from phasepair import Profile,PhaseLaw,independent_growth,worst_growth


def test_rectangle_and_cavalieri_unbiasedness():
    p=Profile([0,2,5],[4,4],kind='slab')
    law=p.phase_law(3)
    assert p.volume==20
    assert law.mean==pytest.approx(1)
    assert sorted(law.lo)==pytest.approx([.6,1.2])
    assert independent_growth(law,law)==pytest.approx(2/9)
    assert worst_growth(law,law)==pytest.approx(1/3)


def test_triangle_exact_phase_distribution():
    p=Profile([0,1,2],[0,1,0])
    law=p.phase_law(3)
    assert law.mean==pytest.approx(1)
    assert min(law.lo)==0
    assert max(law.hi)==pytest.approx(3)
    # At spacing 1 the shifted hat sums to its exact integral.
    a=p.phase_law(1)
    assert np.allclose(a.lo,1)
    assert np.allclose(a.hi,1)


def test_independent_continuous_analytic():
    p=PhaseLaw([1],[2],[1])
    # Unit-square area y>=1.2*x equals integral_(1)^(5/3) (2-1.2x)dx.
    assert independent_growth(p,p)==pytest.approx(4/15)
    assert worst_growth(p,p)==pytest.approx(2/3)


def test_random_discrete_coupling_against_linear_program():
    rng=np.random.default_rng(1720406)
    for _ in range(100):
        a=np.sort(rng.integers(1,30,7)/10); b=np.sort(rng.integers(1,30,8)/10)
        wa=rng.dirichlet(np.ones(7));wb=rng.dirichlet(np.ones(8));q=float(rng.choice([1,1.2,1.5,2]))
        x=PhaseLaw(a,a,wa);y=PhaseLaw(b,b,wb)
        cost=-(b[None,:]>=q*a[:,None]-1e-10).astype(float)
        ae=np.zeros((15,56))
        for i in range(7):ae[i,i*8:(i+1)*8]=1
        for j in range(8):ae[7+j,j::8]=1
        opt=linprog(cost.ravel(),A_eq=ae,b_eq=np.r_[wa,wb],bounds=(0,None),method='highs')
        assert opt.success
        assert worst_growth(x,y,q)==pytest.approx(-opt.fun,abs=2e-9)
        direct=float(np.sum(-cost*wa[:,None]*wb[None,:]))
        assert independent_growth(x,y,q)==pytest.approx(direct,abs=1e-10)


def test_random_linear_against_monte_carlo():
    rng=np.random.default_rng(42)
    for _ in range(10):
        z=np.r_[0,np.cumsum(rng.uniform(.3,2,15))]
        a=np.r_[0,rng.uniform(1,10,14),0]
        p=Profile(z,a); h=float(rng.uniform(1,4));law=p.phase_law(h)
        assert law.mean==pytest.approx(1,abs=1e-10)
        phases=rng.uniform(0,h,(2,100000))
        v=p.sample_volume(h,phases)
        mc=float(np.mean(v[1]>=1.2*v[0]))
        exact=independent_growth(law,law)
        assert abs(mc-exact)<.006
        assert exact<=worst_growth(law,law)+1e-10


def test_zero_baseline_is_not_a_growth_ratio():
    law=PhaseLaw([0],[1],[1])
    with pytest.raises(ValueError):independent_growth(law,law)
    with pytest.raises(ValueError):worst_growth(law,law)
