# Independent floating-point sanity checks of the lemmas used by the rigorous scripts.
import numpy as np, math
from flint import arb, acb
from rigorous_common import h_upper, h_float, near_critical_constants, SIGMA
from part1_small_alpha import tail_constant
rng=np.random.default_rng(1)
# 1. tail lemma: |h - |eta|^(-2/3)| * |eta|^2 <= C(R) for |eta| >= R
for R in [8,16,100]:
    C=float(tail_constant(R).upper()); worst=0
    for _ in range(20000):
        r=R*math.exp(rng.exponential(1.0)); th=rng.uniform(0,2*math.pi); e=r*complex(math.cos(th),math.sin(th))
        worst=max(worst,abs(h_float(e)-abs(e)**(-2/3))*abs(e)**2)
    print('tail lemma R=%g: C=%.4f, sampled max %.4f, ok=%s'%(R,C,worst,worst<=C))
# 2. near-critical bound: h <= 2K1/sqrt|d| + K3 for |d| <= dmax
for dm in [1e-4,1e-3,0.01,0.05]:
    K1,K3=near_critical_constants(arb(dm)); K1=float(K1.upper()); K3=float(K3.upper()); worst=-1
    for _ in range(20000):
        d=dm*math.sqrt(rng.uniform())*complex(*np.exp(1j*rng.uniform(0,2*math.pi)).view(float)) if False else dm*rng.uniform()**2*np.exp(1j*rng.uniform(0,2*np.pi))
        e=1j*float(SIGMA.mid())+d
        worst=max(worst,h_float(e)-(2*K1/math.sqrt(abs(d))+K3))
    print('near-critical dmax=%g: max(h - bound) = %.3e (must be <= 0)'%(dm,worst))
# 3. h_upper encloses h_float
bad=0
for _ in range(3000):
    c=complex(rng.normal(0,2),rng.normal(0,2)); rad=10**rng.uniform(-8,-2)
    ball=acb(arb(c.real,rad),arb(c.imag,rad)); U=h_upper(ball)
    if U is None: continue
    for _ in range(5):
        p=c+complex(rng.uniform(-rad,rad),rng.uniform(-rad,rad))
        if h_float(p)>float(U.upper())+1e-12: bad+=1
print('h_upper enclosure violations:',bad)
