# float prototype: G(alpha) = (L0 - L)/sqrt(a) via scaled tau-formulation
import mpmath as mp
mp.mp.dps=20
sig=2/(3*mp.sqrt(3))
def h(eta):
    rs=mp.polyroots([1,0,1,-eta],maxsteps=200,extraprec=60)
    return sum(1/abs(3*r**2+1) for r in rs)
def G(al):
    al=mp.mpf(al); b=mp.sqrt(1-4*al**2/27); eps=(4*al/27)/(1+b)
    def F(tau):
        e0=1j*tau/(1-1j*al*tau/2)
        return (abs(e0)**(-mp.mpf(2)/3)-h(e0+eps))/(1+al**2*tau**2/4)
    # critical tau
    tc=((1j*sig-eps)/(1j-sig*al/2-1j*eps*al/2))
    tc=mp.re(tc)
    return 2*mp.quad(F,[0,tc,1,10,100,mp.inf])
for al in [0,1e-4,1e-2,0.1,1,2.598]:
    print(al,G(al))
