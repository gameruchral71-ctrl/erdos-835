import mpmath as mp
mp.mp.dps=25
def params(theta):
    s=mp.sin(theta/2); b=mp.cos(theta/2); a=3*mp.cbrt(s**2/4); return a,b
def integrand_sum(a,b,t):
    w=mp.expj(t)
    rs=mp.polyroots([1,0,a,b-w],maxsteps=400,extraprec=80)
    return sum(1/abs(3*r**2+a) for r in rs)
def L(theta):
    a,b=params(theta)
    def arc(t1,t2):
        D=t2-t1
        return mp.quad(lambda ph: D*mp.sin(2*ph)*integrand_sum(a,b,t1+D*mp.sin(ph)**2),[0,mp.pi/4,mp.pi/2])
    if theta==0:
        return arc(0,2*mp.pi)
    return arc(-theta/2,theta/2)+arc(theta/2,2*mp.pi-theta/2)
L0=2**(mp.mpf(1)/3)*mp.beta(mp.mpf(1)/6,0.5)
print('L0',L0)
for th in ['1e-9','1e-7','1e-5','1e-3','1e-2','0.1','0.5','1','2','3',mp.pi]:
    th=mp.mpf(th); a,b=params(th)
    v=L(th); print(mp.nstr(th,6), 'a=',mp.nstr(a,6),'L=',mp.nstr(v,15),'deficit/sqrt(a)=',mp.nstr((L0-v)/mp.sqrt(a),10))
