#!/usr/bin/env python3
"""Independent P3 replay: interval caps for dual potentials and Q.
No project verifier is imported; formulas are re-derived from the source note.
"""
import json, math, hashlib, time, os
from pathlib import Path
from fractions import Fraction as F
from flint import arb, acb, ctx
ctx.dps=90

_HERE = Path(__file__).resolve()
ROOT = Path(__import__('os').environ.get(
    'P3_REFERENCE', str(_HERE.parent / 'reference.json')))
CAPX=F(-4626,1000); CAPV=F(-4480,1000); CAPZ=F(-90,1000)
R=F(3,5)

def frac(s):
    if isinstance(s,F): return s
    if isinstance(s,int): return F(s)
    if '/' in s:
        a,b=s.split('/'); return F(int(a),int(b))
    return F(s)
def A(q):
    q=frac(q); return arb(q.numerator)/q.denominator
def C(q):
    q=frac(q); return acb(A(q),arb(0))
def flabs(q):
    q=frac(q)
    if q==0: return arb(0)
    x=A(q)
    return abs(x).log()
def ffun(q):
    q=frac(q)
    if q==0: return arb(0)
    return A(q)*flabs(q)-A(q)

def phi(q,a,b):
    q,a,b=map(frac,(q,a,b)); L=b-a
    return (ffun(q-a)-ffun(q-b))/A(L)

def gfun(q):
    q=frac(q)
    if q==0: return arb(0)
    x=A(q)
    return A(F(1,2))*x*x*flabs(q)-A(F(3,4))*x*x

def parse_list(D,key,complex_=False):
    out=[]
    for q in D[key]:
        if complex_:
            out.append((frac(q['re']),frac(q['im']),frac(q['w'])))
        else:
            out.append((frac(q['a']),frac(q['b']),frac(q['w'])))
    return out

D=json.loads(ROOT.read_text())
X=parse_list(D,'x'); Y=parse_list(D,'y'); ZR=parse_list(D,'z_real'); ZP=parse_list(D,'z_poisson',True)

def upper_max(vals):
    """Return a candidate whose outward upper endpoint bounds all values."""
    z=vals[0]
    for q in vals[1:]:
        if q.upper() > z.upper():
            z=q
    return z

def max2(vals):
    return upper_max(vals)

def circle_pot_real(q,cr,ci):
    q=frac(q); cr=frac(cr); ci=frac(ci)
    if abs(q) <= R:
        u=A(R*R)-A(q)*A(cr)
        v=A(q)*A(ci)
        return (u*u+v*v).log()/2 - A(R).log()
    u=A(q-cr); v=A(ci)
    return (u*u+v*v).log()/2

def d_one(q,a,b):
    q,a,b=map(frac,(q,a,b));
    if q==0: return arb(2).log()
    u1=A(2+q*a); u2=A(2+q*b)
    # integral log(2+q t)dt / (b-a)
    def FF(u): return u*u.log()-u
    return (FF(u2)-FF(u1))/(A(q)*A(b-a))

def dpot(q, segs):
    return sum((A(w)*d_one(q,a,b) for a,b,w in segs),arb(0))

def umax_real(q0,q1,segs):
    # Per component endpoint maxima; valid since a uniform log potential has one minimum.
    ans=arb(0)
    for a,b,w in segs:
        p0=phi(q0,a,b); p1=phi(q1,a,b)
        ans += A(w)*upper_max((p0,p1))
    return ans

def uz_real_endpoint(q,segs):
    return sum((A(w)*phi(q,a,b) for a,b,w in segs),arb(0))

def uz_poisson_real(q0,q1):
    # endpoint maxima for each Poisson circle; split caller at R.
    ans=arb(0)
    for cr,ci,w in ZP:
        p0=circle_pot_real(q0,cr,ci); p1=circle_pot_real(q1,cr,ci)
        # convex quadratic in q on each branch; add r if crossing
        vals=[p0,p1]
        if q0<R<q1: vals.append(circle_pot_real(R,cr,ci))
        ans += A(w)*max2(vals)
    return ans

def uz_real_qmax(q0,q1):
    # absolute endpoint maxima; componentwise safe
    return umax_real(q0,q1,ZR)+uz_poisson_real(q0,q1)

def uv_qmax(q0,q1): return umax_real(q0,q1,Y)
def ux_qmax(q0,q1): return umax_real(q0,q1,X)

def fx_upper(q0,q1):
    # q0,q1 in [0,1], q0<=q1. The +log q term is at q1, -D at q0.
    return 2*ux_qmax(q0,q1) - dpot(q0,Y) + A(3*F(1,4))*flabs(q1)

def fv_upper(q0,q1):
    return 2*uv_qmax(q0,q1) + uz_real_qmax(q0,q1) - dpot(q0,X) + A(F(1,4))*flabs(q1)

def fz_real_upper(q0,q1):
    # componentwise endpoint bounds on chain real path.
    uz=uz_real_qmax(q0,q1)
    uv=uv_qmax(q0,q1)
    # -5/4 log|q| largest at smallest modulus endpoint
    qm=min(abs(q0),abs(q1))
    return 2*uz + uv - A(F(5,4))*flabs(qm)

def sqrt_frac(x): return math.sqrt(float(x))

def max_circle_logdist(cr,ci,t0,t1):
    # Certified endpoint/cross-product test.  The only interior maximum is
    # the antipodal point -c/|c|.  For an upper semicircle, its angle lies in
    # [t0,t1] iff ci<=0 and both endpoint cross products are nonnegative.
    # All signs are decided by Arb; an interval touching zero falls back to
    # the full-circle bound, so no floating atan2 decision is used.
    from flint import fmpq
    rr=A(R); vals=[]
    for tt in (t0,t1):
        tq=fmpq(tt.numerator,tt.denominator)
        co,sn=arb.cos_pi_fmpq(tq),arb.sin_pi_fmpq(tq)
        dx=rr*co-A(cr); dy=rr*sn-A(ci)
        vals.append((dx*dx+dy*dy).log()/2)
    ciA=A(ci); crA=A(cr)
    if ci <= 0:
        q0=fmpq(t0.numerator,t0.denominator); q1=fmpq(t1.numerator,t1.denominator)
        c0,s0=arb.cos_pi_fmpq(q0),arb.sin_pi_fmpq(q0)
        c1,s1=arb.cos_pi_fmpq(q1),arb.sin_pi_fmpq(q1)
        x0=(-ciA)*c0 + crA*s0
        x1=(-crA)*s1 + ciA*c1
        if x0.lower() >= 0 and x1.lower() >= 0:
            cc=(crA*crA+ciA*ciA).sqrt(); vals.append((rr+cc).log())
        elif x0.upper() >= 0 and x1.upper() >= 0:
            # Ambiguous endpoint: the outward fallback is certified.
            cc=(crA*crA+ciA*ciA).sqrt(); vals.append((rr+cc).log())
    return max2(vals)


def uniform_pot_circle(t,segs):
    # Potential of a real uniform segment at z=r exp(i pi t), via a complex antiderivative.
    from flint import fmpq
    tq=fmpq(t.numerator,t.denominator)
    co,sn=arb.cos_pi_fmpq(tq),arb.sin_pi_fmpq(tq)
    z=acb(A(R)*co,A(R)*sn)
    ans=arb(0)
    for a,b,w in segs:
        # int_a^b log|z-u| du: alpha=z, beta=-1
        ans += A(w)*acb_linear_integral(z,(-1,0),a,b)/A(b-a)
    return ans


def trig_bounds_pi(t0,t1):
    # rigorous real bounds for sin/cos(pi*t) on 0<=t<=1 from endpoints and extrema.
    from flint import fmpq
    q0=fmpq(t0.numerator,t0.denominator); q1=fmpq(t1.numerator,t1.denominator)
    c0=arb.cos_pi_fmpq(q0); c1=arb.cos_pi_fmpq(q1)
    s0=arb.sin_pi_fmpq(q0); s1=arb.sin_pi_fmpq(q1)
    clo=min(c0.lower(),c1.lower()); chi=max(c0.upper(),c1.upper())
    slo=min(s0.lower(),s1.lower()); shi=max(s0.upper(),s1.upper())
    if t0 <= F(1,2) <= t1: shi=max(shi,arb(1).upper())
    if t0 <= F(0) <= t1: slo=min(slo,arb(0).lower())
    if t0 <= F(1) <= t1: slo=min(slo,arb(0).lower())
    def make(lo,hi): return arb((lo+hi)/2,(hi-lo)/2)
    return make(clo,chi),make(slo,shi)

def poisson_sum_upper(t0,t1):
    co,sn=trig_bounds_pi(t0,t1)
    rr=A(R); ans=arb(0)
    for cr,ci,w in ZP:
        # d^2=r^2+|c|^2-2r(cr cos+ci sin)
        d2=rr*rr+A(cr*cr+ci*ci)-2*rr*(A(cr)*co+A(ci)*sn)
        ans += A(w)*d2.log()/2
    return ans

def fz_circle_upper(t0,t1):
    # z=r exp(i theta), upper semicircle theta/pi=t.
    # Negative real z reference decreases in theta -> theta left; positive y ref increases -> theta right.
    uzr=uz_real_endpoint(-R,ZR) if t0==F(1) else uniform_pot_circle(t0,ZR)  # decreases with angle
    uvy=uniform_pot_circle(t1,Y)  # positive real potentials increase with angle
    uzp=arb(0)
    for cr,ci,w in ZP:
        uzp += A(w)*max_circle_logdist(cr,ci,t0,t1)
    return 2*(uzr+uzp)+uvy-A(F(5,4))*flabs(R)

# Q exact formulas ------------------------------------------
def rect_log_int(a,b,c,d):
    # integral_a^b integral_c^d log|u-v| du dv, using F''=log|.|.
    a,b,c,d=map(frac,(a,b,c,d))
    return gfun(a-d)+gfun(b-c)-gfun(a-c)-gfun(b-d)

def I_uniform(P,Q):
    out=arb(0)
    for a,b,w in P:
        for c,d,v in Q:
            out += A(w*v)*rect_log_int(a,b,c,d)/A((b-a)*(d-c))
    return out

def acb_linear_integral(alpha,beta,a,b):
    # int_a^b log|alpha+beta*x| dx, via Re((z log z-z)/beta)
    a,b=map(frac,(a,b));
    if isinstance(alpha,acb): alpha_c=alpha
    elif isinstance(alpha,tuple): alpha_c=acb(A(alpha[0]),A(alpha[1]))
    else: alpha_c=acb(A(alpha),arb(0))
    if isinstance(beta,tuple): beta_c=acb(A(beta[0]),A(beta[1]))
    else: beta_c=acb(A(beta),arb(0))
    aa=alpha_c+beta_c*A(a); bb=alpha_c+beta_c*A(b)
    beta=beta_c
    def L(z): return z*(z.log()-1)
    return ((L(bb)-L(aa))/beta).real

def I_uniform_poisson(P,Q):
    # P real uniform; Q list (center, mass), using Poisson potential branch; split at r.
    out=arb(0)
    for a,b,w in P:
        # z-real P intervals are <= -r, x/y may span r
        cuts=[a]
        if a<R<b: cuts.append(R)
        cuts.append(b)
        for aa,bb in zip(cuts,cuts[1:]):
            mid=(aa+bb)/2
            for cr,ci,v in Q:
                if abs(mid) <= R:
                    # log|r^2 - x*conj(c)| - log r; alpha=r^2, beta=-conj(c)
                    integ=acb_linear_integral(R*R,(-cr,ci),aa,bb) - A(bb-aa)*A(R).log()
                else:
                    integ=acb_linear_integral((-cr,-ci),(1,0),aa,bb)
                out += A(w*v)*A(bb-aa)/A(b-a)*integ/A(bb-aa)
    return out

def I_poisson_poisson(P,Q):
    out=arb(0)
    for cr,ci,w in P:
        for dr,di,v in Q:
            # potential of circle centered c at point d inside disk
            val=((A(R*R)-A(dr)*A(cr))**2 + (A(dr)*A(ci))**2).log()/2 - A(R).log()
            # formula with complex centers: |r^2-d*conj(c)|
            # above assumes d real; correct general expression below
            # z = r^2 - d * conj(c), with d=dr+i di
            z=C(R*R)-acb(A(dr),A(di))*acb(A(cr),-A(ci))
            val=(z.real*z.real+z.imag*z.imag).log()/2-A(R).log()
            out += A(w*v)*val
    return out

def D_uniform(P,Q,N=420):
    # log(2+x y)=log2 + sum_{n>=1}(-1)^(n+1)(xy/2)^n/n.
    out=A(2).log()*sum((A(w*v) for a,b,w in P for c,d,v in Q),arb(0))
    for a,b,w in P:
        for c,d,v in Q:
            sx=arb(0); sy=arb(0)
            for n in range(1,N+1):
                mx=(A(b)**(n+1)-A(a)**(n+1))/(A(n+1)*A(b-a))
                my=(A(d)**(n+1)-A(c)**(n+1))/(A(n+1)*A(d-c))
                term=mx*my/(A(n)*A(2)**n)
                sx += term if n%2 else -term
            pair=A(w*v)
            out += pair*sx
            rho=A(max(abs(a),abs(b))*max(abs(c),abs(d)))/A(2)
            tail=2*(rho**(N+1))/A(N+1)
            # directed remainder; make it symmetric outward by adding a radius
            out += arb(0, pair*tail)
    return out

def compute_Q():
    Ixx=I_uniform(X,X)
    Iyy=I_uniform(Y,Y)
    Izr=I_uniform(ZR,ZR)
    Izp=I_uniform_poisson(ZR,ZP)
    Ipp=I_poisson_poisson(ZP,ZP)
    Izz=Izr+2*Izp+Ipp
    Iyz=I_uniform(Y,ZR)+I_uniform_poisson(Y,ZP)
    Dxy=D_uniform(X,Y)
    Q=Ixx+Iyy+Izz+Iyz-Dxy
    return Q,dict(Ixx=Ixx,Iyy=Iyy,Izr=Izr,Izp=Izp,Ipp=Ipp,Iyz=Iyz,Dxy=Dxy)

# adaptive leaves ------------------------------------------
def cert_cap(kind,q0,q1):
    if kind=='x': return fx_upper(q0,q1)
    if kind=='v': return fv_upper(q0,q1)
    if kind=='zn': return fz_real_upper(q0,q1)
    if kind=='zp': return fz_real_upper(q0,q1)
    if kind=='zc': return fz_circle_upper(q0,q1)
    raise ValueError(kind)

def bisect_real(kind,target,maxdepth=30):
    todo=[(F(0),F(1),0)] if kind in ('x','v') else ([(F(-2),-R,0)] if kind=='zn' else [(R,F(1),0)])
    accepted=[]; maxd=0; rejects=0
    while todo:
        a,b,d=todo.pop()
        val=cert_cap(kind,a,b)
        maxd=max(maxd,d)
        if val.upper() < A(target): accepted.append((a,b,d,val)); continue
        if d>=maxdepth:
            raise RuntimeError(('not certified',kind,a,b,d,val,A(target)))
        mid=(a+b)/2
        todo.append((mid,b,d+1));todo.append((a,mid,d+1));rejects+=1
    return accepted,maxd,rejects

def bisect_circle(target,maxdepth=30):
    todo=[(F(0),F(1),0)];accepted=[];maxd=0;rejects=0
    while todo:
        a,b,d=todo.pop(); val=fz_circle_upper(a,b);maxd=max(maxd,d)
        if val.upper() < A(target): accepted.append((a,b,d,val));continue
        if d>=maxdepth: raise RuntimeError(('not certified','zc',a,b,d,val,A(target)))
        m=(a+b)/2;todo.append((m,b,d+1));todo.append((a,m,d+1));rejects+=1
    return accepted,maxd,rejects

def main():
    t=time.time()
    print('mass checks',sum((w for a,b,w in X),F(0)),sum((w for a,b,w in Y),F(0)),sum((w for a,b,w in ZR+[(0,0,w) for cr,ci,w in ZP]),F(0)))
    print('Q computing...')
    Q,parts=compute_Q(); print('Q',Q); print('parts',parts)
    print('caps adaptive...')
    allacc={}
    for k,target in [('x',CAPX),('v',CAPV),('zn',CAPZ),('zp',CAPZ)]:
        aa,md,rj=bisect_real(k,target); allacc[k]=aa; print(k,len(aa),md,rj,'worst',max((v for *_,v in aa),key=lambda z:z.upper()))
    aa,md,rj=bisect_circle(CAPZ);allacc['zc']=aa;print('zc',len(aa),md,rj,'worst',max((v for *_,v in aa),key=lambda z:z.upper()))
    print('total',sum(map(len,allacc.values())),'time',time.time()-t)
    import pickle
    # JSON compact report, keep intervals as strings
    out={'dps':ctx.dps,'Q':str(Q),'caps':{k:{'leaves':len(v),'max_depth':max(d for a,b,d,x in v),'max_upper':str(max((x for a,b,d,x in v),key=lambda z:z.upper()))} for k,v in allacc.items()},'total_leaves':sum(map(len,allacc.values())),'seconds':time.time()-t}
    output = Path(os.environ.get(
        'P3_OUTPUT', '/tmp/p3_dual_reaudit_20261005/results/replay.json'))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2))

if __name__=='__main__': main()
