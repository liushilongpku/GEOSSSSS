#!/usr/bin/env python3
"""
Dual-porosity (N=2) Mandel analytical solution — Mehrabian & Abousleiman (2014).

Dimensionless formulation following the paper exactly (eqs. 25-45). Effective
coefficients from Appendix A (eq. A25). Key point: in the dimensionless form
Sum_j S*_jj = 1, which makes eps_zz spatially uniform (eq. 37); the stress BC
(eq. 41) then closes f(t). Robust complex Laplace inversion via mpmath de Hoog.

Outputs p_m, p_f, sigma_zz at the center and u_z at the top, normalized for
comparison to digitized Fig 5c (pressure) and Fig 5d (stress).
"""
from pathlib import Path

import numpy as np
import mpmath as mp

mp.mp.dps = 50

# ----- GOM shale, double-porosity: matrix(1) + macrofracture(3), microfracture absent -----
Ks  = 27.0e9
Kf  = 1.744e9
mu  = 1.0e-3
nu  = 0.22
K1, phi1, k1, v1 = 1.1e9,  0.14, 4.93e-21, 0.97   # matrix
K2, phi2, k2, v2 = 2.25e7, 0.95, 4.93e-15, 0.03   # macrofracture
a  = 0.03
b  = 0.03
Gamma12 = 0.0    # Table 2: matrix<->macrofracture direct exchange = 0 for the N=2 case

# ----- effective coefficients (Appendix A) -----
alpha1 = 1.0 - K1/Ks
alpha2 = 1.0 - K2/Ks
Kbar   = 1.0/(v1/K1 + v2/K2)
abar   = np.array([Kbar*v1*alpha1/K1, Kbar*v2*alpha2/K2])
G      = 3.0*Kbar*(1.0-2.0*nu)/(2.0*(1.0+nu))
cm     = abar*(1.0-2.0*nu)/(2.0*G*(1.0-nu))            # eq.19
B1     = (1.0/K1-1.0/Ks)/(1.0/K1-1.0/Ks+phi1*(1.0/Kf-1.0/Ks))
B2     = (1.0/K2-1.0/Ks)/(1.0/K2-1.0/Ks+phi2*(1.0/Kf-1.0/Ks))
# eq.8 Biot modulus from compliance a_ij via constant-strain Legendre transform:
#   1/Mbar_ij = a_(i+1,j+1) - abar_i abar_j / Kbar   (note the MINUS sign)
# diagonal a_ii = v_i alpha_i/(B_i K_i); off-diagonal a_23 = 0 (Appendix A, eq.A13)
invM = np.array([[v1*alpha1/(B1*K1) - abar[0]*abar[0]/Kbar, -abar[0]*abar[1]/Kbar],
                 [-abar[0]*abar[1]/Kbar,                     v2*alpha2/(B2*K2) - abar[1]*abar[1]/Kbar]])
Sbar  = invM + np.outer(abar, cm)                      # eq.22
kappa = np.array([v1*k1/mu, v2*k2/mu])
Gam   = np.array([[Gamma12, -Gamma12], [-Gamma12, Gamma12]])

trS  = Sbar[0,0]+Sbar[1,1]
sumK = kappa[0]+kappa[1]
t0   = a*a*trS/sumK

# ----- dimensionless matrices (eq.25) : S*=S/trS, kappa*=kappa/sumK, Gamma*=Gamma a^2/sumK -----
Ss = Sbar/trS
ks = kappa/sumK
Gs = Gam*a*a/sumK
cms = cm/trS                       # c*_m = c_m/trS (so eps = sum c*_m p* + f)
# gamma = 1/(2G trS); the stress-BC constant kappa_c = gamma (1-2nu)^2/(1-nu)
kappa_c = (1.0-2.0*nu)**2/((1.0-nu)*2.0*G*trS)

def _laplace(sstar):
    """Dimensionless Laplace solution at s* -> (p1(0), p2(0), sigma_zz(0), uz_top)."""
    s = mp.mpc(sstar)
    A = mp.matrix([[s*Ss[0,0]+Gs[0,0], s*Ss[0,1]+Gs[0,1]],
                   [s*Ss[1,0]+Gs[1,0], s*Ss[1,1]+Gs[1,1]]])
    M = mp.matrix([[A[0,0]/ks[0], A[0,1]/ks[0]],
                   [A[1,0]/ks[1], A[1,1]/ks[1]]])
    tr = M[0,0]+M[1,1]; det = M[0,0]*M[1,1]-M[0,1]*M[1,0]
    disc = mp.sqrt(tr*tr-4*det)
    lam = [(tr+disc)/2, (tr-disc)/2]
    # eigenvectors normalized chi_1i = 1
    X = mp.zeros(2,2)
    for j,lm in enumerate(lam):
        if abs(M[0,1])>1e-300:
            v0, v1_ = M[0,1], lm-M[0,0]
        else:
            v0, v1_ = lm-M[1,1], M[1,0]
        X[0,j]=1.0; X[1,j]=v1_/v0
    psi = -s*(A**-1)*mp.matrix([abar[0],abar[1]])      # Psi* (eq.28)
    sq  = [mp.sqrt(lam[0]), mp.sqrt(lam[1])]
    coshv = [mp.cosh(sq[0]), mp.cosh(sq[1])]
    tanhR = [mp.tanh(sq[0])/sq[0], mp.tanh(sq[1])/sq[1]]
    Xinv = X**-1
    c = Xinv*psi                                        # (X^{-1} Psi)_j = Det[X_psi,j]/DetX
    # g_i(x*) = psi_i - sum_j X_ij c_j cosh(x* sq_j)/cosh(sq_j);  p*_i = f * g_i
    g0   = [psi[i] - (X[i,0]*c[0]/coshv[0] + X[i,1]*c[1]/coshv[1]) for i in range(2)]   # x*=0
    gint = [psi[i] - (X[i,0]*c[0]*tanhR[0] + X[i,1]*c[1]*tanhR[1]) for i in range(2)]   # int_0^1
    # stress BC (eq.41): int_0^1 sigma_zz dx* = -Pc/s,  sigma_zz = 2G f/(1-2nu) - (..)*sum abar_i p*_i
    #   -> f = -Pc / [ s * (2G/(1-2nu)) * (1 - kappa_c * sum abar_i gint_i) ]     (Pc=1)
    Q = abar[0]*gint[0] + abar[1]*gint[1]
    brace = (2.0*G/(1.0-2.0*nu))*(1.0 - kappa_c*Q)
    f = -1.0/(s*brace)
    p1 = f*g0[0]; p2 = f*g0[1]
    # sigma_zz at center (eq.34): 2G f/(1-2nu) - [(1-2nu)/((1-nu) sumS)] * sum abar_i p*_i(0)
    # The coupling coeff is 2G*kappa_c/(1-2nu) = (1-2nu)/((1-nu) sumS); kappa_c (with the
    # extra (1-2nu)) is correct only for the integrated stress BC that closes f.
    sig = 2.0*G*f/(1.0-2.0*nu) - (2.0*G*kappa_c/(1.0-2.0*nu))*(abar[0]*p1+abar[1]*p2)
    # u_z(top) = b * eps_zz, eps_zz uniform = f (1-nu)/(1-2nu)  (eq.37)
    uz = b*f*(1.0-nu)/(1.0-2.0*nu)
    return p1, p2, sig, uz

def invert(comp, tarr_phys):
    out=np.zeros(len(tarr_phys))
    fn=lambda s: _laplace(s)[comp]
    for i,t in enumerate(tarr_phys):
        if t<=0: out[i]=np.nan; continue
        out[i]=float(mp.invertlaplace(fn, t/t0, method='dehoog', dps=30))
    return out

# Directory holding the digitized paper curves for this validation case.
AN_DIR = str(Path(__file__).resolve().parent.parent / "analitical_result")

def solve(tau_min=1e-5, tau_max=1e6, n=70):
    """Evaluate the analytical solution over a dimensionless-time sweep.

    Returns a dict with tau and the p0+-normalized matrix/fracture pressures and
    sig0-normalized center stress.  tau spans [tau_min, tau_max]; default reaches 1e6.
    """
    tau=np.logspace(np.log10(tau_min), np.log10(tau_max), n)
    tt=tau*t0
    pm=invert(0,tt); pf=invert(1,tt); sg=invert(2,tt)
    pm0=invert(0,np.array([1e-7*t0]))[0]
    pf0=invert(1,np.array([1e-7*t0]))[0]
    # pressure normalized by undrained p0+ (Fig5c); stress normalized by Pc=1 and
    # sign-flipped to positive compression (Fig5d is sigma/Pc starting at ~1).
    return dict(tau=tau, pm=pm/pm0, pf=pf/pf0, sig=-sg,
                pm0=pm0, pf0=pf0)

def plot_vs_digitized(sol=None, out=None, tau_max=1e6):
    """Plot the analytical solution against the digitized Fig5c (pressure) and
    Fig5d (stress).  Saves to `out` (default analitical_result/analytical_vs_digitized_full.png).
    The dimensionless-time axis extends to tau_max (default 1e6)."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    if sol is None:
        sol=solve(tau_max=tau_max)
    if out is None:
        out=f"{AN_DIR}/analytical_vs_digitized_full.png"
    am=np.loadtxt(f"{AN_DIR}/fig5c_primary_analitical.csv",delimiter=',',skiprows=1)
    af=np.loadtxt(f"{AN_DIR}/fig5c_secondary_analitical.csv",delimiter=',',skiprows=1)
    ad=np.loadtxt(f"{AN_DIR}/fig5d_analitical.csv",delimiter=',',skiprows=1)
    tau=sol['tau']
    fig,ax=plt.subplots(1,3,figsize=(18,5))
    ax[0].semilogx(am[:,0],am[:,1],'k-',lw=2,label='digitized Fig5c (paper)')
    ax[0].semilogx(tau,sol['pm'],'r.-',ms=4,label='analytical')
    ax[0].set_title('MATRIX (primary) pressure'); ax[0].set_ylabel(r'$p/p_0^+$')
    ax[1].semilogx(af[:,0],af[:,1],'k-',lw=2,label='digitized Fig5c (paper)')
    ax[1].semilogx(tau,sol['pf'],'b.-',ms=4,label='analytical')
    ax[1].set_title('FRACTURE (secondary) pressure')
    ax[2].semilogx(ad[:,0],ad[:,1],'k-',lw=2,label='digitized Fig5d (paper)')
    ax[2].semilogx(tau,sol['sig'],'g.-',ms=4,label='analytical')
    ax[2].set_title(r'STRESS $\sigma_{zz}$ (center)'); ax[2].set_ylabel(r'$\sigma/\sigma_0$')
    for a in ax:
        a.set_xlabel(r'$\tau = t/t_0$,  $t_0$=%.2f s'%t0); a.grid(alpha=.3)
        a.legend(); a.set_xlim(1e-5, tau_max)
    plt.tight_layout(); plt.savefig(out, dpi=120); plt.close()
    print("saved %s"%out)
    return out

if __name__=='__main__':
    import sys
    print("t0=%.4f s  abar=[%.4f,%.4f]  Kbar=%.4e  G=%.4e  trS=%.4e"%(t0,abar[0],abar[1],Kbar,G,trS))
    TAU_MAX=1e6
    sol=solve(tau_max=TAU_MAX)
    am=np.loadtxt(f"{AN_DIR}/fig5c_primary_analitical.csv",delimiter=',',skiprows=1)
    af=np.loadtxt(f"{AN_DIR}/fig5c_secondary_analitical.csv",delimiter=',',skiprows=1)
    ad=np.loadtxt(f"{AN_DIR}/fig5d_analitical.csv",delimiter=',',skiprows=1)
    print("\n tau        p_m/p0 fig5c_m | p_f/p0 fig5c_f | sig/s0 fig5d")
    for tq in [1e-3,1e-2,3e-2,0.1,0.3,1,3,10,30,100,300,1000,3000,1e4,1e5,1e6]:
        i=np.argmin(np.abs(sol['tau']-tq))
        ia=np.argmin(np.abs(am[:,0]-tq)); ib=np.argmin(np.abs(af[:,0]-tq)); ic=np.argmin(np.abs(ad[:,0]-tq))
        print("%9.1e  %6.3f %6.3f | %6.3f %6.3f | %6.3f %6.3f"%(
            tq, sol['pm'][i], am[ia,1], sol['pf'][i], af[ib,1], sol['sig'][i], ad[ic,1]))
    plot_vs_digitized(sol, tau_max=TAU_MAX)
