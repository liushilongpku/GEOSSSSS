#!/usr/bin/env python3
"""
Plot GEOS dual-porosity Mandel results against the validated analytical solution
(dpdp_mandel_analytical.py) and the digitized paper curves (Fig5c pressure,
Fig5d stress).

Usage:
  python plot_GEOS_vs_analytical.py [GEOS_OUTPUT_DIR] [TAU_MAX]

  GEOS_OUTPUT_DIR : directory containing pressure_matrix_history.hdf5,
                    pressure_fracture_history.hdf5 (default: /tmp/n2cs).
                    Pass "none" to plot analytical vs digitized only.
  TAU_MAX         : max dimensionless time on the axis (default 1e7).

Normalization matches the paper:
  pressure  p_i / p_i(0+)   (Skempton initial; p_m0=4.55e5, p_f0=4.88e5 Pa)
  stress    sigma / Pc
  time      tau = t / t0,  t0 = a^2 trace(S) / sum(kappa)  (~10.5 s)
"""
import sys, os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VALID_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "dpdkHMValidation"))
sys.path.insert(0, VALID_DIR)

# import the validated analytical solver (provides t0, solve())
import dpdp_mandel_analytical as an

# Skempton initial pressures used as the GEOS normalization (match the input deck IC)
PM0 = 4.55e5
PF0 = 4.88e5

def load_geos(out_dir):
    def load(fn):
        with h5py.File(os.path.join(out_dir, fn), 'r') as f:
            t = np.array(f['pressure Time'])[:, 0]
            p = np.array(f['pressure'])
            ec = np.array(f['pressure elementCenter'])[0]
        return t, p, ec
    t, pm, ec = load('pressure_matrix_history.hdf5')
    _, pf, _ = load('pressure_fracture_history.hdf5')
    x, z = ec[:, 0], ec[:, 2]
    # element nearest the center axis (x~0, mid-z)
    ix = np.argmin(np.abs(x - x.min()) + np.abs(z - 0.015))
    return t, pm[:, ix], pf[:, ix]

def main():
    geos_dir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/n2cs"
    tau_max = float(sys.argv[2]) if len(sys.argv) > 2 else 1e7

    t0 = an.t0
    print("t0 = %.3f s ; tau_max = %.1e" % (t0, tau_max))

    # ---- digitized paper curves ----
    am = np.loadtxt(os.path.join(HERE, "fig5c_primary_analitical.csv"), delimiter=',', skiprows=1)
    af = np.loadtxt(os.path.join(HERE, "fig5c_secondary_analitical.csv"), delimiter=',', skiprows=1)
    ad = np.loadtxt(os.path.join(HERE, "fig5d_analitical.csv"), delimiter=',', skiprows=1)

    # ---- analytical (validated), tau up to tau_max ----
    sol = an.solve(tau_min=1e-5, tau_max=tau_max, n=80)

    # ---- GEOS results (optional) ----
    geos = None
    if geos_dir.lower() != "none" and os.path.exists(os.path.join(geos_dir, 'pressure_matrix_history.hdf5')):
        tg, pmg, pfg = load_geos(geos_dir)
        m = tg > 0
        geos = dict(tau=tg[m] / t0, pm=pmg[m] / PM0, pf=pfg[m] / PF0)
        print("loaded GEOS from %s (tau up to %.1e)" % (geos_dir, geos['tau'].max()))
    else:
        print("no GEOS data (%s) -> analytical vs digitized only" % geos_dir)

    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    # matrix pressure
    ax[0].semilogx(am[:, 0], am[:, 1], 'k-', lw=2, label='digitized Fig5c (paper)')
    ax[0].semilogx(sol['tau'], sol['pm'], 'g-', lw=1.5, label='analytical (validated)')
    if geos: ax[0].semilogx(geos['tau'], geos['pm'], 'r.--', ms=4, label='GEOS')
    ax[0].set_title('MATRIX (primary) pressure'); ax[0].set_ylabel(r'$p/p_0^+$'); ax[0].set_ylim(0, 1.5)
    # fracture pressure
    ax[1].semilogx(af[:, 0], af[:, 1], 'k-', lw=2, label='digitized Fig5c (paper)')
    ax[1].semilogx(sol['tau'], sol['pf'], 'g-', lw=1.5, label='analytical')
    if geos: ax[1].semilogx(geos['tau'], geos['pf'], 'b.--', ms=4, label='GEOS')
    ax[1].set_title('FRACTURE (secondary) pressure')
    # stress (analytical vs digitized; GEOS stress not collected here)
    ax[2].semilogx(ad[:, 0], ad[:, 1], 'k-', lw=2, label='digitized Fig5d (paper)')
    ax[2].semilogx(sol['tau'], sol['sig'], 'g.-', ms=4, label='analytical')
    ax[2].set_title(r'STRESS $\sigma_{zz}$ (center)'); ax[2].set_ylabel(r'$\sigma/\sigma_0$')
    for a in ax:
        a.set_xlabel(r'$\tau = t/t_0$,  $t_0$=%.2f s' % t0)
        a.grid(alpha=.3); a.legend(); a.set_xlim(1e-5, tau_max)

    plt.tight_layout()
    out = os.path.join(HERE, "GEOS_vs_analytical_vs_digitized_full.png")
    plt.savefig(out, dpi=120)
    print("saved %s" % out)

if __name__ == '__main__':
    main()
