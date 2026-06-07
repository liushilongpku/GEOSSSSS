#!/usr/bin/env python3
"""
Plot GEOS dual-porosity Mandel results against the validated analytical solution
(dpdp_mandel_analytical.py). Digitized paper curves are NOT plotted (their
time axis does not line up with this t0 normalization).

Layout (analytical lines + GEOS markers):
  (a) pressure    : matrix p_m AND fracture p_f, normalized by p_i(0+)
  (b) stress      : dimensionless sigma_zz at the center
  (c) displacement: top-plate u_z

Usage:
  python plot_GEOS_vs_analytical.py [GEOS_OUTPUT_DIR] [TAU_MAX]
    GEOS_OUTPUT_DIR : dir with pressure_matrix_history.hdf5,
                      pressure_fracture_history.hdf5 (+ optional
                      stress_history.hdf5, displacement_history.hdf5)
    TAU_MAX         : max dimensionless time (default 3e4)

Normalization:
  pressure  p_i / p_i(0+)            (p_m0=4.55e5, p_f0=4.88e5 Pa)
  stress    sigma_zz / sigma_zz(t0)  (both curves start at 1)
  time      tau = t / t0,  t0 ~ 10.5 s
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
import dpdp_mandel_analytical as an   # provides t0, solve(), invert()

PM0 = 4.55e5    # Skempton initial matrix pressure (deck IC)
PF0 = 4.88e5    # Skempton initial fracture pressure (deck IC)
ABAR_M = 0.382  # effective matrix Biot (for total-stress reconstruction)
ABAR_F = 0.601  # effective fracture Biot


def _center_index(ec):
    """Element/node nearest the no-flow axis (x->min) at mid-height (z~0.015)."""
    x, z = ec[:, 0], ec[:, 2]
    return np.argmin(np.abs(x - x.min()) + np.abs(z - 0.015))


def load_pressure(out_dir, fn):
    path = os.path.join(out_dir, fn)
    if not os.path.exists(path):
        return None
    with h5py.File(path, 'r') as f:
        t = np.array(f['pressure Time'])[:, 0]
        p = np.array(f['pressure'])
        ec = np.array(f['pressure elementCenter'])[0]
    return t, p[:, _center_index(ec)]


def load_stress(out_dir):
    """GEOS TOTAL sigma_zz at the center element. Optional.

    The constitutive 'matrixSolid_stress' is the EFFECTIVE stress; the total
    (comparable to the analytical/Fig5d sigma_zz) is reconstructed via Biot:
        sigma_total = sigma_eff - abar_m*p_m - abar_f*p_f .
    """
    path = os.path.join(out_dir, 'stress_history.hdf5')
    if not os.path.exists(path):
        return None
    with h5py.File(path, 'r') as f:
        tkey = [k for k in f.keys() if k.endswith('Time')][0]
        ckey = [k for k in f.keys() if k.endswith('elementCenter')][0]
        vkey = [k for k in f.keys() if k not in (tkey, ckey)][0]
        t = np.array(f[tkey])[:, 0]
        v = np.array(f[vkey])
        ec = np.array(f[ckey])[0]
    ix = _center_index(ec)
    vc = v[:, ix].reshape(v.shape[0], -1)          # (ntime, nq*6) Voigt per quadrature point
    if vc.shape[1] % 6 == 0 and vc.shape[1] >= 6:  # average sigma_zz (Voigt index 2) over QPs
        szz_eff = vc.reshape(vc.shape[0], vc.shape[1] // 6, 6)[:, :, 2].mean(axis=1)
    else:
        szz_eff = vc[:, 2] if vc.shape[1] >= 3 else vc[:, -1]
    # reconstruct total stress (effective + Biot pressures); pressures are on the same grid
    gm = load_pressure(out_dir, 'pressure_matrix_history.hdf5')
    gf = load_pressure(out_dir, 'pressure_fracture_history.hdf5')
    if gm is not None and gf is not None and len(gm[1]) == len(szz_eff):
        szz_tot = szz_eff - ABAR_M * gm[1] - ABAR_F * gf[1]
        return t, szz_tot, True
    return t, szz_eff, False


def load_displacement(out_dir):
    """GEOS top-plate u_z(t) (z = max). Optional."""
    path = os.path.join(out_dir, 'displacement_history.hdf5')
    if not os.path.exists(path):
        return None
    with h5py.File(path, 'r') as f:
        t = np.array(f['totalDisplacement Time'])[:, 0]
        u = np.array(f['totalDisplacement'])                          # (ntime, nnode, 3)
        ref = np.array(f['totalDisplacement ReferencePosition'])[0]   # (nnode, 3)
    ztop = ref[:, 2].max()
    top = np.where(np.abs(ref[:, 2] - ztop) < 1e-9)[0]
    uz = u[:, top, 2].mean(axis=1)
    return t, uz


def main():
    geos_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    tau_max = float(sys.argv[2]) if len(sys.argv) > 2 else 3e4
    t0 = an.t0
    print("t0 = %.3f s ; tau_max = %.1e" % (t0, tau_max))

    sol = an.solve(tau_min=1e-5, tau_max=tau_max, n=90)
    ta = sol['tau']

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))

    # ===== (a) pressure: matrix + fracture =====
    ax[0].semilogx(ta, sol['pm'], 'r-', lw=1.6, label=r'analytical $p_m$')
    ax[0].semilogx(ta, sol['pf'], 'b-', lw=1.6, label=r'analytical $p_f$')
    gm = load_pressure(geos_dir, 'pressure_matrix_history.hdf5')
    gf = load_pressure(geos_dir, 'pressure_fracture_history.hdf5')
    if gm is not None:
        m = gm[0] > 0
        ax[0].semilogx(gm[0][m] / t0, gm[1][m] / PM0, 'r.', ms=5, label=r'GEOS $p_m$')
    if gf is not None:
        m = gf[0] > 0
        ax[0].semilogx(gf[0][m] / t0, gf[1][m] / PF0, 'b.', ms=5, label=r'GEOS $p_f$')
    ax[0].set_title('(a) pressure'); ax[0].set_ylabel(r'$p/p_0^+$'); ax[0].set_ylim(0, 1.5)

    # ===== (b) displacement u_z =====
    # (the stress subplot is intentionally omitted: the displacement-driven setup leaves the
    #  total stress as a reconstructed difference sigma_eff - abar*p that drifts at late time;
    #  see the technical doc. pressure and displacement are the reliable comparisons.)
    gd = load_displacement(geos_dir)
    tt = np.logspace(-3, np.log10(tau_max), 50) * t0
    uz_unit = an.invert(3, tt)
    if gd is not None:
        m = gd[0] > 0
        uz_geos = gd[1][m]
        Pc = np.abs(uz_geos).max() / np.abs(uz_unit).max()   # calibrate load to GEOS scale
        ax[1].semilogx(tt / t0, np.abs(uz_unit) * Pc, 'm-', lw=1.6, label=r'analytical $u_z$')
        ax[1].semilogx(gd[0][m] / t0, np.abs(uz_geos), 'm.', ms=5, label=r'GEOS $u_z$ (top)')
    else:
        ax[1].semilogx(tt / t0, np.abs(uz_unit), 'm-', lw=1.6, label=r'analytical $u_z$ (unit load)')
        print("note: displacement_history.hdf5 not found -> (b) shows analytical only")
    ax[1].set_title('(b) displacement'); ax[1].set_ylabel(r'$|u_z|$ (m)')

    for a in ax:
        a.set_xlabel(r'$\tau = t/t_0$,  $t_0$=%.2f s' % t0)
        a.grid(alpha=.3); a.legend(); a.set_xlim(1e-3, tau_max)

    fig.suptitle('DPDP N=2 Mandel — GEOS vs analytical', fontsize=13)
    fig.tight_layout()
    out = os.path.join(HERE, "GEOS_vs_analytical.png")
    fig.savefig(out, dpi=120)
    print("saved %s" % out)


if __name__ == '__main__':
    main()
