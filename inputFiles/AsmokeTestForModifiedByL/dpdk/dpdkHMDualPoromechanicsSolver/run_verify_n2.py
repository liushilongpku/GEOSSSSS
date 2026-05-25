#!/usr/bin/env python3
"""N2 verification: read verification_data.hdf5 -> SP_Mandel_N2_verify.png"""
import h5py, numpy as np, os, argparse, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Pc=1e6; b=0.03

p = argparse.ArgumentParser()
p.add_argument('--h5', type=str, default=os.path.expanduser('~/result/dpdkValidation/DualPoromechanicsSolverSinglePorevalidation/verification_data.hdf5'))
p.add_argument('--out', type=str, default=os.path.expanduser('~/result/dpdkValidation/SP_Mandel_N2_verify.png'))
args = p.parse_args()

with h5py.File(args.h5, 'r') as f:
    tau_m  = np.array(f['tau_pressure']); p_norm = np.array(f['p_norm'])
    tau_f  = np.array(f['tau_pf']);       pf_norm = np.array(f['pf_norm'])
    tau_d  = np.array(f['tau_displacement']); u_norm = np.array(f['u_norm'])
    tau_s  = np.array(f['tau_stress']);   s_norm = np.array(f['s_norm'])
    p_m0  = np.array(f['p_m0']); pf0 = np.array(f['pf0'])
    G = 6.789e8
    uz = u_norm * Pc * b / (2*G)

# Reference uz from CSV
csv = '/home/lsl/codes/GEOSSSSS/inputFiles/AsmokeTestForModifiedByL/dpdk/dpdkHMDualPoromechanicsSolver/u_z_physical.csv'
tau_ref, uz_ref = np.loadtxt(csv, delimiter=',', skiprows=1, unpack=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), facecolor='w')
fig.suptitle('N=2 Dual-Porosity Mandel Verification\n(Mehrabian & Abousleiman 2014)', fontsize=13, fontweight='bold')

def style(ax):
    ax.set_facecolor('w'); ax.grid(True, alpha=0.3)
    for s in ax.spines.values(): s.set_color('k')

# (a) p/p0 vs tau
ax = axes[0]
ax.semilogx(tau_m, p_norm, 'b-o', lw=1.5, label=r'$p_m/p_m^0$ (GEOS)')
ax.semilogx(tau_f, pf_norm, 'r--o', lw=1.2, label=r'$p_f/p_f^0$ (GEOS)')
ax.axhline(y=1.0, color='gray', ls=':', lw=0.7)
ax.set_xlabel(r'$\tau$'); ax.set_ylabel(r'$p/p^0$')
ax.set_title('(a) Pore pressure'); ax.legend(fontsize=8)
ax.set_xlim(1e-5, 1); style(ax)

# (b) sigma/Pc vs tau
ax = axes[1]
has_s = len(tau_s) > 1
if has_s:
    ax.semilogx(tau_s, s_norm, 'k-o', lw=1.2, ms=3, label=r'GEOS')
    ax.axhline(y=1.0, color='gray', ls=':', lw=0.7)
    ax.legend(fontsize=8)
ax.set_xlabel(r'$\tau$'); ax.set_ylabel(r'$\sigma / P_c$')
ax.set_title('(b) Centre confining stress')
ax.set_xlim(1e-5, 1); style(ax)

# (c) Diagnostics
ax = axes[2]; ax.axis('off')
s0_str = f'{s_norm[0]:.4f}' if has_s else 'N/A'
se_str = f'{s_norm[-1]:.4f}' if has_s else 'N/A'
lines = [
    'N=2 Parameters',
    'abar1 = 0.382, abar2 = 0.601',
    'Pc = 1 MPa',
    'pm(0+) = 0.455 Pc, pf(0+) = 0.488 Pc',
    'km = 4.93e-21, kf = 4.93e-18',
    '',
    'GEOS Results',
    '',
    '--- p/p0 ---',
    f'p_m0 = {p_m0/1e5:.4f} x10^5 Pa',
    f'p_f0  = {pf0/1e5:.4f} x10^5 Pa',
    f'p_m/p_m0(peak) = {p_norm.max():.4f}',
    '',
    '--- sigma/Pc ---',
    f'sigma/Pc(t=0) = {s0_str}',
    f'sigma/Pc(end) = {se_str}',
    f'steps = {len(tau_m)}',
]
for i, l in enumerate(lines):
    ax.text(0.05, 0.95-i*0.035, l, fontsize=9, fontfamily='monospace', va='top', transform=ax.transAxes)

plt.tight_layout()
fig.savefig(args.out, dpi=150, bbox_inches='tight', facecolor='w'); plt.close()
print(f'Saved: {args.out}')
