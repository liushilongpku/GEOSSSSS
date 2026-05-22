#!/usr/bin/env python3
"""
Collect GEOS output into a single HDF5 for verification.
All data from TimeHistory HDF5 files (no VTU parsing).
Sources: pressure_history.hdf5, displacement_history.hdf5, stress_history.hdf5
Output:  verification_data.hdf5
"""

import numpy as np
import h5py
import os
import sys

# ===========================================================================
# Physical constants (matching analytical solution)
# ===========================================================================
K_mat  = 1.1e9
nu     = 0.22
G_mod  = 3 * K_mat * (1 - 2*nu) / (2 * (1 + nu))
phi    = 0.14
Ks     = 27e9
Kf     = 1.744e9
k_perm = 4.9346e-21
mu_f   = 1e-3
a_len  = 3e-2
P_conf = 1e6
v1     = 0.97

alpha_BW = 1 - K_mat / Ks
inv_M    = (alpha_BW - phi) / Ks + phi / Kf
M_mod    = 1 / inv_M
K_u      = K_mat + alpha_BW**2 * M_mod
B_sk     = alpha_BW * M_mod / K_u
nu_u     = (3*K_u - 2*G_mod) / (2 * (3*K_u + G_mod))
kappa_1  = v1 * k_perm / mu_f
c_m1     = alpha_BW * (1 - 2*nu) / (2 * G_mod * (1 - nu))
S_dim    = inv_M + alpha_BW * c_m1
c_diff   = kappa_1 / S_dim
P0_ratio = B_sk * (1 + nu_u) / 3
P0_plus  = P0_ratio * P_conf
tau2t    = a_len**2 / c_diff


def main():
    geos_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
        '~/result/dpdkValidation/sp_new')
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        geos_dir, 'verification_data.hdf5')

    print(f'GEOS dir: {geos_dir}')
    print(f'Output:   {out_path}')

    # ---- 1. Pressure ----
    ph5 = os.path.join(geos_dir, 'pressure_history.hdf5')
    with h5py.File(ph5, 'r') as f:
        t_p = np.array(f['pressure Time']).squeeze()
        p_all = np.array(f['pressure'])
        ec_p = np.array(f['pressure elementCenter'])
    ei = np.argmin(np.abs(ec_p[0, :, 0]))
    p_centre = p_all[:, ei]
    tau_p = t_p / tau2t
    p_norm = p_centre / P0_plus
    print(f'  Pressure:  {len(t_p)} steps, tau=[{tau_p[0]:.1e}, {tau_p[-1]:.4f}]')

    # ---- 2. Displacement ----
    dh5 = os.path.join(geos_dir, 'displacement_history.hdf5')
    with h5py.File(dh5, 'r') as f:
        t_d = np.array(f['totalDisplacement Time']).squeeze()
        d_all = np.array(f['totalDisplacement zpos'])
        rp = np.array(f['totalDisplacement ReferencePosition zpos'])
    ni = np.argmin(np.abs(rp[0, :, 0] - a_len / 2))
    uz = d_all[:, ni, 2]
    tau_d = t_d / tau2t
    u_norm = 2 * G_mod * np.abs(uz) / (P_conf * a_len)
    print(f'  Displ:     {len(t_d)} steps, u_norm=[{u_norm[0]:.4f}, {u_norm[-1]:.4f}]')

    # ---- 3. Stress (from HDF5, not VTU) ----
    sh5 = os.path.join(geos_dir, 'stress_history.hdf5')
    with h5py.File(sh5, 'r') as f:
        t_s = np.array(f['matrixSolid_stress Time']).squeeze()
        stress_all = np.array(f['matrixSolid_stress'])       # (t, cells, 48)
        ec_s = np.array(f['matrixSolid_stress elementCenter'])
    ei_s = np.argmin(np.abs(ec_s[0, :, 0]))

    # stress layout: 48 = 6 components × 8 quadrature points
    # [qp0(XX,YY,ZZ,YZ,XZ,XY), qp1(...), ..., qp7(...)]
    # Extract ZZ component (index 2) at each qp, average across 8 qp
    n_steps = stress_all.shape[0]
    sig_eff_zz = np.zeros(n_steps)
    for t in range(n_steps):
        qp_vals = []
        for qp in range(8):
            idx = qp * 6 + 2  # ZZ component at qp
            qp_vals.append(stress_all[t, ei_s, idx])
        sig_eff_zz[t] = np.mean(qp_vals)

    # Total stress: σ_total = σ_eff - α·p
    p_centre_s = p_all[:, ei_s] if np.array_equal(ec_p, ec_s) else p_all[:, np.argmin(np.abs(ec_s[0, :, 0]))]
    sigma_total_zz =- sig_eff_zz + alpha_BW * p_centre_s

    tau_s = t_s / tau2t
    s_norm = np.abs(sigma_total_zz) / P_conf
    print(f'  Stress:    {len(tau_s)} steps, tau=[{tau_s[0]:.1e}, {tau_s[-1]:.4f}]')

    # ---- Write HDF5 ----
    with h5py.File(out_path, 'w') as f:
        f.create_dataset('tau_pressure', data=tau_p)
        f.create_dataset('p_norm', data=p_norm)
        f.create_dataset('tau_displacement', data=tau_d)
        f.create_dataset('u_norm', data=u_norm)
        f.create_dataset('tau_stress', data=tau_s)
        f.create_dataset('s_norm', data=s_norm)
    print(f'\nWritten: {out_path}')


if __name__ == '__main__':
    main()