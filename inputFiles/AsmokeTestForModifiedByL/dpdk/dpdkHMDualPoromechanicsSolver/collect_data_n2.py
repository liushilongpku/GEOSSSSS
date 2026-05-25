#!/usr/bin/env python3
"""Collect N2 GEOS output into verification_data.hdf5 (matching N1 format)."""
import h5py, numpy as np, os, sys, re, meshio

geos_dir = sys.argv[1] if len(sys.argv)>1 else os.path.expanduser('~/result/dpdkValidation/DualPoromechanicsSolverSinglePorevalidation')
out_path = sys.argv[2] if len(sys.argv)>2 else os.path.join(geos_dir, 'verification_data.hdf5')
tau2t = 102804.126
Pc = 1e6
alpha_m = 0.3819   # effective Biot coefficient for matrix
alpha_f = 0.6014   # effective Biot coefficient for fracture

# Pressure matrix — normalize to its own t=0 value
with h5py.File(f'{geos_dir}/pressure_matrix_history.hdf5','r') as f:
    t = np.array(f['pressure Time']).squeeze()
    p = np.array(f['pressure'])
    ec = np.array(f['pressure elementCenter'])
    ei = np.argmin(np.abs(ec[0,:,0]))
    tau_p = t / tau2t
    p_m0 = p[0, ei]
    p_norm = p[:, ei] / p_m0 if p_m0 > 0 else p[:, ei]

# Pressure fracture — normalize to its own t=0 value
with h5py.File(f'{geos_dir}/pressure_fracture_history.hdf5','r') as f:
    p_all = np.array(f['pressure'])
    pf = p_all[:, ei]
    pf0 = pf[0]
    pf_norm = pf / pf0 if pf0 > 0 else pf

# Displacement
with h5py.File(f'{geos_dir}/displacement_history.hdf5','r') as f:
    t_d = np.array(f['totalDisplacement Time']).squeeze()
    uz_all = np.array(f['totalDisplacement'])
    rp = np.array(f['totalDisplacement ReferencePosition'])
    zm = np.abs(rp[0,:,2]-0.03)<0.002
    xc = rp[0,:,0][zm]; ni = zm.nonzero()[0][np.argmin(np.abs(xc-0.015))]
    uz = np.abs(uz_all[:,ni,2])
    tau_d = t_d / tau2t
    G=3.11e8; b=0.03
    u_norm = 2*G*uz/(Pc*b)

# Stress from VTU (total sigma_zz = sigma_eff_zz - alpha_m*p_m - alpha_f*p_f)
vtk_dir = os.path.join(geos_dir, 'vtkOutput')
vtu_dirs = sorted([d for d in os.listdir(vtk_dir) if d.isdigit()], key=lambda x: int(x))
t_s_list, sz_list = [], []
for d in vtu_dirs:
    vtu = os.path.join(vtk_dir, d, 'mesh1', 'Level0', 'matrixRegion', 'rank_0.vtu')
    if not os.path.exists(vtu):
        vtu = os.path.join(vtk_dir, d, 'mesh1', 'Level0', 'Domain', 'rank_0.vtu')
    if not os.path.exists(vtu): continue
    with open(vtu, 'r') as fh: header = fh.read(5000)
    tm = re.search(r'TIME.*?RangeMin="([^"]+)"', header)
    t_val = float(tm.group(1)) if tm else float(d)
    m = meshio.read(vtu)
    stress = np.array(m.cell_data['matrixSolid_stress'])[0]
    pressure = np.array(m.cell_data['pressure'])[0]
    centers = np.array(m.cell_data['elementCenter'])[0]
    ci = np.argmin(np.abs(centers[:, 0]) + np.abs(centers[:, 2]))
    sigma_eff = stress[ci, 2]
    p_cell = pressure[ci]
    sigma_total = sigma_eff - alpha_m * p_cell  # matrix Biot contribution
    # Add fracture Biot contribution: approximate with centre p_f at this time
    idx_t = np.argmin(np.abs(t - t_val))
    p_f_at_t = float(pf[idx_t]) if idx_t < len(pf) else 0.0
    sigma_total -= alpha_f * p_f_at_t
    t_s_list.append(t_val)
    sz_list.append(np.abs(sigma_total) / Pc)

order = np.argsort(t_s_list)
tau_s = np.array(t_s_list)[order] / tau2t
s_norm = np.array(sz_list)[order]

with h5py.File(out_path, 'w') as f:
    f.create_dataset('tau_pressure', data=tau_p)
    f.create_dataset('p_norm', data=p_norm)
    f.create_dataset('p_m0', data=p_m0)           # matrix p(t=0) physical
    f.create_dataset('tau_pf', data=tau_p)
    f.create_dataset('pf_norm', data=pf_norm)
    f.create_dataset('pf0', data=pf0)              # fracture p(t=0) physical
    f.create_dataset('tau_displacement', data=tau_d)
    f.create_dataset('u_norm', data=u_norm)
    f.create_dataset('tau_stress', data=tau_s)
    f.create_dataset('s_norm', data=s_norm)

print(f'Written: {out_path}  ({len(tau_p)} pressure, {len(tau_d)} disp, {len(tau_s)} stress steps)')
print(f'  p_m(0)={p_m0:.4e} Pa, p_f(0)={pf0:.4e} Pa')
