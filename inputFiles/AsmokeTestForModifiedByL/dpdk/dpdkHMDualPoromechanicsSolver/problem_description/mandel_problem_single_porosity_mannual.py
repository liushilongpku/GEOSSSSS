"""Fig5_reproduce_final.py
Reproduce Figure 5(a)-(c) of Mehrabian & Abousleiman (2014)
Single-porosity Mandel problem — Python translation of MATLAB code.

Track B formulas (no gamma in pressure denominator):
  Pressure:  p_tilde/Pc = alpha*(1-sech(sqrt(s))) / [s * D(s)]
  Stress:    sigma_tilde/Pc = (psi - alpha*beta1*sech(sqrt(s))) / [s * D(s)]
  Displ:     2G*u_tilde/(Pc*b) = 1 / [s * D(s)]
  where D(s) = psi - alpha*beta1*tanh(sqrt(s))/sqrt(s)
  and   psi = alpha / (P0+/Pc)
"""

import numpy as np
from math import factorial, log, floor
import csv
import os
import time

# =============================================================================
# Physical parameters (Table 1, primary matrix)
# =============================================================================
K_mat  = 1.1e9       # Drained bulk modulus (Pa)
nu     = 0.22        # Drained Poisson's ratio
phi    = 0.14        # Porosity
Ks     = 27e9        # Solid grain bulk modulus (Pa)
Kf     = 1744e6      # Fluid bulk modulus (Pa)
k_d    = 5e-9        # Permeability (darcy)
mu_f   = 1e-3        # Fluid viscosity (Pa.s)
a_len  = 3e-2        # Half-length (m)
P_conf = 1.0e6       # Confining stress (Pa)
v1     = 0.97        # Matrix volume fraction

# =============================================================================
# Derived poroelastic constants
# =============================================================================
G         = 3 * K_mat * (1 - 2*nu) / (2 * (1 + nu))
alpha_BW  = 1 - K_mat / Ks
inv_M     = (alpha_BW - phi) / Ks + phi / Kf
M_mod     = 1 / inv_M
K_u       = K_mat + alpha_BW**2 * M_mod
B_sk      = alpha_BW * M_mod / K_u
nu_u      = (3*K_u - 2*G) / (2 * (3*K_u + G))

# Permeability & diffusivity
m2_per_darcy = 9.869233e-13
perm_m2  = k_d * m2_per_darcy
kappa_1  = v1 * perm_m2 / mu_f

# Poroelastic moduli
c_m1   = alpha_BW * (1 - 2*nu) / (2 * G * (1 - nu))
S_dim  = inv_M + alpha_BW * c_m1
c_diff = kappa_1 / S_dim

# Initial undrained pore pressure (Eq.46)
P0_ratio = B_sk * (1 + nu_u) / 3

print("========== Single-Porosity Mandel (Python) ==========")
print(f"  alpha   = {alpha_BW:.6f}")
print(f"  B       = {B_sk:.6f}")
print(f"  nu_u    = {nu_u:.6f}")
print(f"  P0+/Pc  = {P0_ratio:.6f}")
print(f"  G       = {G:.4e} Pa")
print(f"  c_diff  = {c_diff:.4e} m^2/s")

# =============================================================================
# N=1 reduction parameters
# =============================================================================
beta_1    = c_m1 / S_dim
gamma_val = 1 / (2 * G * S_dim)

# psi from initial condition (Track B: no gamma in pressure denominator)
psi_val = alpha_BW / P0_ratio

print(f"  beta1   = {beta_1:.6f}")
print(f"  gamma   = {gamma_val:.6f} (not used in pressure denominator)")
print(f"  psi     = {psi_val:.6f} (= alpha / (P0+/Pc))")
print(f"  1 + alpha*beta1 = {1 + alpha_BW*beta_1:.6f} (comparison)")

# =============================================================================
# Stehfest coefficients
# =============================================================================
def stehfest_weights(N):
    """Compute Stehfest inversion coefficients of order N (must be even)."""
    V = np.zeros(N)
    for k in range(1, N+1):
        sm = 0.0
        j_lo = (k + 1) // 2
        j_hi = min(k, N // 2)
        for j in range(j_lo, j_hi + 1):
            num = (j ** (N//2)) * factorial(2*j)
            den = (factorial(N//2 - j) * factorial(j) * factorial(j-1)
                   * factorial(k - j) * factorial(2*j - k))
            sm += num / den
        V[k-1] = ((-1) ** (k + N//2)) * sm
    return V

def stehfest_inv(Fs, t, N, V):
    """Numerical inverse Laplace transform via Stehfest algorithm."""
    if t <= 0:
        return 0.0
    ln2t = log(2) / t
    acc = 0.0
    for k in range(1, N+1):
        s_val = k * ln2t
        acc += V[k-1] * Fs(s_val)
    return ln2t * acc

# =============================================================================
# Laplace-domain functions (Track B: no gamma)
# =============================================================================
def safe_sq_C0_T(s):
    """Return (sq, C0, T) handling large-s limits analytically.
    C0 = 1/cosh(sq) -> 0   for sq > 50
    T  = tanh(sq)/sq -> 1/sq for sq > 30
    """
    sq = np.sqrt(max(abs(s), 1e-12))
    if sq > 50.0:
        C0 = 0.0
        T  = 1.0 / sq
    elif sq > 30.0:
        C0 = 1.0 / np.cosh(sq)  # safe, cosh(50) ~ 2.6e21
        T  = 1.0 / sq           # tanh(30)=1-2e-27, essentially 1
    else:
        C0 = 1.0 / np.cosh(sq)
        T  = np.tanh(sq) / sq
    return sq, C0, T

def phat_N1(s, alpha, beta1, psi):
    """Eq.(43) N=1, r=0 -- WITHOUT gamma factor (Track B)."""
    _, C0, T = safe_sq_C0_T(s)
    D = psi - alpha * beta1 * T
    return alpha * (1.0 - C0) / (s * D)

def sigmahat_N1(s, alpha, beta1, psi):
    """Eq.(44) N=1, r=0 (compression positive)."""
    _, C0, T = safe_sq_C0_T(s)
    D = psi - alpha * beta1 * T
    return (psi - alpha * beta1 * C0) / (s * D)

def uhat_N1(s, alpha, beta1, psi):
    """Eq.(45) N=1."""
    _, _, T = safe_sq_C0_T(s)
    D = psi - alpha * beta1 * T
    return 1.0 / (s * D)

# =============================================================================
# Time discretization
# =============================================================================
tau_vec = np.logspace(-6, 3, 300)
N_stehfest = 14
V_weights = stehfest_weights(N_stehfest)

# =============================================================================
# Compute solutions
# =============================================================================
print("Computing Stehfest inversions...", end=" ", flush=True)
t_start = time.perf_counter()

p_center = np.zeros_like(tau_vec)
sigma_c  = np.zeros_like(tau_vec)
uz_norm  = np.zeros_like(tau_vec)

# Pre-bind parameters for speed
a, b1, pv = alpha_BW, beta_1, psi_val

for i, tau in enumerate(tau_vec):
    if tau <= 0:
        continue
    p_center[i] = P_conf * stehfest_inv(
        lambda s: phat_N1(s, a, b1, pv), tau, N_stehfest, V_weights)
    sigma_c[i] = P_conf * stehfest_inv(
        lambda s: sigmahat_N1(s, a, b1, pv), tau, N_stehfest, V_weights)
    uz_norm[i] = stehfest_inv(
        lambda s: uhat_N1(s, a, b1, pv), tau, N_stehfest, V_weights)

elapsed = time.perf_counter() - t_start
print(f"done ({elapsed:.1f} s)")

# Normalize
p_norm = p_center / (P0_ratio * P_conf)   # p / P0+
s_norm = sigma_c / P_conf                  # sigma / Pc

# =============================================================================
# Diagnostics
# =============================================================================
pk_p = np.max(p_norm)
ipk  = np.argmax(p_norm)
pk_s = np.max(s_norm)
isk  = np.argmax(s_norm)

print("\n========== Results ==========")
print(f"  Pressure:   init={p_norm[0]:.4f}  peak={pk_p:.4f} @ tau={tau_vec[ipk]:.4f}  final={p_norm[-1]:.4f}")
print(f"  Stress:     init={s_norm[0]:.4f}  peak={pk_s:.4f} @ tau={tau_vec[isk]:.4f}  final={s_norm[-1]:.4f}")
print(f"  Displacement: init={uz_norm[0]:.4f}  final={uz_norm[-1]:.4f}")

# =============================================================================
# Load extracted data
# =============================================================================
BASE = os.path.dirname(os.path.abspath(__file__))

def load_csv(filename):
    path = os.path.join(BASE, filename)
    if not os.path.exists(path):
        return None
    x, y = [], []
    with open(path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
    return np.array(x), np.array(y)

fig5a = load_csv('fig5a.csv')
fig5b = load_csv('fig5b.csv')
has_data = fig5a is not None and fig5b is not None

if has_data:
    # Interpolate computed to extracted x-points for RMSE
    p_interp = np.interp(fig5a[0], tau_vec, p_norm)
    s_interp = np.interp(fig5b[0], tau_vec, s_norm)
    p_rmse = np.sqrt(np.nanmean((p_interp - fig5a[1])**2))
    s_rmse = np.sqrt(np.nanmean((s_interp - fig5b[1])**2))
    print(f"\n  RMSE vs extracted:")
    print(f"    Pressure  = {p_rmse:.4f}  (extracted peak={np.max(fig5a[1]):.4f})")
    print(f"    Stress    = {s_rmse:.4f}  (extracted peak={np.max(fig5b[1]):.4f})")

# =============================================================================
# Plot Figure 5(a)-(c)
# =============================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), facecolor='w')
fig.suptitle(
    "Single-Porosity Mandel Response  "
    "(cf. Mehrabian & Abousleiman 2014, Fig. 5a-c)",
    fontsize=12, fontweight='bold')

# ---- (a) Center pore pressure ----
ax = axes[0]
ax.semilogx(tau_vec, p_norm, 'k-', linewidth=1.7)
if has_data:
    ax.semilogx(fig5a[0], fig5a[1], 'ro', markersize=4, markerfacecolor='r')
ax.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
ax.grid(True)
ax.set_xlabel(r'$c t / a^2$')
ax.set_ylabel(r'$p / P_0^+$')
ax.set_title('(a) Center pore pressure')
if has_data:
    ax.legend(['Computed', 'Extracted (Fig.5a)'], loc='upper right', fontsize=8)
ax.set_xlim(1e-5, 1e2)
ax.set_ylim(-0.05, 1.2)
ax.set_facecolor('w')
ax.tick_params(colors='k')

# ---- (b) Center confining stress ----
ax = axes[1]
ax.semilogx(tau_vec, s_norm, 'k-', linewidth=1.7)
if has_data:
    ax.semilogx(fig5b[0], fig5b[1], 'ro', markersize=4, markerfacecolor='r')
ax.axhline(y=1, color='k', linestyle='--', linewidth=0.5)
ax.grid(True)
ax.set_xlabel(r'$c t / a^2$')
ax.set_ylabel(r'$\sigma / P_c$')
ax.set_title('(b) Center confining stress')
if has_data:
    ax.legend(['Computed', 'Extracted (Fig.5b)'], loc='lower right', fontsize=8)
ax.set_xlim(1e-5, 1e3)
ax.set_ylim(0.95, 1.20)
ax.set_facecolor('w')
ax.tick_params(colors='k')

# ---- (c) Normalized top displacement ----
ax = axes[2]
ax.semilogx(tau_vec, uz_norm, 'k-', linewidth=1.7)
ax.grid(True)
ax.set_xlabel(r'$c t / a^2$')
ax.set_ylabel(r'$2 G u_z / (P_c b)$')
ax.set_title('(c) Normalized top displacement')
ax.set_xlim(1e-5, 1e3)
ax.set_facecolor('w')
ax.tick_params(colors='k')

plt.tight_layout()
outpath = os.path.join(BASE, 'Fig5_reproduced.png')
fig.savefig(outpath, dpi=150, bbox_inches='tight', facecolor='w', edgecolor='w')
plt.close(fig)
print(f"\nFigure saved: {outpath}")

# =============================================================================
# Export data for GEOS verification
# =============================================================================

# Physical time: t = tau * a^2 / c_diff
t_phys = tau_vec * a_len**2 / c_diff

# Physical displacement: u_z = uz_norm * Pc * b / (2*G)  (b = a_len for square)
u_z_phys = uz_norm * P_conf * a_len / (2 * G)

# --- Export 1: Displacement BC for GEOS XML (time [s], u_z [m])
bc_path = os.path.join(BASE, 'displacement_bc.csv')
np.savetxt(bc_path, np.column_stack([t_phys, u_z_phys]),
           delimiter=',', header='time,u_z', comments='')
print(f"Displacement BC exported: {bc_path}  ({len(t_phys)} rows)")

# --- Export 2: Dimensionless reference for comparison
ref_path = os.path.join(BASE, 'dimensionless_ref.csv')
np.savetxt(ref_path, np.column_stack([tau_vec, p_norm, s_norm, uz_norm]),
           delimiter=',', header='tau,p_norm,sigma_norm,uz_norm', comments='')
print(f"Dimensionless reference exported: {ref_path}  ({len(tau_vec)} rows)")
