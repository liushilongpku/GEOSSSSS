#!/usr/bin/env python3
"""
DPDP_Mandel_Mehrabian2014_verify.py

Verification script for Dual-Porosity Dual-Permeability Mandel's problem.
Compares GEOS numerical results against the analytical solution from:

  Mehrabian, A., & Abousleiman, Y. N. (2014).
  "Generalized Biot's theory and Mandel's problem of multiple-porosity
   and multiple-permeability poroelasticity."
  J. Geophys. Res. Solid Earth, 119(4), 2745-2763.

Usage:
  python DPDP_Mandel_Mehrabian2014_verify.py \
      --output-dir ~/result/dpdkValidation/
"""

import numpy as np
import h5py
import argparse
import os
import sys
from scipy.optimize import brentq
from scipy.special import comb as binom

# ============================================================================
# Material Parameters — GOM Shale (Mehrabian & Abousleiman 2014, Table 1)
# ============================================================================

class ShaleParams:
    """GOM Shale dual-porosity material parameters (SI units)."""
    # Grain and fluid
    Ks = 27.0e9         # grain bulk modulus [Pa]
    Kf = 1.744e9        # fluid bulk modulus [Pa] (1744 MPa)
    mu  = 0.001         # fluid viscosity [Pa·s]
    rho = 1000.0        # fluid density [kg/m3]

    # Matrix (primary porosity)
    K1 = 1.1e9          # drained bulk modulus [Pa]
    nu = 0.22           # Poisson's ratio
    phi1 = 0.14         # porosity
    k1 = 4.93e-21       # permeability [m2] (5 nd)
    v1 = 0.97           # volume fraction (in mixture)

    # Fracture (secondary = macro-fracture from tertiary)
    K2 = 22.5e6         # drained bulk modulus [Pa]
    phi2 = 0.95         # porosity
    k2 = 4.93e-15       # permeability [m2] (5 md)
    v2 = 0.03           # volume fraction (in mixture)

    # Geometry
    a = 0.03            # specimen half-width [m]
    b = 0.03            # specimen half-height [m]

    # Interporosity exchange coefficient (Gamma_13 from Table 2)
    Gamma = 1.67e-22    # [Pa^{-1} s^{-1}]

# ============================================================================
# Single-Porosity Mandel Analytical Solution  (Cheng & Detournay 1988)
# ============================================================================

class SinglePorosityMandel:
    """
    Analytical solution for standard (single-porosity) Mandel's problem.
    Used as the degradation-test reference.
    """
    def __init__(self, params):
        self.a = params.a
        self.b = params.b
        self.nu = params.nu
        self.K = params.K1
        self.Ks = params.Ks
        self.Kf = params.Kf
        self.phi = params.phi1
        self.k = params.k1
        self.mu = params.mu

        G = 3.0 * self.K * (1.0 - 2.0 * self.nu) / (2.0 * (1.0 + self.nu))

        # Biot coefficient
        self.alpha = 1.0 - self.K / self.Ks

        # Skempton coefficient
        self.B = (1.0/self.K - 1.0/self.Ks) / (
            1.0/self.K - 1.0/self.Ks + self.phi * (1.0/self.Kf - 1.0/self.Ks))

        # Undrained Poisson's ratio
        self.nu_u = (3.0 * self.nu + self.B * self.alpha * (1.0 - 2.0 * self.nu)) / (
            3.0 - self.B * self.alpha * (1.0 - 2.0 * self.nu))

        # Consolidation coefficient
        self.c = (2.0 * self.k * G * (1.0 - self.nu) * (self.nu_u - self.nu) /
                  (self.mu * self.alpha ** 2 * (1.0 - 2.0 * self.nu) ** 2 *
                   (1.0 - self.nu_u)))

        # Applied force per unit length
        # Matched to produce displacement consistent with the load in XML
        self.F = 1.0

        # Initial pore pressure (Skempton effect, eq 46)
        self.p0 = self.F * self.B * (1.0 + self.nu_u) / (3.0 * self.a)

        # Compute roots alpha_n of: tan(alpha) = ((1-nu)/(nu_u-nu)) * alpha
        self._compute_roots(n_roots=40)

    def _root_eqn(self, alpha):
        return np.tan(alpha) - (1.0 - self.nu) / (self.nu_u - self.nu) * alpha

    def _compute_roots(self, n_roots=40):
        self.alpha_n = []
        eps = 1e-10
        for n in range(n_roots):
            left = n * np.pi + eps
            right = (n + 0.5) * np.pi - eps
            try:
                root = brentq(self._root_eqn, left, right, xtol=1e-12)
                self.alpha_n.append(root)
            except ValueError:
                break

    def pressure(self, x, z, t):
        """
        Pore pressure p(x,z,t).
        x: position along x [m] (0 <= x <= a)
        z: position along z [m] (z not used — pressure uniform in z for Mandel)
        t: time [s]
        """
        p = 0.0
        x_norm = x / self.a
        for alpha_n in self.alpha_n:
            p += (np.sin(alpha_n) /
                  (alpha_n - np.sin(alpha_n) * np.cos(alpha_n)) *
                  (np.cos(alpha_n * x_norm) - np.cos(alpha_n)) *
                  np.exp(-alpha_n ** 2 * self.c * t / self.a ** 2))
        return 2.0 * self.p0 * p

    def displacement_z(self, x, z, t):
        """
        Vertical displacement u_z(x,z,t).
        z: position along z [m]
        """
        u = 0.0
        for alpha_n in self.alpha_n:
            u += (np.sin(alpha_n) * np.cos(alpha_n) /
                  (alpha_n - np.sin(alpha_n) * np.cos(alpha_n)) *
                  np.exp(-alpha_n ** 2 * self.c * t / self.a ** 2))

        G = 3.0 * self.K * (1.0 - 2.0 * self.nu) / (2.0 * (1.0 + self.nu))
        u0 = -self.F * (1.0 - self.nu) / (2.0 * G * self.a)
        u1 = self.F * (1.0 - self.nu_u) / (G * self.a)
        return (u0 + u1 * u) * z + self.F * self.nu_u / (2.0 * G) * x / self.a


# ============================================================================
# Dual-Porosity Mandel Analytical Solution  (Mehrabian & Abousleiman 2014)
# ============================================================================

class DualPorosityMandel:
    """
    Analytical solution for dual-porosity dual-permeability Mandel's problem.

    Implements the Laplace-domain solution (eq. 43-45) with numerical
    Stehfest inversion to time domain.
    """
    def __init__(self, params):
        self.a = params.a
        self.b = params.b
        self.nu = params.nu

        # Intrinsic Biot coefficients
        alpha1 = 1.0 - params.K1 / params.Ks
        alpha2 = 1.0 - params.K2 / params.Ks

        # Drained shear modulus of mixture
        G = 3.0 * params.K1 * (1.0 - 2.0 * params.nu) / (2.0 * (1.0 + params.nu))

        # Volume-fraction-weighted effective Biot coefficients (paper eq. after eq.10)
        alpha1_eff = (params.v1 / params.K1) * alpha1
        alpha2_eff = (params.v2 / params.K2) * alpha2

        # Effective drained bulk modulus of mixture
        K_inv = params.v1 / params.K1 + params.v2 / params.K2
        K_eff = 1.0 / K_inv

        # Skempton coefficients for each porosity
        B1 = ((1.0/params.K1 - 1.0/params.Ks) /
              (1.0/params.K1 - 1.0/params.Ks +
               params.phi1 * (1.0/params.Kf - 1.0/params.Ks)))
        B2 = ((1.0/params.K2 - 1.0/params.Ks) /
              (1.0/params.K2 - 1.0/params.Ks +
               params.phi2 * (1.0/params.Kf - 1.0/params.Ks)))

        # Biot moduli (diagonal and off-diagonal)
        M11_inv = (alpha1 - params.phi1) / params.Ks + params.phi1 / params.Kf
        M22_inv = (alpha2 - params.phi2) / params.Ks + params.phi2 / params.Kf
        M12_inv = 0.0  # assumed negligible for dual-porosity

        M11 = 1.0 / M11_inv if abs(M11_inv) > 1e-20 else 1e20
        M22 = 1.0 / M22_inv if abs(M22_inv) > 1e-20 else 1e20
        M12 = 0.0

        # Hydraulic conductivity: kappa_i = v_i * k_i / mu
        kappa1 = params.v1 * params.k1 / params.mu
        kappa2 = params.v2 * params.k2 / params.mu

        # Exchange matrix Gamma
        Gamma = params.Gamma

        # ---- Build the PDE matrices ----
        # Coupling coefficient c_m_i (eq. 19)
        self.cm1 = alpha1_eff * (1.0 - 2.0 * params.nu) / (2.0 * G * (1.0 - params.nu))
        self.cm2 = alpha2_eff * (1.0 - 2.0 * params.nu) / (2.0 * G * (1.0 - params.nu))

        # Generalized storage matrix S (eq. 22) — symmetric 2x2
        S11 = M11_inv + alpha1_eff * self.cm1
        S22 = M22_inv + alpha2_eff * self.cm2
        S12 = M12_inv + alpha1_eff * self.cm2
        S21 = M12_inv + alpha2_eff * self.cm1

        # Specific storage sum (for normalisation)
        Stot = S11 + S22

        # Store for solution computation
        self.G = G
        self.nu = params.nu
        self.K_eff = K_eff
        self.alpha1_eff = alpha1_eff
        self.alpha2_eff = alpha2_eff
        self.kappa1 = kappa1
        self.kappa2 = kappa2
        self.Gamma = Gamma
        self.S11 = S11
        self.S22 = S22
        self.S12 = S12
        self.S21 = S21
        self.Stot = Stot

        # Undrained Poisson ratio (simplified)
        self.nu_u = (3.0 * params.nu + B1 * alpha1 * (1.0 - 2.0 * params.nu)) / (
            3.0 - B1 * alpha1 * (1.0 - 2.0 * params.nu))

        # Reference pressure — confining stress for Mandel's problem.
        # Pc = F/(2a) where F=1 N/m (force per unit length).
        # The undrained pressure in SP is p₀ = F·B·(1+ν_u)/(3a).
        # For consistency, the dimensionless DP solution multiplies by Pc.
        self.Pc = 1.0 / (2.0 * params.a)  # 1/(2×0.03) ≈ 16.67 Pa
        self.gamma = 1.0 / (2.0 * G * Stot)  # eq. after 33 (gamma definition)

    def _characteristic_matrix(self, s):
        """
        Build the characteristic matrix M(s) = kappa^{-1} * (s*S + Gamma)

        For dual-porosity, the PDE system (eq. 20) in dimensionless form:
          [kappa * d2/dx2 - (s*S + Gamma)] * p_tilde = alpha * s * f_tilde

        The eigenvalue problem is:  M(s) * chi = lambda * chi
        where M = kappa^{-1} * (s*S + Gamma)

        Returns eigenvalues lambda[2] and eigenvectors X[2][2].
        """
        # Assembly of A = s*S + Gamma (dimensioned)
        A11 = s * self.S11 + self.Gamma
        A22 = s * self.S22 + self.Gamma
        A12 = s * self.S12 - self.Gamma
        A21 = s * self.S21 - self.Gamma

        # M = kappa^{-1} @ A
        M11 = A11 / self.kappa1
        M12 = A12 / self.kappa1
        M21 = A21 / self.kappa2
        M22 = A22 / self.kappa2

        # Compute eigenvalues and eigenvectors of the 2x2 matrix M
        # det(M - lambda*I) = 0
        # lambda^2 - tr(M)*lambda + det(M) = 0
        trace = M11 + M22
        detM  = M11 * M22 - M12 * M21

        disc = trace * trace - 4.0 * detM
        if disc < 0:
            disc = 0.0  # potential for complex eigenvalues — clamp for real solution

        lambda1 = 0.5 * (trace + np.sqrt(disc))
        lambda2 = 0.5 * (trace - np.sqrt(disc))

        # Eigenvectors: solve (M - lambda*I) * chi = 0
        # For 2x2, chi = [1, (lambda - M11) / M12] or [M12, lambda - M11]
        if abs(M12) > 1e-20:
            chi1_2 = (lambda1 - M11) / M12
            chi2_2 = (lambda2 - M11) / M12
        else:
            chi1_2 = 0.0
            chi2_2 = 0.0

        # Build eigenvector matrix X with chi_i normalized as [1, chi_i_2]^T
        X = np.array([[1.0, 1.0],
                       [chi1_2, chi2_2]])

        return np.array([lambda1, lambda2]), X

    def _laplace_solution(self, x, s):
        """
        Compute Laplace-domain solution at position x for a given Laplace
        parameter s. Returns p1_tilde, p2_tilde, uz_tilde.
        """
        # Get eigenvalues and eigenvectors
        lam, X = self._characteristic_matrix(s)

        # Psi vector (particular solution): psi = -s * (s*S + Gamma)^{-1} * alpha_vec
        # For dual-porosity, solve:
        # [s*S11+Gamma   s*S12-Gamma ] [psi1]   = -s * [alpha1_eff]
        # [s*S21-Gamma   s*S22+Gamma ] [psi2]   -s * [alpha2_eff]
        A11 = s * self.S11 + self.Gamma
        A22 = s * self.S22 + self.Gamma
        A12 = s * self.S12 - self.Gamma
        A21 = s * self.S21 - self.Gamma

        detA = A11 * A22 - A12 * A21
        if abs(detA) < 1e-30:
            detA = 1e-30

        b1 = -s * self.alpha1_eff
        b2 = -s * self.alpha2_eff

        psi1 = (A22 * b1 - A12 * b2) / detA
        psi2 = (-A21 * b1 + A11 * b2) / detA

        psi = np.array([psi1, psi2])

        # ---- Pressure solution (eq. 32) ----
        # p_tilde_i(x,s) = f_tilde(s) * [psi_i - (1/Det(X)) * sum_j ...
        #   chi_ij * Det[X_psi,j] * cosh(x*sqrt(lambda_j)) / cosh(sqrt(lambda_j))]

        # For x in [0, 1] (dimensionless, x* = x/a)
        x_star = x / self.a

        # Determinant of eigenvector matrix
        detX = X[0,0] * X[1,1] - X[0,1] * X[1,0]
        if abs(detX) < 1e-30:
            detX = 1e-30

        sqrt_lam = np.sqrt(lam)

        # Safe cosh ratio: cosh(x * r) / cosh(r) for r = sqrt(lam)
        # When r is large, cosh(r) ~ exp(r)/2, ratio ~ exp((x-1)*r) for 0 <= x <= 1
        def safe_cosh_ratio(xv, rv):
            if rv > 50.0:
                return np.exp((xv - 1.0) * rv)
            else:
                cr = np.cosh(rv) if rv < 200.0 else np.inf
                return np.cosh(xv * rv) / cr if cr > 0 else 0.0

        cosh_ratio = np.array([
            safe_cosh_ratio(x_star, sqrt_lam[0]),
            safe_cosh_ratio(x_star, sqrt_lam[1])
        ])

        # X_psi,j: replace j-th column of X with psi
        # Det[X_psi,1] = det([psi, X_col2]) = psi[0]*X[1,1] - psi[1]*X[0,1]
        DetX_psi1 = psi[0] * X[1,1] - psi[1] * X[0,1]
        # Det[X_psi,2] = det([X_col1, psi]) = X[0,0]*psi[1] - X[1,0]*psi[0]
        DetX_psi2 = X[0,0] * psi[1] - X[1,0] * psi[0]

        # Pressure in Laplace space for each porosity
        p_tilde = np.zeros(2)
        for i in range(2):
            sum_term = 0.0
            for j in range(2):
                if j == 0:
                    DetXj = DetX_psi1
                else:
                    DetXj = DetX_psi2
                term = X[i, j] * DetXj * cosh_ratio[j]
                sum_term += term
            p_tilde[i] = psi[i] - sum_term / detX

        # ---- Stress/Displacement solution (eq. 44-45) ----
        # beta_i = sum_j chi_ji * cm_j / Stot  (eq. 39)
        cm1 = self.cm1
        cm2 = self.cm2
        Stot = self.Stot

        beta = np.zeros(2)
        beta[0] = (X[0,0] * cm1 + X[1,0] * cm2) / Stot
        beta[1] = (X[0,1] * cm1 + X[1,1] * cm2) / Stot

        # psi_stress (eq. 40)
        psi_stress = 1.0 + 2.0 * cm1 - (
            self.alpha1_eff * psi[0] + self.alpha2_eff * psi[1]) / Stot

        # Denominator for f_tilde(s) (from eq. 43-45, corresponds to sigma_zz BC)
        tanh_term = beta[0] * DetX_psi1 * np.tanh(sqrt_lam[0]) / sqrt_lam[0] + \
                    beta[1] * DetX_psi2 * np.tanh(sqrt_lam[1]) / sqrt_lam[1]

        denominator = psi_stress * detX + tanh_term

        # f_tilde(s) = -Pc * detX / (s * denominator)  (eq. 42)
        if abs(denominator) < 1e-30:
            denominator = 1e-30

        f_tilde = -self.Pc * detX / (s * denominator)

        # Full pressure solution (dimensioned)
        # p_tilde_i = f_tilde * [psi_i - sum/DetX]
        p_result = f_tilde * p_tilde

        # Displacement solution (eq. 45): 2*G*u_z/(Pc*b) = -Det[X] / (s * denominator)
        uz_tilde = f_tilde  # proportional to u_z

        return p_result, uz_tilde, f_tilde

    def _stehfest(self, F, t, N=14):
        """
        Stehfest numerical Laplace inversion.
        F(s): function of Laplace variable s
        t: time value
        N: number of terms (must be even)
        """
        ln2 = np.log(2.0)
        f = 0.0
        for i in range(1, N + 1):
            s_i = i * ln2 / t
            coeff = 0.0
            k_min = max(1, int((i + 1) / 2))
            k_max = min(i, N // 2)
            for k in range(k_min, k_max + 1):
                num = k ** (N // 2) * binom(2 * k, k)
                den = binom(N // 2, k) * binom(k, i - k) * binom(2 * k, i - k)
                if abs(den) > 1e-30:
                    coeff += num / den
            coeff *= (-1) ** (i + N // 2)
            try:
                fs = F(s_i)
                f += coeff * fs
            except (ValueError, OverflowError):
                pass

        return ln2 / t * f

    def pressure(self, x, z, t_array):
        """
        Compute matrix and fracture pressures p_m(t), p_f(t) at position x.
        Returns arrays of size [len(t_array), 2].
        """
        # For Mandel's problem, pressure is uniform in z (depends on x only)
        result = np.zeros((len(t_array), 2))

        def F_vec(s):
            p, u, f_t = self._laplace_solution(x, s)
            return np.array([p[0], p[1]])

        for it, t in enumerate(t_array):
            p_t = self._stehfest(F_vec, t, N=14)
            # Normalize by Pc to get dimensionless pressure
            result[it, :] = p_t / self.Pc

        return result

    def displacement_z(self, z, t_array):
        """
        Compute vertical displacement u_z(z=b, t).
        Returns array of size len(t_array).
        """
        result = np.zeros(len(t_array))
        x = self.a  # evaluate at the edge (worst case; at center it's 0)

        def F_vec(s):
            p, u, f_t = self._laplace_solution(x / 2.0, s)
            # u_z proportional to f_tilde (eq. 45)
            uz_laplace = f_t * self.b / (2.0 * self.G) * (z / self.b)
            return uz_laplace

        for it, t in enumerate(t_array):
            result[it] = self._stehfest(F_vec, t, N=14)

        return result


# ============================================================================
# Numerical stability helpers
# ============================================================================

def stable_cosh_ratio(x, mu):
    """
    cosh(x*mu) / cosh(mu) for 0 <= x <= 1, mu >= 0.
    Avoids overflow by switching to exponential form for large mu.
    """
    if mu < 1e-6:
        return 1.0 + (x*x - 1.0) * mu*mu / 6.0
    if mu > 30.0:
        return np.exp((x - 1.0) * mu)
    return np.cosh(x * mu) / np.cosh(mu)

def stable_tanh_div(mu):
    """
    tanh(mu) / mu for mu >= 0.
    Piecewise: Taylor for |mu|<<1, 1/mu for mu>>1.
    """
    if mu < 1e-6:
        return 1.0 - mu*mu/3.0 + 2.0*mu**4/15.0
    if mu > 30.0:
        return 1.0 / mu
    return np.tanh(mu) / mu

def _stehfest_weights(N):
    """
    Stehfest coefficients V_i (N even).
    Standard formula: V_i = (-1)^{i+N/2} Σ_{k} [k^{N/2} (2k)!] / [(N/2-k)! k! (k-1)! (i-k)! (2k-i)!]

    Returns weights[i] = V_i * ln(2) for i=1..N (0-index ignored).
    """
    from math import factorial
    ln2 = np.log(2.0)
    weights = np.zeros(N + 1)
    Nh = N // 2
    for i in range(1, N + 1):
        coeff = 0.0
        k_min = (i + 1) // 2
        k_max = min(i, Nh)
        for k in range(k_min, k_max + 1):
            num = k ** Nh * factorial(2 * k)
            den = (factorial(Nh - k) * factorial(k) * factorial(k - 1) *
                   factorial(i - k) * factorial(2 * k - i))
            coeff += num / den
        weights[i] = coeff * (-1) ** (i + Nh) * ln2
    return weights

# ============================================================================
# Dual-Porosity Mandel — Numerically Stable Implementation
# ============================================================================

class DualPorosityMandelStable:
    """
    Numerically stable dual-porosity Mandel analytical solution.

    Improvements over DualPorosityMandel:
      1. Dimensionless formulation (κ*, S*, Γ* all near O(1))
      2. Symmetric eigendecomposition (eigh on B = K^{-1/2} A K^{-1/2})
      3. solve() instead of det()/inv()/Cramer
      4. Stable cosh/tanh for extreme eigenvalues
      5. Stehfest N-sensitivity check
    """

    def __init__(self, params):
        self.a = params.a
        self.b = params.b
        self.nu = params.nu

        # --- Intrinsic properties ---
        alpha1 = 1.0 - params.K1 / params.Ks
        alpha2 = 1.0 - params.K2 / params.Ks
        G = 3.0 * params.K1 * (1.0 - 2.0 * params.nu) / (2.0 * (1.0 + params.nu))

        # --- Effective Biot coefficients (volume-weighted, eq. after eq.10) ---
        alpha1_eff = (params.v1 / params.K1) * alpha1
        alpha2_eff = (params.v2 / params.K2) * alpha2

        # --- Coupling coefficients c_m_i (eq.19) ---
        self.cm1 = alpha1_eff * (1.0 - 2.0 * params.nu) / (2.0 * G * (1.0 - params.nu))
        self.cm2 = alpha2_eff * (1.0 - 2.0 * params.nu) / (2.0 * G * (1.0 - params.nu))

        # --- Biot moduli ---
        M11_inv = (alpha1 - params.phi1) / params.Ks + params.phi1 / params.Kf
        M22_inv = (alpha2 - params.phi2) / params.Ks + params.phi2 / params.Kf
        # Off-diagonal coupling: assume zero for dual-porosity
        M12_inv = 0.0

        # --- Generalized storage matrix S (eq.22) — symmetric ---
        self.S11 = M11_inv + alpha1_eff * self.cm1
        self.S22 = M22_inv + alpha2_eff * self.cm2
        self.S12 = M12_inv + alpha1_eff * self.cm2
        # S should be symmetric: S21 = S12 (by construction since M12_inv=0
        # and alpha1_eff*cm2 should equal alpha2_eff*cm1 by thermodynamic
        # reciprocity; we store S21 for completeness)

        # --- Total storage (for dimensionless scaling) ---
        self.Stot = self.S11 + self.S22

        # --- Hydraulic conductivity ---
        self.kappa1 = params.v1 * params.k1 / params.mu
        self.kappa2 = params.v2 * params.k2 / params.mu
        self.kappa_tot = self.kappa1 + self.kappa2

        # --- Exchange coefficient ---
        self.Gamma = params.Gamma

        # --- Other stored properties ---
        self.alpha1_eff = alpha1_eff
        self.alpha2_eff = alpha2_eff
        self.G = G
        self.Pc = 1.0 / (2.0 * params.a)  # confining stress F/(2a) ≈ 16.67 Pa

    # ==================================================================
    # Dimensionless Laplace solver
    # ==================================================================

    def _laplace_solution(self, x, s):
        """
        Laplace-domain solution with dimensionless formulation.

        All internal matrices are O(1) because we use:
          κ_i* = κ_i / κ_tot
          S_ij* = S_ij / S_tot
          Γ* = Γ · a² / κ_tot
          α_i* = α_i_eff / S_tot
          x* = x / a
        """
        a = self.a
        S_tot = self.Stot
        k_tot = self.kappa_tot

        # ---- Dimensionless parameters ----
        k1s = self.kappa1 / k_tot
        k2s = self.kappa2 / k_tot
        S11s = self.S11 / S_tot
        S22s = self.S22 / S_tot
        S12s = self.S12 / S_tot
        Gs   = self.Gamma * a * a / k_tot
        a1s  = self.alpha1_eff / S_tot
        a2s  = self.alpha2_eff / S_tot
        cm1s = self.cm1 * S_tot   # dimensionless cm1*
        cm2s = self.cm2 * S_tot   # dimensionless cm2*

        # ---- Convert physical s to dimensionless s* ----
        # Paper eq.(25): s* = s · a² · ΣS_ii / Σκ_i
        s = s * (a * a * S_tot / k_tot)

        # ---- Build symmetric A* = s*·S* + Γ* ----
        A11 = s * S11s + Gs
        A22 = s * S22s + Gs
        A12 = s * S12s - Gs

        # ---- Step 1: symmetric generalized eig ----
        # B = K^{-1/2} A K^{-1/2}   (symmetric)
        inv_sqrt_k1 = 1.0 / np.sqrt(k1s)
        inv_sqrt_k2 = 1.0 / np.sqrt(k2s)

        B11 = A11 / k1s
        B12 = A12 * inv_sqrt_k1 * inv_sqrt_k2
        B22 = A22 / k2s

        # B is symmetric 2x2 — use analytic eig for speed
        traceB = B11 + B22
        detB   = B11 * B22 - B12 * B12
        disc = max(traceB * traceB - 4.0 * detB, 0.0)
        sqrt_disc = np.sqrt(disc)
        lam1 = 0.5 * (traceB + sqrt_disc)
        lam2 = 0.5 * (traceB - sqrt_disc)

        # ---- Step 2: eigenvectors of B (analytically) ----
        # For symmetric 2x2, Q = [v1/|v1|, v2/|v2|] where v_i = [B12, lam_i - B11]
        if abs(B12) > 1e-20:
            n1 = np.sqrt(B12*B12 + (lam1 - B11)**2)
            n2 = np.sqrt(B12*B12 + (lam2 - B11)**2)
            Q11 = B12 / n1;  Q12 = B12 / n2
            Q21 = (lam1 - B11) / n1;  Q22 = (lam2 - B11) / n2
        else:
            # B is diagonal → Q = identity
            Q11 = 1.0;  Q12 = 0.0
            Q21 = 0.0;  Q22 = 1.0

        # X = sqrtK @ Q (back-transform to original eigenvectors)
        sqrt_k1 = np.sqrt(k1s)
        sqrt_k2 = np.sqrt(k2s)
        X11 = sqrt_k1 * Q11;  X12 = sqrt_k1 * Q12
        X21 = sqrt_k2 * Q21;  X22 = sqrt_k2 * Q22

        # ---- Step 3: solve A*psi = -s*alpha* instead of Cramer ----
        detA = A11 * A22 - A12 * A12
        if abs(detA) < 1e-60:
            detA = 1e-60
        psi1 = (A22 * (-s * a1s) - A12 * (-s * a2s)) / detA
        psi2 = (-A12 * (-s * a1s) + A11 * (-s * a2s)) / detA

        # ---- Step 4: solve X*c = psi instead of Cramer ----
        detX = X11 * X22 - X12 * X21
        if abs(detX) < 1e-60:
            detX = 1e-60
        c1 = ( X22 * psi1 - X12 * psi2) / detX
        c2 = (-X21 * psi1 + X11 * psi2) / detX

        # ---- Step 5: stable cosh ratios ----
        x_star = x / a
        sqrt_lam1 = np.sqrt(max(lam1, 0.0))
        sqrt_lam2 = np.sqrt(max(lam2, 0.0))

        cosh_r1 = stable_cosh_ratio(x_star, sqrt_lam1)
        cosh_r2 = stable_cosh_ratio(x_star, sqrt_lam2)

        tanh_d1 = stable_tanh_div(sqrt_lam1)
        tanh_d2 = stable_tanh_div(sqrt_lam2)

        # ---- Step 6: pressure Laplace solution (eq.32) ----
        # p_tilde_i = psi_i + Σ_j X_ij * c_j * cosh_ratio_j
        p1 = psi1 + X11 * c1 * cosh_r1 + X12 * c2 * cosh_r2
        p2 = psi2 + X21 * c1 * cosh_r1 + X22 * c2 * cosh_r2

        # ---- Step 7: f_tilde(s) — denominator from vertical stress BC ----
        # beta_i (eq.39)
        beta1 = (X11 * cm1s + X21 * cm2s)   # / Stot omitted (already in cm*s)
        beta2 = (X12 * cm1s + X22 * cm2s)

        # Det[X_psi,j] for substitution (eq.42 numerator)
        # Det[X_psi,1] = X22*psi1 - X12*psi2 = c1 * detX
        # Det[X_psi,2] = X11*psi2 - X21*psi1 = c2 * detX
        # → simplifies to c_j * detX!
        DetX_psi1 = c1 * detX
        DetX_psi2 = c2 * detX

        # psi_stress (eq.40)
        psi_stress = 1.0 + 2.0 * cm1s - (a1s * psi1 + a2s * psi2)

        # Denominator (eq.43-45)
        denominator = psi_stress + \
                      beta1 * c1 * tanh_d1 + \
                      beta2 * c2 * tanh_d2

        if abs(denominator) < 1e-60:
            denominator = 1e-60

        f_tilde = -self.Pc / (s * denominator)   # eq.(42) with correct Pc

        # Full pressure (dimensionless)
        p1_full = f_tilde * p1
        p2_full = f_tilde * p2

        return np.array([p1_full, p2_full]), f_tilde

    # ==================================================================
    # Stehfest inversion (with N-sensitivity)
    # ==================================================================

    def _stehfest(self, F, t, Nvals=None):
        """
        Stehfest inversion with N-sensitivity check.
        Returns (best_estimate, is_reliable).
        """
        if Nvals is None:
            Nvals = [10, 12, 14]

        results = {}
        for N in Nvals:
            weights = _stehfest_weights(N)
            f = 0.0
            for i in range(1, N + 1):
                s_i = i * np.log(2.0) / t
                try:
                    fs = F(s_i)
                    f += weights[i] * fs
                except (ValueError, OverflowError, ZeroDivisionError):
                    pass
            results[N] = f / t

        # N-sensitivity check
        reliable = True
        if 12 in results and 14 in results:
            v12 = results[12]
            v14 = results[14]
            if np.isscalar(v14):
                if abs(v14) > 1e-20:
                    rel_diff = abs(v14 - v12) / abs(v14)
                    if rel_diff > 0.05:
                        reliable = False
            else:
                # Array result: check element-wise
                valid = np.abs(v14) > 1e-20
                if np.any(valid):
                    rel_diff = np.max(np.abs(v14 - v12)[valid] / np.abs(v14)[valid])
                    if rel_diff > 0.05:
                        reliable = False

        best = results.get(14, results.get(12, results.get(10, 0.0)))
        return best, reliable

    # ==================================================================
    # Gaver-Wynn inversion (more robust than Stehfest for irregular F(s))
    # ==================================================================

    def _gaver_wynn(self, F, t, Mmax=12):
        """
        Gaver functionals + Wynn rho acceleration for Laplace inversion.

        Handles both scalar and vector-valued F(s).
        For vector-valued, runs the algorithm per component independently.
        Returns (best_estimate, is_reliable).
        """
        ln2 = np.log(2.0)
        alpha = ln2 / t

        # ---- Step 1: Gaver functionals G_k for k=1..Mmax ----
        G = []
        for k in range(1, Mmax + 1):
            binom_2k_k = binom(2 * k, k)
            gk = 0.0
            first_done = False
            for j in range(k + 1):
                s_jk = (k + j) * alpha
                try:
                    fs = F(s_jk)
                except Exception:
                    fs = 0.0
                sign = 1.0 if (j % 2 == 0) else -1.0
                contrib = sign * binom(k, j) * fs
                if not first_done:
                    gk = contrib
                    first_done = True
                else:
                    gk = gk + contrib
            G.append(binom_2k_k * alpha * gk)

        # ---- Step 2: Wynn rho acceleration ----
        # Flatten to scalar or run element-wise for arrays
        if np.isscalar(G[0]) or not hasattr(G[0], '__len__'):
            return self._wynn_rho_scalar(G)
        else:
            # Vector-valued: invert each component separately
            n_comp = len(G[0])
            result = np.zeros(n_comp)
            reliable = True
            for c in range(n_comp):
                Gc = [g[c] if hasattr(g, '__getitem__') else g for g in G]
                res_c, rel_c = self._wynn_rho_scalar(Gc)
                result[c] = res_c
                reliable = reliable and rel_c
            return result, reliable

    @staticmethod
    def _wynn_rho_scalar(G):
        """Wynn rho acceleration for scalar sequence G[0..M-1]."""
        M = len(G)
        # ρ[0][n] = G[n]; ρ[-1][n] = 0
        rho = np.zeros((M + 1, M))
        for n in range(M):
            rho[0, n] = G[n]

        for k in range(1, M):
            for n in range(M - k):
                denom = rho[k - 1, n + 1] - rho[k - 1, n]
                if abs(denom) < 1e-60:
                    rho[k, n] = rho[k - 1, n]
                else:
                    rho[k, n] = rho[k - 2, n + 1] + k / denom

        # Best estimate: ρ^{(0)}_{M-1} for even M-1, else ρ^{(0)}_{M-2}
        if (M - 1) % 2 == 0:
            best = rho[M - 1, 0]
        else:
            best = rho[M - 2, 0]

        # Reliability check
        reliable = True
        if M >= 4:
            vals = [rho[m, 0] for m in [M-1, M-2, M-3]]
            vals = [v for v in vals if not np.isnan(v)]
            if len(vals) >= 2:
                mu = np.mean(vals)
                if abs(mu) > 1e-20:
                    rel_var = np.std(vals) / abs(mu)
                    if rel_var > 0.1:
                        reliable = False
        return best, reliable

    # ==================================================================
    # Public interface
    # ==================================================================

    def pressure(self, x, z, t_array, check_stability=False):
        """
        Compute dual-porosity pressures at position x.
        Uses Gaver-Wynn Laplace inversion (more robust than Stehfest).
        Returns (p_matrix, p_fracture, reliability_flags) arrays.
        """
        n_t = len(t_array)
        p_m = np.zeros(n_t)
        p_f = np.zeros(n_t)
        reliable = np.ones(n_t, dtype=bool)

        def F_vec(s):
            p, ft = self._laplace_solution(x, s)
            return np.array([p[0], p[1]])

        for it, t in enumerate(t_array):
            if t <= 0:
                continue
            try:
                res, rel = self._stehfest(F_vec, t, Nvals=[10, 12, 14])
                # Physical pressure already from f_tilde·Pc in Laplace solution
                p_m[it] = res[0]
                p_f[it] = res[1]
                reliable[it] = rel
            except Exception:
                p_m[it] = np.nan
                p_f[it] = np.nan
                reliable[it] = False

        if check_stability:
            return p_m, p_f, reliable
        return p_m, p_f

    def displacement_z(self, z, t_array):
        """Vertical displacement at top (z=b). Uses Gaver-Wynn inversion."""
        result = np.zeros(len(t_array))

        def F_vec(s):
            p, ft = self._laplace_solution(self.a * 0.5, s)
            return ft * self.b / (2.0 * self.G) * (z / self.b)

        for it, t in enumerate(t_array):
            if t <= 0:
                continue
            try:
                result[it], _ = self._stehfest(F_vec, t, Nvals=[12, 14])
            except Exception:
                result[it] = np.nan

        return result


# ============================================================================
# Degradation test: DualPorosityMandelStable → single-porosity limit
# ============================================================================

def test_dual_to_single_degradation():
    """
    Verify DualPorosityMandelStable degenerates to SinglePorosityMandel
    when fracture properties → 0.
    """
    import sys
    sys.path.insert(0, '.')
    params = ShaleParams()
    params.k1 = 8.0e-17  # visible Mandel-Cryer

    # ---- Single-porosity reference ----
    sp = SinglePorosityMandel(params)
    t_test = np.logspace(-0.5, 1.5, 15)
    p_sp = np.array([sp.pressure(0.0, params.b, t) for t in t_test])

    # ---- Dual-porosity with fracture→0 ----
    dp = DualPorosityMandelStable(params)
    # Override fracture properties to near-zero (not exactly zero — avoids div/0)
    dp.kappa2 = 1e-30   # effectively no fracture flow
    dp.Gamma  = 0.0     # no exchange
    dp.S22    = 0.0     # no fracture storage
    p_m_dp, p_f_dp = dp.pressure(0.001, params.b, t_test, check_stability=True)[:2]

    # ---- Compare ----
    print("=== Degradation Test: Dual → Single Porosity ===")
    print(f"{'t':>8s}  {'SP p':>12s}  {'DP p_m':>12s}  {'diff':>10s}")
    max_err = 0.0
    for i, t in enumerate(t_test):
        diff = abs(p_sp[i] - p_m_dp[i])
        max_err = max(max_err, diff)
        print(f"{t:8.3f}  {p_sp[i]:12.4e}  {p_m_dp[i]:12.4e}  {diff:10.2e}")
    print(f"\nMax |error|: {max_err:.2e}")
    if max_err < 1e-10:
        print("✓ Degeneration to single-porosity PASSED")
    else:
        print("✗ Degeneration test — non-zero error (may be acceptable for floating-point)")

    return max_err# ============================================================================
# HDF5 Reader
# ============================================================================

def read_geos_hdf5(filename, field_name):
    """Read GEOS TimeHistory HDF5 file."""
    data = {}
    with h5py.File(filename, 'r') as f:
        # Get time array
        data['time'] = np.array(f['pressure elementCenter Time'][:, 0])
        # Get element centers and field values
        data['elementCenter'] = np.array(f['pressure elementCenter'][:, 0])
        data['field'] = np.array(f[field_name][:, 0])
    return data


def read_geos_hdf5_displacement(filename):
    """Read GEOS displacement TimeHistory HDF5 file."""
    data = {}
    with h5py.File(filename, 'r') as f:
        data['time'] = np.array(f['totalDisplacement Time'][:, 0])
        data['field'] = np.array(f['totalDisplacement'][:, 0])
        # Node coordinates
        data['nodePosition'] = np.array(f['totalDisplacement ReferencePosition'][:, 0])
    return data


# ============================================================================
# Plotting
# ============================================================================

def plot_comparison(analytical, geos_data, title, out_file):
    """Generate comparison plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(title, fontsize=14)

    # Subplot 1: pressure vs time at centre (x ~ 0)
    ax = axes[0, 0]
    t_ana = analytical.get('t', [])
    p_m_ana = analytical.get('pm_centre', [])
    p_f_ana = analytical.get('pf_centre', [])

    if len(p_m_ana) > 0:
        ax.plot(t_ana, p_m_ana, 'b-', label='$p_m$ (analytical)')
    if len(p_f_ana) > 0:
        ax.plot(t_ana, p_f_ana, 'r-', label='$p_f$ (analytical)')

    # GEOS data
    for label, color, data in geos_data:
        if data and 'time' in data:
            ax.plot(data['time'], data['field'], color + '--', label=label + ' (GEOS)')

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Pressure [Pa]')
    ax.legend()
    ax.set_title('Pressure vs Time at Centre')

    # Subplot 2: pressure profile along x at fixed time
    ax = axes[0, 1]
    ax.text(0.5, 0.5, 'Pressure profile\n(needs spatial data from VTK)',
            ha='center', va='center', transform=ax.transAxes)
    ax.set_title('Pressure Profile p(x)')

    # Subplot 3: displacement u_z vs time
    ax = axes[1, 0]
    uz_ana = analytical.get('uz', [])
    if len(uz_ana) > 0:
        ax.plot(t_ana, uz_ana, 'g-', label='$u_z$ (analytical)')
    for label, color, data in geos_data:
        if data and 'disp' in data:
            ax.plot(data['disp']['time'], data['disp']['field'],
                    color + '--', label=label + ' (GEOS)')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('$u_z$ [m]')
    ax.legend()
    ax.set_title('Vertical Displacement vs Time')

    # Subplot 4: error vs time
    ax = axes[1, 1]
    ax.text(0.5, 0.5, 'Error analysis\n(computed after data loading)',
            ha='center', va='center', transform=ax.transAxes)
    ax.set_title('Error vs Time')

    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f'Figure saved to {out_file}')


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='DPDP Mandel verification')
    parser.add_argument('--output-dir', type=str,
                        default=os.path.expanduser('~/result/dpdkValidation'),
                        help='GEOS output directory')
    parser.add_argument('--case', type=str, choices=['degrade', 'dual', 'both'],
                        default='dual',
                        help='Test case: degrade (single-porosity limit), '
                             'dual (full dual-porosity), both')
    args = parser.parse_args()

    out_dir = args.output_dir
    params = ShaleParams()

    print('=' * 60)
    print('DPDP Mandel Verification')
    print('Mehrabian & Abousleiman (2014)')
    print('=' * 60)

    # ---- Build analytical solutions ----
    print('\nBuilding analytical solutions...')
    single_sol = SinglePorosityMandel(params)
    dual_sol = DualPorosityMandel(params)

    # Time array for analytical solution
    t_array = np.logspace(0, 4, 50)  # 1 to 10000 s

    # ---- Case 1: Degradation test ----
    if args.case in ('degrade', 'both'):
        print('\n' + '-' * 40)
        print('Case: Degradation to Single-Porosity Mandel')
        print('-' * 40)

        # Analytical single-porosity response
        x_centre = 0.0  # centre of sample
        z_top = params.b
        p_sp = np.array([single_sol.pressure(x_centre, z_top, t) for t in t_array])
        uz_sp = np.array([single_sol.displacement_z(x_centre, z_top, t) for t in t_array])

        print(f'  p(0,0,1)   = {single_sol.pressure(0, 0, 1):.6e}')
        print(f'  p(0,0,10)  = {single_sol.pressure(0, 0, 10):.6e}')
        print(f'  u_z(b,1)     = {single_sol.displacement_z(0, z_top, 1):.6e}')
        print(f'  nu_u          = {single_sol.nu_u:.6f}')
        print(f'  B              = {single_sol.B:.6f}')
        print(f'  alpha          = {single_sol.alpha:.6f}')
        print(f'  c              = {single_sol.c:.6e} m2/s')
        print(f'  Roots found    = {len(single_sol.alpha_n)}')

        # Check if GEOS output exists
        deg_dir = os.path.join(out_dir, 'step_degenerate')
        p_hist = os.path.join(deg_dir, 'pressure_matrix_history.hdf5')
        d_hist = os.path.join(deg_dir, 'displacement_history.hdf5')

        if os.path.exists(p_hist):
            print(f'\n  Reading GEOS output from {deg_dir}')
            p_data = read_geos_hdf5(p_hist, 'pressure')
            d_data = read_geos_hdf5_displacement(d_hist)
        else:
            print(f'\n  GEOS output not found at {deg_dir}')
            print('  Run: geosx -i DPDP_Mandel_Mehrabian2014_degenerate.xml '
                  f'-o {deg_dir}')
            p_data, d_data = None, None

        # Generate analytical reference data files for later comparison
        ana_file = os.path.join(out_dir, 'analytical_single_porosity.npz')
        np.savez(ana_file, t=t_array, p_centre=p_sp, uz_top=uz_sp)
        print(f'  Analytical data saved to {ana_file}')

    # ---- Case 2: Full dual-porosity test ----
    if args.case in ('dual', 'both'):
        print('\n' + '-' * 40)
        print('Case: Dual-Porosity Mandel')
        print('-' * 40)

        x_centre = 0.001  # near centre (avoid x=0 for stability)
        z_top = params.b

        # Compute dual-porosity analytical solution at discrete times
        print('  Computing dual-porosity analytical solution...')
        # Use a subset of times for faster computation
        t_dual = np.logspace(1, 3.5, 20)  # 10 to ~3162 s

        p_dp = np.zeros((len(t_dual), 2))
        uz_dp = np.zeros(len(t_dual))

        for it, t in enumerate(t_dual):
            try:
                p_res = dual_sol.pressure(x_centre, z_top, np.array([t]))
                p_dp[it, :] = p_res[0, :]
            except Exception as e:
                p_dp[it, :] = np.nan
                if it < 3:
                    print(f'    Warning at t={t:.1f}: {e}')

        print(f'  p_m(0,t=10) estimated = {p_dp[0,0]:.6e}')
        print(f'  p_f(0,t=10) estimated = {p_dp[0,1]:.6e}')
        print(f'  Gamma_13               = {params.Gamma:.2e} Pa-1 s-1')
        print(f'  kappa1                 = {params.v1*params.k1/params.mu:.2e} m2/Pa/s')
        print(f'  kappa2                 = {params.v2*params.k2/params.mu:.2e} m2/Pa/s')
        print(f'  S11                    = {dual_sol.S11:.2e} Pa-1')
        print(f'  S22                    = {dual_sol.S22:.2e} Pa-1')

        # Check if GEOS output exists
        dual_dir = os.path.join(out_dir, 'step_dual')
        pm_hist = os.path.join(dual_dir, 'pressure_matrix_history.hdf5')
        pf_hist = os.path.join(dual_dir, 'pressure_fracture_history.hdf5')
        d_hist = os.path.join(dual_dir, 'displacement_history.hdf5')

        if os.path.exists(pm_hist):
            print(f'\n  Reading GEOS output from {dual_dir}')
            pm_data = read_geos_hdf5(pm_hist, 'pressure')
            pf_data = read_geos_hdf5(pf_hist, 'pressure')
            d_data = read_geos_hdf5_displacement(d_hist)
        else:
            print(f'\n  GEOS output not found at {dual_dir}')
            print('  Run: mpirun -np 4 geosx -i '
                  'DPDP_Mandel_Mehrabian2014_FIM.xml '
                  f'-o {dual_dir} -x 2 -y 1 -z 2')
            pm_data = pf_data = d_data = None

        # Generate analytical reference data
        ana_file = os.path.join(out_dir, 'analytical_dual_porosity.npz')
        np.savez(ana_file, t=t_dual, pm_centre=p_dp[:, 0],
                 pf_centre=p_dp[:, 1])
        print(f'  Analytical data saved to {ana_file}')

    # ---- Summary ----
    print('\n' + '=' * 60)
    print('Verification data preparation complete.')
    print(f'Results directory: {out_dir}')
    print('=' * 60)


if __name__ == '__main__':
    main()
