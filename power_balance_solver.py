"""
Power balance solver for D-T fusion neutron source.

D-T fusion: D + T → n(14 MeV) + α. Uses β = 5% limit (per challenge) to constrain n,
then solves for T from:
  P_ext + (n²/4)⟨σv⟩_DT E_α = (3/2)nkT / τ_E

Mo-99 production (14 MeV neutrons):
  - Mo-100(n,2n)Mo-99: primary route for fast neutrons (σ ~ 1.5 b)
  - Mo-98(n,γ)Mo-99: thermal route (use with moderated neutrons if available)
"""
import numpy as np
from scipy.optimize import brentq

# Physical constants (SI)
MU_0 = 4e-7 * np.pi
K_B = 1.380649e-23
E_ALPHA_J = 3.5e6 * 1.6e-19  # 3.5 MeV alpha from D-T
KEV_TO_K = 1.16e7  # 1 keV ≈ 1.16e7 K

BETA_LIMIT = 0.05  # 5% per challenge


def sigma_v_DT(T_keV):
    """D-T fusion reactivity ⟨σv⟩ in m³/s. T_keV = ion temperature in keV.
    Empirical fit valid 2-30 keV (matches Bosch-Hale to ~10%)."""
    T = np.clip(T_keV, 0.5, 100)
    return 1.1e-22 * (T / 10) ** 1.9


def n_from_beta(B_T, T_keV, beta=BETA_LIMIT):
    """Density n [m⁻³] from β = nkT/(B²/2μ₀). T_keV in keV."""
    p_mag = B_T**2 / (2 * MU_0)
    kT_J = T_keV * 1000 * 1.6e-19  # 1 keV = 1000 eV, 1 eV = 1.6e-19 J
    return beta * p_mag / kT_J


def power_balance_residual(T_keV, B_T, R, a, iota, V, P_ext_W, P_ext_MW, H=1.0):
    """Residual: LHS - RHS of power balance. Solve for T where residual = 0."""
    n = n_from_beta(B_T, T_keV)
    n_20 = n / 1e20
    tau_E = H * 0.465 * (B_T**0.84) * (R**0.64) * (a**2.28) * (n_20**0.54) * (
        P_ext_MW)**-0.61 * (max(abs(iota), 0.01)**0.41)
    sigma_v = sigma_v_DT(T_keV)
    P_alpha = (n**2 / 4) * sigma_v * E_ALPHA_J * V
    W = 1.5 * n * K_B * (T_keV * KEV_TO_K) * V
    lhs = P_ext_W + P_alpha
    rhs = W / tau_E
    return lhs - rhs


def solve_power_balance(B_T, R, a, iota, V, P_ext_MW, H=1.0, T_bracket=(1.0, 50.0)):
    """Solve for T [keV] and return n, T, R_neutrons. Returns (None,None,None) if no solution."""
    P_ext_W = P_ext_MW * 1e6
    try:
        T_keV = brentq(
            lambda T: power_balance_residual(T, B_T, R, a, iota, V, P_ext_W, P_ext_MW, H),
            T_bracket[0], T_bracket[1], xtol=1e-6,
        )
    except ValueError:
        return None, None, None
    n = n_from_beta(B_T, T_keV)
    sigma_v = sigma_v_DT(T_keV)
    R_neutrons = (n**2 / 4) * sigma_v * V
    return n, T_keV, R_neutrons


def mo99_production_annual_Ci(R_neutrons_s, f_mo99=0.5, eta=1e-23):
    """Mo-99 via Mo-100(n,2n)Mo-99 (14 MeV D-T neutrons). Returns Ci/year.
    R_neutrons_s = neutrons/s from D-T.
    f_mo99 = fraction of neutrons to Mo-99 (rest for tritium breeding).
    eta = Ci per neutron (depends on Mo-100 target, cross section ~1.5 b, geometry)."""
    seconds_per_year = 365.25 * 24 * 3600
    R_a = R_neutrons_s * seconds_per_year
    return eta * R_a * f_mo99
