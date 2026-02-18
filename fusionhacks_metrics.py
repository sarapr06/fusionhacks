"""
FusionHacks 2026 metrics — temperature and neutron fluence from DESC equilibria.

From https://github.com/IssraAli/fusionhacks2026/blob/main/fusionhacks_metrics.py

- temp_from_eq(eq, P_ext): solve power balance for T [keV]
- neutron_fluence(eq, T): total neutron production rate [neutrons/s]
"""
from desc.equilibrium import Equilibrium
from scipy.optimize import fsolve
import numpy as np
from desc.grid import LinearGrid
from desc.integrals import Bounce2D

mu0 = (4 * np.pi) * 1e-7  # permeability of free space
keV_to_J = 1.6022e-16
E_alpha = 3.5e3 * keV_to_J  # in J
BETA = 0.05


def sigmav(T):
    """D-T reactivity ⟨σv⟩ in m³/s. T in keV. Bosch-Hale fit."""
    T = np.asarray(T, dtype=float)
    BG = 34.3827
    mr_c2 = 1124656
    C1, C2, C3 = 1.17302e-9, 1.51361e-2, 7.61886e-2
    C4, C5, C6, C7 = 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5
    numerator = T * (C2 + T * (C4 + T * C6))
    denominator = 1.0 + T * (C3 + T * (C5 + T * C7))
    theta = T / (1.0 - (numerator / denominator))
    xi = (BG**2 / (4.0 * theta)) ** (1.0 / 3.0)
    prefactor = C1 * theta * np.sqrt(xi / (mr_c2 * T**3))
    sigma_v = prefactor * np.exp(-3.0 * xi)
    return 1e-6 * sigma_v


def eps_avg(eq):
    """RMS effective ripple ε for H-factor model."""
    rho = np.linspace(0.01, 1, 5)
    grid = LinearGrid(rho=rho, M=eq.M, N=eq.N, NFP=eq.NFP, sym=False)
    X, Y = 16, 32
    theta = Bounce2D.compute_theta(eq, X, Y, rho)
    Y_B = 32
    num_transit = 20
    num_well = 10 * num_transit
    num_quad = 32
    num_pitch = 45
    data = eq.compute(
        "effective ripple",
        grid=grid,
        theta=theta,
        Y_B=Y_B,
        num_transit=num_transit,
        num_well=num_well,
        num_quad=num_quad,
        num_pitch=num_pitch,
    )
    eps = grid.compress(data["effective ripple"])
    return np.sqrt(np.mean(eps**2))


def h_factor(eq, verbose=False):
    """H-factor from ε. h_min=0.7, h_max=2, μ=2, k=2."""
    h_max, h_min = 2, 0.7
    mu, k = 2, 2
    eps = eps_avg(eq)
    q = eps / 0.1
    h = h_min + (h_max - h_min) / (1 + np.exp(k * (np.log(q) + mu)))
    if verbose:
        print(f"eps = {eps}")
        print(f"H factor = {h}")
    return float(h)


def tau_E_iss04(h, B, R0, a, n, iota_23, P_ext, verbose=False):
    """τ_E [s]. P_ext in [W]."""
    n0 = n / 1e19
    p0 = P_ext / 1e6
    tau_E = h * 0.465 * (B**0.84) * (R0**0.64) * (a**2.28) * (n0**0.54) * (p0**-0.61) * (float(iota_23) ** 0.41)
    if verbose:
        print("B =", B, "R =", R0, "a =", a, "n19 =", n / 1e19, "PMW =", P_ext / 1e6, "iota =", iota_23, "h =", h, "tau_E =", tau_E)
    return tau_E


def pb_res(T, beta, B_on_axis, R0, a, iota_23, h, vol, P_ext):
    """Power balance residual [MW]. T in keV."""
    n = (beta * B_on_axis**2) / (2 * mu0 * T * keV_to_J)
    tau_E = tau_E_iss04(h, B_on_axis, R0, a, n, iota_23, P_ext)
    res = P_ext + (n**2 / 4) * sigmav(T) * vol * E_alpha - 3 * (n * T * vol * keV_to_J) / tau_E
    return res / 1e6


def _scalar(x):
    """Extract scalar from possibly array-valued DESC compute output."""
    arr = np.asarray(x)
    return float(arr.flat[0]) if arr.size else float("nan")


def temp_from_eq(eq, P_ext, verbose=False):
    """
    Solve power balance for T [keV] given equilibrium and external power.

    :param eq: DESC Equilibrium
    :param P_ext: External power [W]
    :param verbose: Pass to h_factor and fsolve prints
    :returns: Temperature [keV] (scalar)
    """
    a = _scalar(eq.compute("a")["a"])
    R0 = _scalar(eq.compute("R0")["R0"])
    B_on_axis = _scalar(eq.compute("<|B|>_axis")["<|B|>_axis"])
    vol = _scalar(eq.compute("V")["V"])
    rho = [2.0 / 3.0]
    iota_23 = _scalar(eq.compute("iota", grid=LinearGrid(rho=rho))["iota"])
    h = h_factor(eq, verbose=verbose)
    args = (BETA, B_on_axis, R0, a, iota_23, h, vol, P_ext)
    t_solved = fsolve(pb_res, x0=5.0, args=args)
    T_keV = float(t_solved[0])
    if verbose:
        print(f"T: {T_keV} keV")
    return T_keV


def neutron_fluence(eq, T, verbose=False):
    """
    Total neutron production rate [neutrons/s] from D-T fusion.

    :param eq: DESC Equilibrium
    :param T: Temperature [keV]
    :param verbose: Print density
    :returns: Fluence [neutrons/s]
    """
    B_on_axis = _scalar(eq.compute("<|B|>_axis")["<|B|>_axis"])
    vol = _scalar(eq.compute("V")["V"])
    n = (BETA * B_on_axis**2) / (2 * mu0 * T * keV_to_J)
    if verbose:
        print(f"density: {n}")
    fluence = 0.25 * (n**2) * sigmav(T) * vol
    return fluence
