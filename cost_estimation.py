"""
FusionHacks 2026 official cost function.

From https://github.com/IssraAli/fusionhacks2026/blob/main/cost_estimation.py

cost = baseline * (c_v * (V/V0)^1.1 + c_coils * (c_l*|L/L0|^1.2 + c_kt*(κτ term)^1.4) + c_p * P_ext/p0)
"""
import numpy as np


def _scalar(x):
    arr = np.asarray(x)
    return float(arr.flat[0]) if arr.size else float("nan")


def total_reactor_cost(eq, coilset, P_ext_MW):
    """
    Compute initial reactor cost for FusionHacks 2026.

    :param eq: DESC Equilibrium
    :param coilset: DESC CoilSet (iterable of coils)
    :param P_ext_MW: Heating power [MW]
    :returns: Cost [$]
    """
    vol = _scalar(eq.compute("V")["V"])
    total_current_length = 0
    mean_curvature_list = []
    mean_torsion_list = []
    max_curvature_list = []
    max_torsion_list = []

    for coil in coilset:
        L = _scalar(coil.compute("length")["length"])
        total_current_length += float(coil.current) * L
        kappa = np.asarray(coil.compute("curvature")["curvature"]).flatten()
        tau = np.asarray(coil.compute("torsion")["torsion"]).flatten()
        mean_curvature_list.append(float(np.mean(kappa**2)))
        mean_torsion_list.append(float(np.mean(tau**2)))
        max_curvature_list.append(float(np.max(kappa**2)))
        max_torsion_list.append(float(np.max(tau**2)))

    mean_k2 = np.mean(mean_curvature_list)
    mean_t2 = np.mean(mean_torsion_list)
    max_k2 = np.max(max_curvature_list)
    max_t2 = np.max(max_torsion_list)

    vol0 = 50
    l0 = 1e7
    kt0 = 5
    p0 = 10e6  # 10 MW

    c_v = 0.4
    c_coils = 0.5
    c_l = 0.6
    c_kt = 0.4
    f_k = 0.99
    f_t = 1e-3
    c_p = 0.1

    baseline = 5e8  # $500M

    cost_kt = np.sqrt(f_k * mean_k2 + f_t * mean_t2) + np.sqrt(
        f_k * max_k2 + f_t * max_t2
    )
    cost = baseline * (
        c_v * (vol / vol0) ** 1.1
        + c_coils
        * (
            c_l * np.abs(total_current_length / l0) ** 1.2
            + c_kt * (cost_kt / kt0) ** 1.4
        )
        + c_p * (P_ext_MW * 1e6) / p0  # p0=10e6 W; P_ext_MW in MW
    )
    return float(cost)
