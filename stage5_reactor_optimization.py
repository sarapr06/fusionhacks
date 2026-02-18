"""
Stage 5: Reactor optimization — τ_E, power balance (β=5%), D-T neutron rate, cost.

Uses FusionHacks 2026 metrics (temp_from_eq, neutron_fluence) when equilibrium available.
Falls back to power_balance_solver for designs without equilibrium files.
Cost: official FusionHacks cost_estimation when eq+coilset available; else placeholder.
D-T fusion → 14 MeV neutrons → Mo-100(n,2n)Mo-99 for medical isotope.
"""
import json
import os

import matplotlib.pyplot as plt
import numpy as np

from power_balance_solver import solve_power_balance, mo99_production_annual_Ci

try:
    from desc.coils import MixedCoilSet
    from desc.equilibrium import Equilibrium
    from fusionhacks_metrics import temp_from_eq, neutron_fluence
    from cost_estimation import total_reactor_cost
    FUSIONHACKS_AVAILABLE = True
except ImportError:
    FUSIONHACKS_AVAILABLE = False

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
STAGE3_DIR = os.path.join(SCRIPT_DIR, "stage3_outputs")
STAGE4_DIR = os.path.join(SCRIPT_DIR, "stage4_outputs")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "stage5_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Operational assumptions for τ_E (placeholders; planner will refine)
B_T = 1.0       # Tesla
n_20 = 1.0      # 10^20 m^-3
P_ext_MW = 10.0 # MW
H_factor = 1.0  # confinement enhancement

# Placeholder when eq+coilset not available
COST_MIN_B = 0.7
COST_MAX_B = 5.0

P_EXT_MW = 10.0


def cost_placeholder_B(V, V_min, V_max):
    """Placeholder reactor cost in $B. Scales with V^1.2 (matches challenge formula)."""
    v12_min = V_min ** 1.2
    v12_max = V_max ** 1.2
    v12 = V ** 1.2
    frac = (v12 - v12_min) / (v12_max - v12_min) if v12_max > v12_min else 0
    return COST_MIN_B + (COST_MAX_B - COST_MIN_B) * frac


def tau_E_ISS04(R, a, iota, B=B_T, n_20=n_20, P_ext_MW=P_ext_MW, H=H_factor):
    """τ_E in seconds. Per challenge: R,a in m; B in T; n_20 in 10^20 m^-3; P_ext_MW in MW."""
    iota_safe = max(abs(iota), 0.01)  # avoid singularity
    return H * 0.465 * (B ** 0.84) * (R ** 0.64) * (a ** 2.28) * (n_20 ** 0.54) * (P_ext_MW ** -0.61) * (iota_safe ** 0.41)


def main():
    stage3 = os.path.join(STAGE3_DIR, "sweep_summary.json")
    if not os.path.isfile(stage3):
        print("Run stage3 sweep first.")
        return

    if FUSIONHACKS_AVAILABLE:
        print("Using FusionHacks 2026 metrics (temp_from_eq, neutron_fluence) for designs with equilibrium files.\n")
    else:
        print("FusionHacks metrics unavailable; using power_balance_solver fallback.\n")

    with open(stage3) as f:
        s3 = json.load(f)

    V_vals = [r["V"] for r in s3]
    V_min, V_max = min(V_vals), max(V_vals)

    P_ext_W = P_ext_MW * 1e6
    results = []
    for r in s3:
        R0, a = r["R0"], r["a"]
        R0_t, a_t = r.get("R0_target", R0), r.get("a_target", a)
        A = R0 / a
        iota = r.get("iota_2_3") or 0.42
        V = r["V"]
        tau_E = tau_E_ISS04(R0, a, iota)
        cost_B = cost_placeholder_B(V, V_min, V_max)

        T_keV, R_neutrons, n_20, Q_mo99 = None, None, None, None
        eq_path = os.path.join(STAGE3_DIR, f"sweep_R{R0_t:.2f}_a{a_t:.3f}", "equilibrium.h5")

        if FUSIONHACKS_AVAILABLE and os.path.isfile(eq_path):
            try:
                eq = Equilibrium.load(eq_path)
                T_keV = temp_from_eq(eq, P_ext_W, verbose=False)
                R_neutrons = neutron_fluence(eq, T_keV, verbose=False)
                if R_neutrons is not None:
                    B = float(eq.compute("<|B|>_axis")["<|B|>_axis"])
                    n = (0.05 * B**2) / (2 * 4e-7 * np.pi * T_keV * 1.602e-16)
                    n_20 = n / 1e20
                    Q_mo99 = mo99_production_annual_Ci(R_neutrons, f_mo99=0.5)
            except Exception as ex:
                print(f"  FusionHacks failed for R={R0:.2f} a={a:.2f}: {ex}")

        if T_keV is None or R_neutrons is None:
            sol = solve_power_balance(B_T, R0, a, iota, V, P_ext_MW, H_factor)
            n, T_keV, R_neutrons = sol
            if n is not None and R_neutrons is not None:
                n_20 = n / 1e20
                Q_mo99 = mo99_production_annual_Ci(R_neutrons, f_mo99=0.5)
            else:
                n_20, T_keV, R_neutrons, Q_mo99 = None, None, None, None

        # Cost: official FusionHacks when eq+coilset available
        cost_B = cost_placeholder_B(V, V_min, V_max)
        if FUSIONHACKS_AVAILABLE and os.path.isfile(eq_path):
            coil_path = os.path.join(STAGE4_DIR, f"sweep_R{R0_t:.2f}_a{a_t:.3f}", "coilset_optimized.txt")
            if os.path.isfile(coil_path):
                try:
                    eq_for_cost = Equilibrium.load(eq_path)
                    coilset = MixedCoilSet.from_makegrid_coilfile(coil_path, ignore_groups=True)
                    cost_USD = total_reactor_cost(eq_for_cost, coilset, P_EXT_MW)
                    cost_B = cost_USD / 1e9
                except Exception:
                    pass

        results.append({
            "R0": R0,
            "a": a,
            "A": A,
            "V": V,
            "iota_2_3": iota,
            "tau_E_s": tau_E,
            "cost_B": round(cost_B, 3),
            "n_20": round(n_20, 4) if n_20 is not None else None,
            "T_keV": round(T_keV, 2) if T_keV is not None else None,
            "R_neutrons_s": round(R_neutrons, 2) if R_neutrons is not None else None,
            "Q_Mo99_Ci_year": round(Q_mo99, 2) if Q_mo99 is not None else None,
        })

    with open(os.path.join(OUTPUT_DIR, "geometry_for_planner.json"), "w") as f:
        json.dump(results, f, indent=2)

    # τ_E vs A plot — swept (R,a), A = R/a, τ_E from ISS04
    fig, ax = plt.subplots(figsize=(8, 5))
    A_vals = [x["A"] for x in results]
    tau_vals = [x["tau_E_s"] for x in results]
    ax.plot(A_vals, tau_vals, "o-")
    for r in results:
        ax.annotate(f"R={r['R0']:.2f}\na={r['a']:.2f}", (r["A"], r["tau_E_s"]), fontsize=7, xytext=(5, 5),
                    textcoords="offset points", ha="left")
    ax.set_xlabel("A = R/a (from R,a sweep; a varies widely)")
    ax.set_ylabel("τ_E (s)")
    ax.set_title(f"τ_E vs A: sweep (R,a) → A=R/a → τ_E (a grows faster than R)\n(B={B_T} T, n={n_20}e20 m⁻³, P_ext={P_ext_MW} MW)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "tau_E_vs_A.png"), dpi=150)
    plt.close(fig)

    # Cost vs A (placeholder; R vs C tradeoff proxy)
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    A_vals = [x["A"] for x in results]
    cost_vals = [x["cost_B"] for x in results]
    ax2.plot(A_vals, cost_vals, "s-", color="C1")
    for r in results:
        ax2.annotate(f"${r['cost_B']:.2f}B", (r["A"], r["cost_B"]), fontsize=7, xytext=(5, 5),
                    textcoords="offset points", ha="left")
    ax2.set_xlabel("A = R/a")
    ax2.set_ylabel("Cost ($B)")
    ax2.set_title(f"Reactor cost vs A (FusionHacks 2026 cost_estimation)")
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(OUTPUT_DIR, "cost_vs_A.png"), dpi=150)
    plt.close(fig2)

    best_tau = max(results, key=lambda x: x["tau_E_s"])
    worst_tau = min(results, key=lambda x: x["tau_E_s"])

    print("τ_E sweep (ISS04) complete.")
    print(f"Assumptions: B={B_T} T, n={n_20}e20 m⁻³, P_ext={P_ext_MW} MW, H={H_factor}")
    print(f"\nSwept (R,a); a sampled widely. A = R/a.")
    print(f"Max τ_E = {best_tau['tau_E_s']:.2f} s at R={best_tau['R0']:.2f} m, a={best_tau['a']:.2f} m → A={best_tau['A']:.2f}")
    print(f"Min τ_E = {worst_tau['tau_E_s']:.2f} s at R={worst_tau['R0']:.2f} m, a={worst_tau['a']:.2f} m → A={worst_tau['A']:.2f}")
    with_R = [r for r in results if r.get("R_neutrons_s") is not None]
    if with_R:
        fig3, ax3 = plt.subplots(figsize=(8, 5))
        ax3.plot([x["A"] for x in with_R], [x["R_neutrons_s"] for x in with_R], "o-", color="C2")
        ax3.set_xlabel("A = R/a")
        ax3.set_ylabel("R (neutrons/s, D-T)")
        ax3.set_title("D-T neutron rate vs A (β=5%, power balance)")
        ax3.set_yscale("log")
        ax3.grid(True, alpha=0.3)
        fig3.tight_layout()
        fig3.savefig(os.path.join(OUTPUT_DIR, "R_neutrons_vs_A.png"), dpi=150)
        plt.close(fig3)

    print(f"\nOutput: {OUTPUT_DIR}/geometry_for_planner.json, tau_E_vs_A.png, cost_vs_A.png" +
          (", R_neutrons_vs_A.png" if with_R else ""))
    print(f"\nPower balance: D-T fusion, β=5%, Mo-100(n,2n)Mo-99")
    print(f"Cost: FusionHacks 2026 total_reactor_cost (V, L, κτ, P_ext)\n")
    for r in results:
        ext = f", R={r['R_neutrons_s']:.2e}/s, Q_Mo99={r['Q_Mo99_Ci_year']} Ci/yr" if r.get("R_neutrons_s") else ""
        print(f"  R={r['R0']:.2f}, a={r['a']:.2f} → A={r['A']:.2f}: τ_E={r['tau_E_s']:.4f} s, V={r['V']:.4f} m³, cost≈${r['cost_B']:.2f}B{ext}")


if __name__ == "__main__":
    main()
