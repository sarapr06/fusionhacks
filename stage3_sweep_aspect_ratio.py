"""
Phase 3 R–a Sweep (geometry for QA + ISS04).

Sweep over (R, a) — a grows WAY faster than R, so we sample a more densely/widely.
For each (R, a): create surface, solve MHD, optimize QA first, then non-QA geometry.

  Order: QA geometry first, then LOCK IT. Coils and V vs A are optimized later (stage4+).
  1. Create initial surface with R0=R, minor radius a
  2. Solve MHD equilibrium
  3. QA optimization (QuasisymmetryTwoTerm primary) — geometry for confinement
  4. LOCK: use QA-optimized geometry as-is. No further boundary deformation.
  5. Save V, R0, a, ι. Coils (stage4) and V/A tradeoffs come after.

ISS04: τ_E ∝ B^0.84 R^0.64 a^2.28 ... — geometry from this (R,a) sweep.
"""
import matplotlib
matplotlib.use("Agg")

import os
import json

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "stage3_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

from desc import set_device
set_device("cpu")

import numpy as np
from desc.continuation import solve_continuation_automatic
from desc.equilibrium import EquilibriaFamily, Equilibrium
from desc.geometry import FourierRZToroidalSurface
from desc.grid import LinearGrid
from desc.objectives import (
    AspectRatio,
    FixBoundaryR,
    FixBoundaryZ,
    FixCurrent,
    FixPressure,
    FixPsi,
    ForceBalance,
    ObjectiveFunction,
    QuasisymmetryTwoTerm,
    RotationalTransform,
)
from desc.optimize import Optimizer
from desc.plotting import plot_surfaces, plot_boundary, plot_boozer_surface
import matplotlib.pyplot as plt

# Sweep R and a: R varies slowly, a varies more (a grows way faster than R)
R_VALS = [1.0]  # major radius (m); keep narrow
A_VALS = np.linspace(0.083, 0.52, 12)  # minor radius (m); a from 0.083 (A~12) to 0.52 (A~2), 12 points

# Build (R,a) pairs
SWEEP_POINTS = [(R, a) for R in R_VALS for a in A_VALS]
results = []

for R0_target, a_target in SWEEP_POINTS:
    A_target = R0_target / a_target
    print(f"\n{'#'*60}")
    print(f"# R={R0_target:.3f} m, a={a_target:.3f} m → A=R/a={A_target:.2f}")
    print("="*60)

    sweep_dir = os.path.join(OUTPUT_DIR, f"sweep_R{R0_target:.2f}_a{a_target:.3f}")
    os.makedirs(sweep_dir, exist_ok=True)

    # Initial surface: R0 = major radius, a = minor radius (base coeffs ~0.2 give a~0.17)
    scale_a = a_target / 0.17
    surf = FourierRZToroidalSurface(
        R_lmn=[R0_target, 0.166 * scale_a, 0.1 * scale_a],
        Z_lmn=[-0.166 * scale_a, -0.1 * scale_a],
        modes_R=[[0, 0], [1, 0], [0, 1]],
        modes_Z=[[-1, 0], [0, -1]],
        NFP=2,
    )

    # 1. Solve MHD equilibrium
    print("\n[1/4] Solving MHD equilibrium...")
    eq = Equilibrium(M=4, N=4, Psi=0.087, surface=surf)
    eq = solve_continuation_automatic(eq, objective="force", bdry_step=0.5, verbose=1)[-1]

    # 2. QA optimization — geometry for confinement first (per hackathon: QA first)
    #    Drive ∇f_QA → 0 before any volume/flux terms. QS is primary (weight 1e2).
    print("\n[2/4] QA optimization (QuasisymmetryTwoTerm primary)...")
    eqfam = EquilibriaFamily(eq)
    grid = LinearGrid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=np.array([0.6, 0.8, 1.0]), sym=True)

    for k in range(1, min(eq.M + 1, 4)):
        # QA primary, but ι must not be ≈0 (τ_E ∝ ι^0.41). Strong enough RotationalTransform
        # weight so ι≈0.42 is jointly optimized with QA before lock. (200 caused geometry drift)
        objective = ObjectiveFunction(
            (
                QuasisymmetryTwoTerm(eq=eqfam[-1], helicity=(1, 0), grid=grid, normalize=False, weight=1e2),
                AspectRatio(eq=eqfam[-1], target=A_target, weight=5, normalize=False),
                RotationalTransform(eq=eqfam[-1], target=0.42, weight=50, normalize=False),
            ),
        )
        R_modes = np.vstack(([0, 0, 0], eq.surface.R_basis.modes[np.max(np.abs(eq.surface.R_basis.modes), 1) > k, :]))
        Z_modes = eq.surface.Z_basis.modes[np.max(np.abs(eq.surface.Z_basis.modes), 1) > k, :]
        constraints = (
            ForceBalance(eq=eqfam[-1]),
            FixBoundaryR(eq=eqfam[-1], modes=R_modes),
            FixBoundaryZ(eq=eqfam[-1], modes=Z_modes),
            FixPressure(eq=eqfam[-1]),
            FixCurrent(eq=eqfam[-1]),
            FixPsi(eq=eqfam[-1]),
        )
        eq_new, _ = eqfam[-1].optimize(
            objective=objective, constraints=constraints,
            optimizer=Optimizer("proximal-lsq-exact"),
            maxiter=50, verbose=1, copy=True,
            options={"initial_trust_radius": 0.5},
        )
        eqfam.append(eq_new)

    eq = eqfam[-1]

    # 3. LOCK — QA geometry from step 2 is final. No further boundary deformation.
    #    Coils (stage4) and V vs A tradeoffs are optimized afterward.
    print("\n[3/4] QA geometry locked (no refinement).")

    # Metrics (some quantities are arrays on grid; take first value)
    def _scalar(x):
        arr = np.asarray(x)
        return float(arr.flat[0]) if arr.size else float("nan")

    data = eq.compute(["V", "a_major/a_minor", "R0", "a", "iota"])
    V = _scalar(data["V"])
    A_actual = _scalar(data["a_major/a_minor"])
    R0 = _scalar(data["R0"])
    a = _scalar(data["a"])
    # ι at rho=2/3 for ISS04: τ_E ∝ ι_2/3^0.41
    iota_data = eq.compute("iota", grid=LinearGrid(rho=np.array([2.0 / 3.0]), M=2, N=2, NFP=eq.NFP))
    iota_2_3 = float(np.asarray(iota_data["iota"]).flat[0])

    results.append({
        "R0_target": R0_target,
        "a_target": a_target,
        "A": R0_target / a_target,
        "A_actual": A_actual,
        "V": V,
        "R0": R0,
        "a": a,
        "iota_2_3": iota_2_3,
    })

    # Save
    # 4. Save
    eq.save(os.path.join(sweep_dir, "equilibrium.h5"))
    for name, fn in [
        ("flux_surfaces", lambda: plot_surfaces(eq, rho=6, theta=8)),
        ("boundary", lambda: plot_boundary(eq)),
        ("B_boozer", lambda: plot_boozer_surface(eq, rho=1.0, fill=True, ncontours=25)),
    ]:
        fig, ax = fn()
        fig.suptitle(f"R={R0_target:.2f} m, a={a_target:.2f} m, A=R/a={A_target:.2f}")
        fig.tight_layout()
        fig.savefig(os.path.join(sweep_dir, f"{name}.png"), dpi=150)
        plt.close(fig)

    print(f"\nR={R0_target:.2f}, a={a_target:.2f}, A=R/a={A_target:.2f}: V={V:.6f} m³, R0={R0:.4f}, a={a:.4f}, ι={iota_2_3:.4f}")

# Summary
summary_path = os.path.join(OUTPUT_DIR, "sweep_summary.json")
with open(summary_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{'#'*60}")
print("R–a SWEEP COMPLETE (a sampled widely; R fixed)")
print("="*60)
print(f"\nResults saved to {OUTPUT_DIR}/")
for r in results:
    print(f"  R={r['R0_target']:.2f}, a={r['a_target']:.2f} → A={r['A']:.2f}: V={r['V']:.6f} m³")
print(f"\nSummary: {summary_path}")
print("\nNext: stage5 computes A=R/a and τ_E for each (R,a).")
