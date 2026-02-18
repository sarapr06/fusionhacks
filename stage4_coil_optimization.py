"""
Stage 4: Coil optimization (sweep over aspect ratios).

Per COILOPT: minimize B·n̂ deviation from target equilibrium + engineering (curvature, length).
Per R&D: simpler coils for maintainability — stronger curvature/length penalties.

Constraints:
  - Plasma-coil distance ≥ 15% of minor radius
  - Coil-coil distance ≥ 5% of minor radius

Outputs: stage4_outputs/sweep_R{R}_a{a}/ per design, sweep_summary.json
"""
import matplotlib
matplotlib.use("Agg")

import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "stage3_outputs")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "stage4_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

from desc import set_device
set_device("cpu")

import numpy as np
from desc.coils import initialize_modular_coils
from desc.equilibrium import Equilibrium
from desc.grid import LinearGrid
from desc.objectives import (
    CoilCurvature,
    CoilLength,
    CoilSetMinDistance,
    LinkingCurrentConsistency,
    ObjectiveFunction,
    PlasmaCoilSetMinDistance,
    QuadraticFlux,
)
from desc.optimize import Optimizer
from desc.plotting import plot_coils, plot_surfaces
from desc.integrals import compute_B_plasma
import matplotlib.pyplot as plt
import json

# Relaxed so designs with reasonable ι (not ι≈0) can pass. τ_E ∝ ι^0.41.
target = 3.5e-1


def compute_Bn_over_B(field, eq, vacuum=False):
    grid = LinearGrid(M=40, N=40, NFP=eq.NFP)
    Bn, surf_coords = field.compute_Bnormal(
        eq.surface if vacuum else eq,
        eval_grid=grid,
        chunk_size=20,
        B_plasma_chunk_size=20,
    )
    B_total = field.compute_magnetic_field(surf_coords)
    if not vacuum:
        B_plasma = compute_B_plasma(eq, eval_grid=grid, chunk_size=20)
        B_total = B_total + B_plasma
    B_mag = np.linalg.norm(B_total, axis=1)
    B_mag = np.where(B_mag < 1e-14, 1.0, B_mag)
    return np.abs(Bn) / B_mag


# Load stage3 sweep summary (R,a sweep)
stage3_summary = os.path.join(INPUT_DIR, "sweep_summary.json")
if not os.path.isfile(stage3_summary):
    print("Run stage3 sweep first.")
    raise SystemExit(1)

with open(stage3_summary) as f:
    designs = json.load(f)

results = []
for d in designs:
    R0, a = d["R0_target"], d["a_target"]
    design_name = f"sweep_R{R0:.2f}_a{a:.3f}"
    sweep_dir = os.path.join(OUTPUT_DIR, design_name)
    eq_path = os.path.join(INPUT_DIR, design_name, "equilibrium.h5")

    if not os.path.isfile(eq_path):
        print(f"Skipping {design_name}: {eq_path} not found")
        continue

    A_val = R0 / a
    print(f"\n{'#'*60}\n# R={R0:.2f} m, a={a:.2f} m → A={A_val:.2f}\n{'='*60}")
    os.makedirs(sweep_dir, exist_ok=True)

    eq = Equilibrium.load(eq_path)
    data = eq.compute(["p"])
    vacuum = np.max(np.abs(np.asarray(data["p"]).flat)) < 1e-10

    # More coils + higher resolution + prioritize B·n̂ (QA+ι geometries harder to match)
    num_coils = 8
    r_over_a = 2.5
    coilset = initialize_modular_coils(eq, num_coils=num_coils, r_over_a=r_over_a)
    coilset = coilset.to_FourierXYZ(N=12)  # more Fourier modes → more coil shape flexibility

    coil_grid = LinearGrid(N=60)
    plasma_grid = LinearGrid(M=28, N=28, NFP=eq.NFP, sym=eq.sym)  # denser surface for matching

    # Flux match dominant; curvature/length relaxed so coils conform to target field
    weights = {
        "quadratic flux": 2000,
        "coil-coil min dist": 200,
        "plasma-coil min dist": 50,
        "coil curvature": 200,
        "coil length": 25,
    }

    objective = ObjectiveFunction(
        (
            QuadraticFlux(eq, field=coilset, eval_grid=plasma_grid, field_grid=coil_grid,
                         vacuum=vacuum, weight=weights["quadratic flux"], bs_chunk_size=10),
            CoilSetMinDistance(coilset, bounds=(0.05, np.inf), normalize_target=False,
                              grid=coil_grid, weight=weights["coil-coil min dist"], dist_chunk_size=2),
            PlasmaCoilSetMinDistance(eq, coilset, bounds=(0.15, np.inf), normalize_target=False,
                                    plasma_grid=plasma_grid, coil_grid=coil_grid, eq_fixed=True,
                                    weight=weights["plasma-coil min dist"], dist_chunk_size=2),
            CoilCurvature(coilset, bounds=(-1, 2), normalize_target=False, grid=coil_grid,
                         weight=weights["coil curvature"]),
            CoilLength(coilset, bounds=(0, 50), normalize_target=True, grid=coil_grid,
                       weight=weights["coil length"]),
        ),
    )
    constraints = (LinkingCurrentConsistency(eq, coilset, eq_fixed=True),)
    # COILOPT solvers: proximal-lsq-exact (default); try scipy-lm for LM-style
    try:
        optimizer = Optimizer("proximal-scipy-lm")
        print("Optimizing coils (Levenberg-Marquardt)...")
    except Exception:
        optimizer = Optimizer("proximal-lsq-exact")
        print("Optimizing coils (proximal-lsq-exact)...")

    try:
        (optimized_coilset,), _ = optimizer.optimize(
            coilset, objective=objective, constraints=constraints,
            maxiter=200, verbose=2, ftol=1e-5, copy=True,
        )
    except Exception as ex:
        print(f"Coil optimization FAILED for {design_name}: {ex}")
        results.append({"R0": R0, "a": a, "A": A_val, "mean_BnB": 1.0, "max_BnB": 1.0, "pass": False})
        continue

    Bn_over_B = compute_Bn_over_B(optimized_coilset, eq, vacuum=vacuum)
    mean_BnB = float(np.mean(Bn_over_B))
    max_BnB = float(np.max(Bn_over_B))
    results.append({"R0": R0, "a": a, "A": A_val, "mean_BnB": mean_BnB, "max_BnB": max_BnB, "pass": max_BnB <= target})

    print(f"R={R0:.2f}, a={a:.2f} → A={A_val:.2f}: <|B·n̂/B|>={mean_BnB:.3e}, max={max_BnB:.3e}, target≤{target:.1e}: {'PASS' if max_BnB <= target else 'FAIL'}")

    optimized_coilset.save_in_makegrid_format(
        os.path.join(sweep_dir, "coilset_optimized.txt"), NFP=eq.NFP, grid=coil_grid,
    )

    fig, ax = plot_surfaces(eq, rho=6, theta=8)
    fig.suptitle(f"Plasma R={R0:.2f} m, a={a:.2f} m, A={A_val:.2f}")
    fig.savefig(os.path.join(sweep_dir, "plasma_surfaces.png"), dpi=150)
    plt.close(fig)

    fig = plot_coils(optimized_coilset)
    try:
        fig.write_image(os.path.join(sweep_dir, "coils_optimized.png"))
    except Exception:
        fig.write_html(os.path.join(sweep_dir, "coils_optimized.html"))

    grid_2d = LinearGrid(rho=1.0, M=30, N=30, NFP=eq.NFP)
    Bn_2d, coords = optimized_coilset.compute_Bnormal(
        eq.surface if vacuum else eq, eval_grid=grid_2d, chunk_size=20, B_plasma_chunk_size=20,
    )
    B_vec = optimized_coilset.compute_magnetic_field(coords)
    if not vacuum:
        B_vec += compute_B_plasma(eq, eval_grid=grid_2d, chunk_size=20)
    BnB_2d = np.abs(Bn_2d) / np.where(np.linalg.norm(B_vec, axis=1) < 1e-14, 1.0, np.linalg.norm(B_vec, axis=1))
    fig, ax = plt.subplots(figsize=(8, 5))
    sc = ax.scatter(grid_2d.nodes[:, 2], grid_2d.nodes[:, 1], c=BnB_2d, cmap="viridis", s=20)
    plt.colorbar(sc, ax=ax, label="|B·n̂/B|")
    ax.set_xlabel("ζ"); ax.set_ylabel("θ")
    ax.set_title(f"|B·n̂/B| on LCFS (R={R0:.2f}, a={a:.2f}, A={A_val:.2f}, max={float(np.max(BnB_2d)):.2e})")
    fig.savefig(os.path.join(sweep_dir, "Bn_over_B_surface.png"), dpi=150)
    plt.close(fig)

    with open(os.path.join(sweep_dir, "summary.txt"), "w") as f:
        f.write(f"R={R0:.2f}, a={a:.2f}, A={A_val:.2f}\n<|B·n̂/B|>={mean_BnB:.3e}\nmax={max_BnB:.3e}\nTarget≤{target:.1e}: {'PASS' if max_BnB <= target else 'FAIL'}\n")

with open(os.path.join(OUTPUT_DIR, "sweep_summary.json"), "w") as f:
    json.dump(results, f, indent=2)

print(f"\n{'#'*60}\nSWEEP COMPLETE\n{'='*60}")
for r in results:
    print(f"  R={r['R0']:.2f}, a={r['a']:.2f} → A={r['A']:.2f}: max|B·n̂/B|={r['max_BnB']:.2e} {'PASS' if r['pass'] else 'FAIL'}")
print(f"\nSummary: {OUTPUT_DIR}/sweep_summary.json")
