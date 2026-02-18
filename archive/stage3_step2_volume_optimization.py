"""
Phase 3 Step 2: Volume optimization.

Loads QA-optimized equilibrium from stage3_outputs, adds volume objective,
re-optimizes while preserving QA quality. Saves to stage3_outputs.

f(x) = QA error + Volume (target higher) + AspectRatio
"""
import matplotlib
matplotlib.use("Agg")

import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "stage3_qa_outputs")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "stage3_qa_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

from desc import set_device
set_device("cpu")

import numpy as np
from desc.equilibrium import EquilibriaFamily, Equilibrium
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
    Volume,
)
from desc.optimize import Optimizer
from desc.plotting import plot_surfaces, plot_boundary, plot_boozer_surface
import matplotlib.pyplot as plt

# Load QA-optimized equilibrium from stage3_step1
eq_path = os.path.join(INPUT_DIR, "stage3_qa_optimized.h5")
print(f"Loading {eq_path} ...")
eq = Equilibrium.load(eq_path)

# Current volume (for target)
data = eq.compute("V")
V_current = float(data["V"])
V_target = V_current * 1.15  # Target 15% higher volume
print(f"Current plasma volume: {V_current:.6f} m^3")
print(f"Target plasma volume:  {V_target:.6f} m^3 (+15%)")

eqfam = EquilibriaFamily(eq)
grid = LinearGrid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=np.array([0.6, 0.8, 1.0]), sym=True)
A_target = 6.0

k = 2
objective = ObjectiveFunction(
    (
        QuasisymmetryTwoTerm(eq=eq, helicity=(1, 0), grid=grid, normalize=False),
        Volume(eq=eq, target=V_target, weight=0.5, normalize=False),
        AspectRatio(eq=eq, target=A_target, weight=1e1, normalize=False),
        RotationalTransform(eq=eq, target=0.42, weight=10, normalize=False),
    ),
)
R_modes = np.vstack(([0, 0, 0], eq.surface.R_basis.modes[np.max(np.abs(eq.surface.R_basis.modes), 1) > k, :]))
Z_modes = eq.surface.Z_basis.modes[np.max(np.abs(eq.surface.Z_basis.modes), 1) > k, :]
constraints = (
    ForceBalance(eq=eq),
    FixBoundaryR(eq=eq, modes=R_modes),
    FixBoundaryZ(eq=eq, modes=Z_modes),
    FixPressure(eq=eq),
    FixCurrent(eq=eq),
    FixPsi(eq=eq),
)
eq_new, _ = eq.optimize(
    objective=objective, constraints=constraints,
    optimizer=Optimizer("proximal-lsq-exact"),
    maxiter=50, verbose=1, copy=True,
    options={"initial_trust_radius": 0.5},
)
eq = eq_new

data = eq.compute("V")
V_final = float(data["V"])
eq.save(os.path.join(OUTPUT_DIR, "stage3_qa_volume_optimized.h5"))
print(f"\nSaved {OUTPUT_DIR}/stage3_qa_volume_optimized.h5")
print(f"Volume change: {(V_final - V_current) / V_current * 100:.2f}%")

# Plots
fig, ax = plot_surfaces(eq, rho=6, theta=8)
fig.suptitle("Flux surfaces (volume optimized)")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stage3_step2_flux_surfaces.png"), dpi=150)
plt.close(fig)
print(f"Saved {OUTPUT_DIR}/stage3_step2_flux_surfaces.png")

fig, ax = plot_boundary(eq)
fig.suptitle("Boundary (LCFS)")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stage3_step2_boundary_lcfs.png"), dpi=150)
plt.close(fig)
print(f"Saved {OUTPUT_DIR}/stage3_step2_boundary_lcfs.png")

fig, ax = plot_boozer_surface(eq, rho=1.0, fill=True, ncontours=25)
fig.suptitle("|B| in Boozer coords")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stage3_B_boozer.png"), dpi=150)
plt.close(fig)
print(f"Saved {OUTPUT_DIR}/stage3_B_boozer.png")
