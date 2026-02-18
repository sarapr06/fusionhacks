"""
Phase 3 Step 1: QA optimization.

Loads equilibrium from stage1_outputs, minimizes quasi-symmetry error
while preserving force balance. Saves to stage3_outputs.
"""
import matplotlib
matplotlib.use("Agg")

import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "stage1_outputs")
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
)
from desc.optimize import Optimizer
from desc.plotting import plot_surfaces, plot_boundary, plot_boozer_surface
import matplotlib.pyplot as plt

eq_path = os.path.join(INPUT_DIR, "stage1_qa_equilibrium.h5")
if not os.path.isfile(eq_path):
    print(f"Run stage1 first to create {eq_path}")
    raise SystemExit(1)

print(f"Loading {eq_path} ...")
eq = Equilibrium.load(eq_path)

eqfam = EquilibriaFamily(eq)
grid = LinearGrid(M=eq.M, N=eq.N, NFP=eq.NFP, rho=np.array([0.6, 0.8, 1.0]), sym=True)
A_target = 6.0  # aspect ratio target

# Optimize in steps (like precise_QA.py): gradually free more boundary modes
for k in range(1, min(eq.M + 1, 4)):
    objective = ObjectiveFunction(
        (
            QuasisymmetryTwoTerm(eq=eqfam[-1], helicity=(1, 0), grid=grid, normalize=False),
            AspectRatio(eq=eqfam[-1], target=A_target, weight=1e1, normalize=False),
            RotationalTransform(eq=eqfam[-1], target=0.42, weight=10, normalize=False),
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

eqfam.save(os.path.join(OUTPUT_DIR, "stage3_qa_family.h5"))
eq_out = os.path.join(OUTPUT_DIR, "stage3_qa_optimized.h5")
eqfam[-1].save(eq_out)
print(f"\nSaved {eq_out}")
