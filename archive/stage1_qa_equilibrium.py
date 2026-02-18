"""
Run DESC, get one converged QA equilibrium.

Solves MHD (J x B = grad p) by minimizing the force residual using
a Fourier-Zernike spectral expansion.

Output: equilibrium + plots (saved to stage1_outputs/).
"""
import os
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "stage1_outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

from desc import set_device
set_device("cpu")
import numpy as np
from desc.continuation import solve_continuation_automatic
from desc.equilibrium import Equilibrium
from desc.geometry import FourierRZToroidalSurface
from desc.plotting import plot_surfaces, plot_boundary, plot_boozer_surface

# Initial surface: aspect ratio ~6 (Landreman & Paul), circular cross section with slight torsion
surf = FourierRZToroidalSurface(
    R_lmn=[1, 0.166, 0.1],
    Z_lmn=[-0.166, -0.1],
    modes_R=[[0, 0], [1, 0], [0, 1]],
    modes_Z=[[-1, 0], [0, -1]],
    NFP=2,
)

# M=4, N=4 for speed (increase to 8,8 for higher accuracy)
eq = Equilibrium(M=4, N=4, Psi=0.087, surface=surf)
print("Solving MHD equilibrium (minimizing force residual via Fourier-Zernike)...")
eq = solve_continuation_automatic(eq, objective="force", bdry_step=0.5, verbose=2)[-1]
print("Equilibrium converged.")
# Save equilibrium
eq_path = os.path.join(OUTPUT_DIR, "stage1_qa_equilibrium.h5")
eq.save(eq_path)
print(f"Saved {eq_path}")

# Plot
import matplotlib.pyplot as plt
# 1. Flux surfaces (default phi = 6 slices for non-axisymmetric)
fig, ax = plot_surfaces(eq, rho=6, theta=8)
fig.suptitle("Flux surfaces (nested pressure surfaces)")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stage1_flux_surfaces.png"), dpi=150)
plt.close(fig)
print(f"Saved {OUTPUT_DIR}/stage1_flux_surfaces.png")
# 2. Boundary (LCFS)
fig, ax = plot_boundary(eq)
fig.suptitle("Boundary (LCFS)")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stage1_boundary_lcfs.png"), dpi=150)
plt.close(fig)
print(f"Saved {OUTPUT_DIR}/stage1_boundary_lcfs.png")
# 3. |B| contours in Boozer coords
fig, ax = plot_boozer_surface(eq, rho=1.0, fill=True, ncontours=25)
fig.suptitle("|B| in Boozer coords (QA check)")
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stage1_B_boozer.png"), dpi=150)
plt.close(fig)
print(f"Saved {OUTPUT_DIR}/stage1_B_boozer.png")
