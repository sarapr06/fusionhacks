# Design Principles (NCSX, COILOPT, R&D Priorities)

Reference slides: NCSX "Designed for Attractive Properties," COILOPT coil optimization flowchart, ARIES-CS/ReNeW "R&D Needs."

---

## 1. Quasi-axisymmetric (QA) geometry

- **Principle:** Design for quasi-axisymmetric magnetic field (B ≈ B(ρ, φ)).
- **Our implementation:** `QuasisymmetryTwoTerm(helicity=(1,0))` in stage3; QA optimized first, then non-QA geometry.
- **Status:** ✅ Aligned

---

## 2. High β (4–6%)

- **Principle:** NCSX passively stable at β≈4.1%; can reach β>6% with coil current adjustments.
- **Our implementation:** β=5% limit in power balance; vacuum equilibria (leader instruction); operating point uses β=5% for n.
- **Status:** ✅ Power balance uses β=5%. Equilibrium is vacuum (geometry design). For finite-β stability we would need a separate stability study.

---

## 3. Steady state

- **Principle:** Stellarators are inherently steady state; no inductive current drive.
- **Our implementation:** Fixed-boundary MHD equilibria; vacuum fields; time-independent.
- **Status:** ✅ Aligned

---

## 4. Coil optimization (COILOPT)

- **Principle:** Minimize deviation of B-normal on plasma surface from target equilibrium; include engineering (curvature, length, spacing).
- **Our implementation:** Stage4: QuadraticFlux, CoilCurvature, CoilLength, CoilSetMinDistance, PlasmaCoilSetMinDistance; target B·n̂/B ≤ 5×10⁻³.
- **Status:** ✅ Aligned
- **Slide solvers (LM, Differential Evolution, Genetic):** DESC uses proximal-lsq-exact; could explore alternatives if needed.

---

## 5. Simpler coils / maintainability

- **Principle:** ARIES-CS: "Simpler coils, improving maintenance access"; "Systematic study to simplify designs."
- **Our implementation:** Modular coils (4 coils); CoilCurvature and CoilLength in objectives to favor simpler shapes.
- **Status:** ⚠️ Partial — could add explicit simplicity/maintenance metrics (e.g., fewer coils, lower curvature).

---

## 6. High-performance integration

- **Principle:** ReNeW: "Validate high β at low collisionality, high confinement on QS experiment."
- **Our implementation:** QA for confinement; ISS04 scaling; power balance with β=5%.
- **Status:** ✅ Aligned for the current scope

---

## 7. Stability (kink, ballooning, vertical, Mercier, etc.)

- **Principle:** NCSX designed for passive stability against kink, ballooning, vertical, Mercier, neoclassical-tearing modes.
- **Our implementation:** Not run; would require stability modules (e.g., DESC stability tools).
- **Status:** ❌ Not implemented — future work

---

## 8. Flexible coils

- **Principle:** "By adjusting currents can control stability, transport, shape: iota, shear."
- **Our implementation:** Coil currents optimized implicitly; fixed equilibrium boundary. Free-boundary with coil currents as DoFs would go further.
- **Status:** ⚠️ Partial — fixed-boundary only

---

## Checklist for ongoing work

| Principle           | In mind? | Where                          |
|--------------------|----------|--------------------------------|
| QA first           | ✅       | stage3                         |
| High β (5%)        | ✅       | power_balance_solver           |
| Steady state       | ✅       | Vacuum MHD equilibria          |
| B·n̂ match          | ✅       | stage4 QuadraticFlux           |
| Engineering (κ, L) | ✅       | stage4 CoilCurvature, CoilLength |
| Mercier stability  | ✅       | stage3, stage6                 |
| Simpler coils      | ⚠️       | Consider fewer coils / simpler shapes |
| Ballooning         | ⚠️       | Optional in stage6 (REFINE=1, STAGE6_BALLOONING=1) |
| Flexible coils     | ⚠️       | Fixed-boundary only            |

## Geometry finalization (stage6)

- **Final design**: R≈1 m, a≈0.45 m, A≈2.3 — coils pass B·n̂/B ≤ 5×10⁻³.
- **Output**: `final/equilibrium.h5`, `final/coilset.txt`, `final/final_summary.json`
- **Run**: `python3 stage6_finalize_geometry.py` (default: export best as-is)
- **Refinement**: `REFINE=1` to attempt ι improvement (may degrade coil match)
