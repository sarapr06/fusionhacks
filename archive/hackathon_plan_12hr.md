# 12-Hour Hackathon Plan
## QA Stellarator Volumetric Neutron Source for Mo-99 Production

**Problem:** Optimize a quasi-axisymmetric stellarator design as a volumetric neutron source to produce medical isotope Molybdenum-99.

**Weekend theme: Optimizing volumes and physics.** Optimize **flux surface geometry** via **stage-one optimization**. x = boundary Fourier coefficients (R^b_mn, Z^b_mn). **Order:** optimize for quasi-symmetry first, *then* add surface/flux and volume considerations. At the optimum, ∇_x f(x) = 0.

---

## Theory & Equations (from Stellarator Physics Primer)

Use these to ground your analysis and justify design choices. Each equation is explained below.

### Fusion & Confinement

**Lawson criterion** (Eq. 2): `nT τ_E ≥ 3×10²¹ keV·s/m³`
- *What it does:* n = plasma density, T = temperature, τ_E = energy confinement time. The product nTτ_E must exceed ~3×10²¹ for net fusion gain.
- *How it helps:* In Phase 2, use this for neutron scaling. More fusion → more neutrons from higher n, T, τ_E. Cite it when arguing QA improves τ_E (Eq. 16) and thus neutron yield.

**MHD equilibrium** (Eq. 7): `J × B = ∇p`
- *What it does:* Force balance: Lorentz force balances pressure gradient. B and J lie on nested flux surfaces; B·∇p = 0.
- *How it helps:* This is what DESC solves. In Phase 1, flux surfaces come from Eq. 7. Say "We solve MHD equilibrium to obtain nested flux surfaces" in your summary.

**Gyroradius** (Eq. 5): `ρ_c = mv⊥/(qB)`
- *What it does:* Radius of particle orbit around a field line. Hotter or lighter particles, weaker B → larger ρ_c.
- *How it helps:* Confinement requires ρ_c ≪ machine size. Use when sizing: low B or high T means particles drift out. Mention in business case for why compact designs need sufficient B.

### DESC / Stellarator Optimization

**DESC solve** (Eq. 9): `min ||J×B − ∇p||`
- *What it does:* DESC finds Fourier-Zernike coefficients (Eq. 8) that minimize the force residual.
- *How it helps:* Every `solve_continuation_automatic` call minimizes this. Report "force residual < X" as a quality metric in your technical summary.

**Stage-one optimization** (Eq. 11–13)
- *What it does:* Minimize **f(x)** where **x** = boundary coefficients (R^b_mn, Z^b_mn). f(x) measures QA quality (e.g., quasi-symmetry error) or volume. Constraints: force balance J×B=∇p, B·n̂=0 on boundary. At optimum, **∇_x f(x) = 0** — zero gradient of the target w.r.t. the degrees of freedom that control the boundary shape. You optimize flux surface geometry by varying x.
- *How it helps:* `QuasisymmetryTwoTerm` defines f(x) as QA error; volume/aspect-ratio terms add to f(x). The optimizer drives ∇_x f → 0 while satisfying constraints. Core of the weekend theme: volumes + physics optimization.

**Stage-two optimization** (Eq. 14)
- *What it does:* Find coils that match the target boundary (B·n̂≈0 on LCFS) with engineering penalties (curvature, length).
- *How it helps:* `initialize_modular_coils` starts coil design. Cite Eq. 14 when discussing manufacturability.

### Quasi-Axisymmetry (QA)

**Three QA objectives in DESC (Part 3 paper):** (1) Boozer f_B — expensive, global. (2) Two-term f_C — `QuasisymmetryTwoTerm`, local, needs helicity. (3) Triple product f_T — `QuasisymmetryTripleProduct`, no helicity, local, volume; **Part 3: ~10× better confinement** in their test.

**Quasisymmetry**: `B = B(ρ, η)` with `η = Mϑ_B − Nφ_B`
- *What it does:* In Boozer coords, B depends on one angle η. QA = M=1, N=0 → B = B(ρ, φ_B), tokamak-like.
- *How it helps:* `QuasisymmetryTwoTerm(helicity=(1,0))` enforces QA. Plot |B| in Boozer coords; contours should look axisymmetric. Use in slides to show "our design is QA."

**Neoclassical diffusion** (Eq. 16): `D ∼ ε_eff^(3/2) v_d²/ν`
- *What it does:* In 1/ν regime, transport scales with v_d (drift speed) and ν (collisionality). Higher v_d → worse confinement.
- *How it helps:* QA minimizes v_d → smaller D → better τ_E → better Lawson → more neutrons. Main argument: "We optimize QA to reduce neoclassical transport and improve neutron yield." Connect QA → τ_E → Lawson.

**Trapped particles** (Eq. 15): `E = ½mv∥² + μB`, `μ = mv⊥²/(2B)`
- *What it does:* Particles can't enter B > E/μ; they bounce (trapped). In non-QA stellarators, trapped orbits drift radially (v_d ≠ 0).
- *How it helps:* Background for Eq. 16. QA removes the drift. If asked "why QA?", say: "Trapped particles normally drift out; QA confines them, reducing D and improving τ_E."

### Neutron Source → Mo-99

D–D or D–T fusion → neutrons. Thermalize in blanket → irradiate Mo-98 → Mo-99. V and τ_E (via Lawson) set neutron yield.

---

## Time Allocation (12 hours total)

| Phase | Hours | Focus |
|-------|-------|-------|
| **1. QA Equilibrium** | 2–3 hr | Get a working QA design |
| **2. Neutron Model** | 2 hr | Simple scaling + Mo-99 estimate |
| **3. Optimization** | 2–3 hr | Tweak for neutron-relevant metrics |
| **4. Business Case** | 2–3 hr | Economics + value proposition |
| **5. Deliverables** | 2 hr | Slides, summary, plots |

---

## Phase 1: QA Equilibrium (2–3 hr)

**Goal:** Run DESC and obtain one converged QA stellarator equilibrium.

**Physics:** DESC solves the MHD equilibrium (Primer Eq. 7: J×B=∇p) by minimizing force residual (Eq. 9) via the Fourier-Zernike ansatz (Eq. 8). The result gives nested flux surfaces where B·∇p=0.

**Steps:**
1. Copy and run `precise_QA.py` from DESC examples.
2. Use CPU if no GPU: `set_device("cpu")` instead of `"gpu"`.
3. Reduce resolution if needed: try `M=4, N=4` for speed.
4. Save output: `eqfam.save("qa_base.h5")`.
5. Generate 2–3 plots: flux surfaces (nested pressure surfaces), boundary (LCFS), |B| contours in Boozer coords (QA quality check).

**Reference:** `desc/examples/precise_QA.py`, [Basic QS Optimization](https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/basic_optimization.html), Primer §3.2–3.2.1

---

## Phase 2: Neutron Model (2 hr)

**Goal:** Connect stellarator design to Mo-99 production with a simple model.

**Theory link:** Neutron yield scales with fusion reactions. Use **Lawson criterion** (Primer Eq. 2): nTτ_E ≥ 3×10²¹ keV·s/m³. For a volumetric source, neutron rate ∝ plasma volume × (fusion reactivity) ∝ V·n²·⟨σv⟩. Better confinement (τ_E) from QA (Primer §4.2) enables higher nT product.

**Mo-99 basics:**
- Reaction: Mo-98(n,γ) → Mo-99
- Needs thermal/epithermal neutrons (≈10¹³–10¹⁴ n/cm²/s for useful yield)
- ~1 Ci/g Mo-99 possible with optimized targets in research reactors

**Simple model:**
1. Compute plasma volume V from equilibrium (DESC `eq.compute()` or flux-surface geometry).
2. Use scaling: neutron source ≈ V × n² × ⟨σv⟩ (or literature values for D-D/D-T yield per m³).
3. Assume Mo-98 targets in blanket; thermal flux scales with source and moderator.
4. Estimate Mo-99 (Ci/day) from flux × target area × σ(Mo-98) × time.

**Output:** One number for Mo-99 production rate + short list of assumptions. Reference Primer Eq. 2 when stating n, T, τ_E.

---

## Phase 3: Optimization (2–3 hr)

**Goal:** Stage-one optimization — flux surface geometry. **Order:** (1) optimize for quasi-symmetry first, (2) then add volume and flux surface considerations.

**Formulation:** **min f(x)** s.t. constraints. **x** = (R^b_mn, Z^b_mn). Constraints: force balance J×B=∇p, B·n̂=0 on LCFS.

**Step 1 — Quasi-symmetry first:** f(x) = QA error only. Options: `QuasisymmetryTwoTerm` (helicity=(1,0) for QA) or `QuasisymmetryTripleProduct` (Part 3: ~10× better confinement, no helicity). Run to convergence; save equilibrium.

**Step 2 — Volume & flux surfaces:** From QA equilibrium, add volume and flux-surface terms. f(x) = QA error + volume + aspect ratio. Don't sacrifice QA quality.

**DESC optimization (Part 3):** c = boundary coeffs R^b_mn, Z^b_mn. Equilibrium x* from force balance; optimizer c* = arg min |g(x*,c)|². Gauss–Newton trust-region + AD (single equilibrium per step, exact derivatives).

**Reference:** [Advanced QS Optimization](https://desc-docs.readthedocs.io/en/stable/notebooks/tutorials/advanced_optimization.html), Primer §4, §4.2

---

## Phase 4: Business Case (2–3 hr)

**Goal:** One-page summary and comparison vs existing Mo-99 supply.

**Content:**
1. **Problem:** Mo-99 supply relies on HEU fission reactors; supply risks, proliferation.
2. **Solution:** Compact QA stellarator → no HEU, distributed production.
3. **Technical summary:** Design params, plasma volume, neutron yield, Mo-99 (Ci/week). *Cite Lawson (Eq. 2) and QA confinement (Eq. 16) in the narrative.*
4. **Economics:** Capex vs fission; cost per Ci (order-of-magnitude).
5. **Value:** No HEU, scalable, regional, reduced regulatory complexity.

---

## Phase 5: Deliverables (2 hr)

1. **DESC outputs:** Saved equilibrium(s), 3–5 key plots.
2. **1-page technical summary:** Design, neutron model (with Primer eqs.), Mo-99 estimate, assumptions.
3. **1-page business case:** Problem, solution, economics, value proposition.
4. **Pitch deck (5–7 slides):** Problem → Solution → Design → Results → Business case → Next steps.

---

## Quick Reference: DESC Commands

```python
# Solve equilibrium (minimizes Eq. 9: ||J×B - ∇p||)
from desc.continuation import solve_continuation_automatic
eq = solve_continuation_automatic(eq, objective="force", bdry_step=0.5)[-1]

# QA optimization (enforces B=B(ρ,η), M=1,N=0)
from desc.objectives import QuasisymmetryTwoTerm, AspectRatio, ObjectiveFunction
objective = ObjectiveFunction((QuasisymmetryTwoTerm(helicity=(1,0), ...), AspectRatio(...)))

# Coils (stage-two, Eq. 14)
from desc.coils import initialize_modular_coils
coils = initialize_modular_coils(eq, n_coils=6, ...)
```

---

## Hackathon Priorities

1. **Working equilibrium first** – One converged QA design is enough.
2. **Simple neutron model** – Use scaling + Lawson (Eq. 2); skip full neutronics.
3. **Story over precision** – Plausible numbers + clear narrative + Primer equations.
4. **Presentation-ready** – Save plots as PNG/PDF; prepare slides early.

---

**References:** (2)stellarator-physics-primer.pdf, the-desc-stellarator-code-suite-part-3-quasi-symmetry-optimization.pdf. DESC docs: https://desc-docs.readthedocs.io/*
