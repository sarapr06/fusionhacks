# FusionHacks 2026 Challenge Compliance Check

Reference: [fusionhacks2026_challenge.ipynb](https://github.com/IssraAli/fusionhacks2026/blob/main/fusionhacks2026_challenge.ipynb)

---

## ✅ What we have correctly

### 1. Equilibrium
| Requirement | Status |
|-------------|--------|
| QA stellarator (QuasisymmetryTwoTerm) | ✅ stage3 uses `QuasisymmetryTwoTerm(helicity=(1,0))` |
| Do NOT optimize ε directly | ✅ Use QuasisymmetryTwoTerm as proxy |
| Geometric constraints (aspect ratio) | ✅ AspectRatio objective, sweep (R,a) |
| Don't optimize R directly | ✅ Pick B, R, a, P_ext after optimization |
| MHD equilibrium (J×B=∇p) | ✅ solve_continuation_automatic, objective="force" |
| ι at ρ=2/3 for ISS04 | ✅ Saved as `iota_2_3` |

### 2. Coils
| Requirement | Status |
|-------------|--------|
| Minimize B·n̂/B on ρ=1 | ✅ QuadraticFlux objective |
| Target ~5×10⁻³ (penalized if exceeded) | ✅ target=5e-3, reported in sweep |
| Plasma-coil distance ≥ 15% of minor radius | ✅ PlasmaCoilSetMinDistance bounds=(0.15, ∞) |
| Coil-coil distance ≥ 5% of minor radius | ⚠️ We use 8% (stricter); challenge says 5% |

### 3. ISS04
| Requirement | Status |
|-------------|--------|
| τ_E = H × 0.465 × B^0.84 × R^0.64 × a^2.28 × n^0.54 × P_ext^-0.61 × ι^0.41 | ✅ Implemented in stage5 |

---

## ❌ Gaps (need hackathon planner / starter code)

### 1. Reaction rate R
- **Challenge:** R = (n/4)⟨σv⟩V
- **Note:** Challenge may use D-T; we had D-D. Formula differs (D-T: R ∝ n²⟨σv⟩V for 50-50 fuel).
- **Status:** Deferred – need challenge starter code for ⟨σv⟩(T) and power-balance solver.

### 2. Power balance & self-consistent T, n
- **Challenge:** P_ext + (n²/4)⟨σv⟩E_α = (3/2)nkT / τ_E
- **β limit:** n from β = nkT/(B²/2μ_0) = 5%
- **Status:** Not implemented – challenge says "We implement a function in python below"; we need that starter code.

### 3. H factor from ε
- **Challenge:** H = H_min + (H_max−H_min) / (1 + exp(−k(log(ε/ε_0)−μ))), ε_0=0.1, μ=4, H_min=0.7, H_max=2
- **Status:** We use H=1.0 fixed. Need ε (effective ripple) from equilibrium to compute H.

### 4. Cost formula
- **Challenge:** C = w_L L + w_V V^1.2 + w_κ(⟨κ⟩ + w_κ,m|κ_max|)^1.7
- **L** = total **current-length** [A·m] (NOT just length!)
- **Status:** We deferred cost. Formula differs from our earlier C = w_L L + w_V V^1.2 + curvature.

### 5. Business case
- Market analysis, Mo-99 supply/demand, Q_Mo99 = η × R_a × f
- 10-yr financial model, NPV
- **Status:** Not started.

---

## Units check for ISS04

Challenge states:
- n in **10^20 m⁻³**
- P_ext in **MW**

**Fixed:** stage5 now uses `n_20` (dimensionless) and `P_ext_MW` directly in the formula, per ISS04 convention.

---

## Action items

1. **Obtain starter code** from fusionhacks2026 repo (power balance solver, ⟨σv⟩, H(ε)).
2. **Implement H(ε)** from effective ripple (or proxy from QA error).
3. **Implement cost** per C = w_L L + w_V V^1.2 + w_κ(⟨κ⟩ + w_κ,m|κ_max|)^1.7 with L = current×length.
4. **Coil-coil 5%**: Change from 8% to 5% if challenge strictly requires 5%.
5. **Reaction rate & power balance**: Wire in starter code to solve T, then R.
