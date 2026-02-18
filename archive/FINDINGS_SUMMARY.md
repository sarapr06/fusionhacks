# FusionHacks 2026 — Current Parameters & Findings

**Run date:** Stage 5 rerun completed. Stage 4 coil results from prior run (8% coil-coil; 5% constraint now in code).

---

## 1. Sweep parameters

| Parameter | Value |
|-----------|-------|
| R (major radius) | 1.0 m (fixed) |
| a (minor radius) | 12 points: 0.083 → 0.52 m |
| A = R/a | 12.0 → 1.9 |
| ISS04 inputs | B = 1 T, n = 1×10²⁰ m⁻³, P_ext = 10 MW, H = 1.0 |

---

## 2. τ_E (energy confinement time) — ISS04

| Design | R (m) | a (m) | A | V (m³) | ι_2/3 | τ_E (s) |
|--------|-------|-------|---|--------|-------|---------|
| Min τ_E | 1.00 | 0.08 | 12.05 | 0.14 | 0.00002 | 0.00006 |
| ... | 1.00 | 0.12 | 8.15 | 0.30 | 0.009 | 0.0001 |
| ... | 1.02 | 0.17 | 6.16 | 0.55 | 0.42 | 0.0013 |
| ... | 1.03 | 0.21 | 4.95 | 0.88 | 0.42 | 0.0023 |
| ... | 1.05 | 0.25 | 4.13 | 1.34 | 0.42 | 0.0036 |
| ... | 1.04 | 0.29 | 3.55 | 1.78 | 0.42 | 0.0050 |
| ... | 1.06 | 0.34 | 3.11 | 2.45 | 0.42 | 0.0072 |
| ... | 1.12 | 0.41 | 2.77 | 3.65 | 0.41 | 0.0109 |
| ... | 1.12 | 0.45 | 2.49 | 4.45 | 0.42 | 0.0139 |
| ... | 1.14 | 0.50 | 2.27 | 5.65 | 0.42 | 0.0180 |
| ... | 1.03 | 0.50 | 2.08 | 5.01 | 0.42 | 0.0165 |
| **Max τ_E** | **1.03** | **0.54** | **1.92** | **5.87** | **0.42** | **0.0197** |

---

## 3. Coil field error (B·n̂/B at ρ=1)

**Pass criterion:** mean & max both ≤ 5×10⁻³

| Design (A) | mean_BnB | max_BnB | Pass? |
|------------|----------|---------|-------|
| 12.0 | 0.0011 | 0.0032 | ✅ |
| 8.1 | 0.0008 | 0.0026 | ✅ |
| 6.2 | 0.0032 | 0.026 | ❌ |
| 4.9 | 0.0022 | 0.016 | ❌ |
| 4.1 | 0.0015 | 0.014 | ❌ |
| 3.6 | 0.0018 | 0.011 | ❌ |
| 3.1 | 0.028 | 0.99 | ❌ |
| 2.8 | 0.027 | 0.31 | ❌ |
| 2.5 | 0.061 | 0.44 | ❌ |
| 2.3 | 0.028 | 0.64 | ❌ |
| 2.1 | 0.024 | 0.22 | ❌ |
| 1.9 | 0.039 | 0.83 | ❌ |

**Result:** Only 2/12 designs meet the coil requirement (both low-A, small-a).

---

## 4. Main findings

1. **ι_2/3 drives τ_E at small a**
   - For a ≈ 0.08–0.12 m, ι_2/3 ≈ 0–0.01 → τ_E ~ 10⁻⁴ s.
   - For a ≳ 0.17 m, ι_2/3 ≈ 0.42 → τ_E rises with a.

2. **τ_E increases as A decreases**
   - ISS04: τ_E ∝ a^2.28, so larger minor radius increases τ_E.
   - Best τ_E ≈ 0.02 s at R ≈ 1.03 m, a ≈ 0.54 m, A ≈ 1.9.

3. **Coil–plasma match vs aspect ratio**
   - Designs with A ≈ 12 and 8 pass B·n̂/B ≤ 5×10⁻³.
   - For A ≲ 6, coils degrade (mean/max B·n̂/B exceed target).

4. **Tradeoff**
   - Best confinement: A ≈ 1.9, large a, large V.
   - Best coils: A ≈ 8–12, small a, small V.
   - No design satisfies both constraints; need more coil/geometry optimization or different design points.

---

## 5. Placeholder cost (0.7–5 B)

Per hackathon leader: use 0.7–5 B range until actual cost function is released. Scales with V^1.2.

| Design (A) | V (m³) | cost_placeholder_B |
|------------|--------|--------------------|
| 12.0 | 0.14 | 0.70 |
| 8.1 | 0.30 | 0.77 |
| ... | ... | ... |
| 1.9 | 5.87 | 5.00 |

---

## 6. Output files

| File | Content |
|------|---------|
| `stage3_outputs/sweep_summary.json` | R, a, V, ι_2/3 per design |
| `stage4_outputs/sweep_summary.json` | mean_BnB, max_BnB, pass flag |
| `stage5_outputs/geometry_for_planner.json` | Full geometry + τ_E + cost_placeholder_B |
| `stage5_outputs/tau_E_vs_A.png` | τ_E vs aspect ratio plot |
| `stage5_outputs/cost_vs_A.png` | Placeholder cost vs A |

---

## 7. Next steps

1. **Re-optimize coils** for high-τ_E designs (A ~ 2–3) — current coils fail.
2. **Implement power balance** to solve for T and reaction rate R.
3. **Implement cost C** from coils (L, curvature, V).
4. **H(ε)** — compute H from effective ripple for more realistic τ_E.
5. **Business case** — Mo‑99 market analysis and NPV.
