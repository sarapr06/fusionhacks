# FusionHacks 2026 — QA Stellarator Volumetric Neutron Source

Design of a quasi-axisymmetric (QA) stellarator for Mo-99 medical isotope production. Final design: R=2 m, a=0.22 m, ι≈0.1, max|B·n̂/B|=4.08×10⁻³ (passes), cost **$0.71B**.

---

## Final Design Summary

| Parameter | Value |
|-----------|-------|
| R₀ | 2.00 m |
| a | 0.22 m |
| Aspect ratio A | 9.1 |
| Volume V | 1.92 m³ |
| ι₂/₃ | 0.10 |
| max\|B·n̂/B\| | 4.08×10⁻³ ✓ PASS |
| Coils | Modular, B·n̂ ≤ 5×10⁻³ |
| Cost | $0.71B |
| P_ext | 10 MW |
| R_neutrons | 8.97×10¹³ /s |

---

## Outputs (Complete Inventory)

### Final design (`final/`)

| File | Description |
|------|-------------|
| `equilibrium.h5` | DESC vacuum MHD equilibrium |
| `coilset.txt` | Coil filaments (makegrid format) |
| `submission.h5` | Bundled submission (eq + coilset + P_ext) — run `python3 submission.py 3` for `group_3.h5` |
| `final_summary.json` | Design metrics (R, a, cost, ι, B·n̂, neutron rate, Mo-99, tritium) |
| `REFINED_DESIGN_RESULTS.md` | ι refinement rationale, Harris (2004) citations |
| `ECONOMIC_INPUTS.md` | Economic parameters, COGS, references |
| `plasma_surfaces.png` | Flux surfaces plot |
| `coils.html` | Interactive 3D coil viewer |
| `README.md` | This directory overview |

### Stage outputs

| Directory | Contents |
|-----------|----------|
| `stage3_outputs/` | QA equilibria per (R,a) — regenerate with `stage3_sweep_aspect_ratio.py` (not in git) |
| `stage4_outputs/` | Coil optimization outputs — regenerate with `stage4_coil_optimization.py` (not in git) |
| `stage4_outputs/sweep_summary.json` | R, a, A, mean/max B·n̂, pass/fail for all designs |
| `stage5_outputs/` | tau_E_vs_A.png, cost_vs_A.png, R_neutrons_vs_A.png, geometry_for_planner.json |
| `best_qa/` | Original QA design before ι refinement |
| `final_refined/` | Alternative output from ι refinement runs |

### Scripts

| Script | Purpose |
|--------|----------|
| `run_pipeline.sh` | Full pipeline: stage3 → stage4 → stage5 → stage6 |
| `run_fast.sh` | Fast pipeline (A∈[3,9]): stage3 → stage4 → stage5 |
| `run_rerun_final.sh` | ι refinement (ι→0.1) on best QA |
| `submission.py` | Bundle eq + coilset + P_ext → submission.h5 |
| `stage3_sweep_aspect_ratio.py` | R–a sweep, QA optimization |
| `stage4_coil_optimization.py` | Coil optimization for each design |
| `stage5_reactor_optimization.py` | τ_E, cost, neutron rate vs A |
| `stage6_finalize_geometry.py` | Export to final/, optional ι refinement |
| `power_balance_solver.py` | T from power balance (ISS04 scaling) |
| `fusionhacks_metrics.py` | temp_from_eq, neutron_fluence |
| `cost_estimation.py` | total_reactor_cost (V, coils, P_ext) |

---

## References

- [FusionHacks 2026 challenge](https://github.com/IssraAli/fusionhacks2026)
- Harris (2004), *Small to mid-sized stellarator experiments* — [ANU](https://people.physics.anu.edu.au/~jnh112/AIIM/My%20publications/Harris_eps.pdf)
