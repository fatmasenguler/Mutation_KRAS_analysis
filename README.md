# Mutation_KRAS_analysis

Spanning-tree thermostatistics (STT) of the KRAS G12D mutation: global
thermodynamics, exact heat-capacity decomposition, and allosteric channel
analysis of the active (GPPNHP-bound) state.

This repository contains the code and Jupyter notebooks used to produce the
tables and figures of the accompanying manuscript.

---

## Overview

The protein graph is built from its Cα contact network. Assigning each contact a
Boltzmann weight `w_e(β) = exp(-β d_e)` turns the weighted Laplacian into a
partition function over spanning trees (Kirchhoff's theorem), from which global
thermodynamic quantities follow exactly:

```
ln Z   = ln det L*(β)          (grounded Laplacian cofactor)
⟨E⟩    = -d lnZ / dβ  = -Tr(L⁻¹ L_β)
Var(E) =  d²lnZ / dβ²
C      =  Var(E) / (kT)²
F      = -kT lnZ
S      = (⟨E⟩ - F) / kT
```

Heat capacity is evaluated from **exact β-derivative trace identities**; no
finite-difference step enters the reported result. A finite-difference audit is
included only to show where naive numerical differentiation suffers from
catastrophic cancellation.

## System studied

| Role | PDB | Description |
|------|-----|-------------|
| Wild type | **6GOD** | KRAS wild type, GPPNHP (active) |
| Mutant    | **6GOF** | KRAS G12D, GPPNHP (active) |

Chain **A** is used for both structures.

## Key parameters

- Cα contact cutoff: **7.8 Å** (used consistently across all tables and figures)
- Reduced temperature: **kT = 1.0**
- Grounded node: common index for both structures (grounding-invariant)

---

## Repository contents

| File | Produces | Description |
|------|----------|-------------|
| `Table_1_Global_results.ipynb` | Table 1 | Global thermodynamics (F, ⟨E⟩, S, C) for WT vs G12D — exact trace-identity result with finite-difference, edge-order, and grounding-invariance audits |
| `Table_2.ipynb` | Table 2 | Channel (allosteric corridor) thermodynamics with exact heat-capacity decomposition |
| `Table_3.ipynb` | Table 3 | Residue occupancy / path-probability analysis |
| `Table_3.ipynb` | Figure 2 | Channel convergence: `r(L) = R_ab(L)/R_ab(full)` vs max path length L, for all 10 channels (WT vs G12D) |
| `6GOD.pdb`, `6GOF.pdb` | — | Input structures |

## Allosteric channels

Figure 2 validates the truncated path ensemble for **10 channels** by comparing the truncated
effective resistance `R_ab(L)` against the exact full-graph value `R_ab(full)`.

**Convergence criterion at L = 9:** `r < 1.01` → PASS · `recovery ≥ 93%` →
MARGINAL (report exact %) · else FAIL.

| Channel | Group | Convergence @ L=9 (recovery) |
|---------|-------|------------------------------|
| 6 → 11    | original | PASS (99.6%) |
| 55 → 60   | original | PASS (99.4%) |
| 110 → 117 | original | MARGINAL (98.4%) |
| 141 → 146 | original | MARGINAL (98.9%) |
| 19 → 142  | original | MARGINAL (98.9%) |
| 12 → 35   | new      | PASS (99.6%) |
| 12 → 61   | new      | MARGINAL (98.5–98.9%) |
| 12 → 156  | new      | MARGINAL (97.5%) |
| 12 → 170  | new      | **FAIL (87–88%)** — does not converge by L=10; extend L or exclude from occupancy/heat-capacity tables |
| 35 → 61   | new      | PASS (99.5%) |

Recovery values are computed at cutoff 7.8 Å, kT = 1.0. Any occupancy or
heat-capacity number derived from a non-converged channel should not be
trusted at face value.

## Reproducing the results

### Requirements

```bash
pip install numpy networkx biopython matplotlib
```

Google Colab:

```python
!pip install -q numpy networkx biopython matplotlib
```

### Run

1. Place `6GOD.pdb` and `6GOF.pdb` in the working directory.
2. Open the notebook of interest and run all cells.

Each notebook prints its result and the associated numerical audits. `Figure_2`
also writes `channel_convergence.txt` (the per-channel r(L) table with verdicts)
and `channel_convergence.png`.

> Note: the convergence path enumeration up to L = 10 on the dense 7.8 Å graph
> is computationally heavy and takes time.

To explore cutoff sensitivity, change the `CUTOFF` value at the top of a
notebook (default `7.8`).

---

## Notes on numerical reproducibility

- Coordinates and matrices use `float64`.
- Global-thermodynamics `ln det` uses `numpy.linalg.slogdet` (no eigenvalue
  thresholding); channel resistances use the Moore–Penrose pseudoinverse.
- The **exact** heat capacity is the reported value; finite-difference columns
  are diagnostic only.

---

## Citation

If you use this code, please cite:

> Şengüler Çiftçi, F. Erman, B. *Spanning-Tree Thermostatistics of Protein Allostery: An Exact Kirchhoff Framework with Application to Oncogenic KRAS.*
> Bioarxiv
> DOI:  https://doi.org/10.64898/2026.04.29.721570.

## Contact

Maintained by [@fatmasenguler](https://github.com/fatmasenguler).
