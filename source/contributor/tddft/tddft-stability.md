# SCF wave-function stability analysis

Function path: `ri_tddft::stability::stability` (`src/ri_tddft/stability.rs`)

## Overview

An SCF solution may be a **saddle point** of the energy surface rather than a minimum (the classic example being the unstable RHF solution for $\mathrm{O_2}$, $\mathrm{O_3}$, etc.). The stability analysis expands the energy to second order around the converged solution and checks whether the Hessian in the orbital-rotation space is positive definite: a negative eigenvalue signals a downhill orbital-rotation direction, i.e. an unstable SCF solution.

REST's stability analysis is built entirely on the AO-mode TDDFT (A/B) matrix-vector machinery (see [tddft-ao](tddft-ao.md)): the stability operator is the **(A+B) orbital Hessian**, and its lowest eigenvalues are obtained with the batched Davidson solver. The analysis is **check-only**: it reports instabilities and their directions (the eigenvectors), but does not automatically rotate the orbitals and re-run the SCF (orbital following is not implemented).

## Theoretical background

The stability operator is the second derivative of the energy with respect to occupied→virtual orbital rotations (Seeger & Pople, JCP 66, 3045 (1977)), including the B-type (de-excitation / symmetrized-density) response, expressed through the TDDFT A and B operators on the same reference:

$$
\mathbf{H}_{\mathrm{stab}} = f \cdot (\mathbf{A}^{x} + \mathbf{B}^{x})
$$

The channels and factors used by REST (matching PySCF `rhf_internal`/`uhf_internal`/`rhf_external`):

| Check | Hessian | Channel `$x$` | Factor `$f$` |
|--|--|--|--|
| RHF/RKS internal | $4(\mathbf{A}^{\mathrm{S}} + \mathbf{B}^{\mathrm{S}})$ | `'S'` (singlet) | 4 |
| RHF/RKS external (RHF→UHF) | $\mathbf{A}^{\mathrm{T}} + \mathbf{B}^{\mathrm{T}}$ | `'T'` (triplet) | 1 |
| UHF/UKS internal | $2(\mathbf{A} + \mathbf{B})$ (concatenated $[\alpha;\beta]$ rotation space) | `'R'` | 2 |

**Stability criterion**: stable iff $\lambda_{\min} \ge -10^{-5}$ (`STABILITY_THRESHOLD`, the PySCF threshold); otherwise an instability is reported, and the eigenvector of a negative eigenvalue spans the downhill orbital-rotation direction.

References: Seeger & Pople, JCP 66, 3045 (1977); Bauernschmitt & Ahlrichs, JCP 104, 9047 (1996); PySCF `scf/stability.py`.

## Usage

```toml
[tddft]
     stability =                 "auto"    # off | internal | external | full | auto
     stability_nroots =          3         # number of lowest eigenvalues to solve for
     stability_tol =             1.0e-8    # Davidson convergence threshold
```

- Values of `stability`:

| Value | Behavior |
|--|--|
| `"off"` (default) | no stability analysis |
| `"internal"` | RHF/RKS singlet or UHF/UKS orbital Hessian |
| `"external"` | RHF→UHF (triplet channel) check, RHF/RKS only; for UHF/UKS a note is printed and the check is skipped (UHF→GHF not implemented) |
| `"full"` | internal + external |
| `"auto"` (recommended) | every check applicable to the reference type; currently resolves to `"full"` |

- The top-level `[ctrl]` keyword `check_stab` accepts the same values and serves as the fallback when `stability = "off"`; when both are given, `[tddft] stability` wins. A stability-only deck needs no further excited-state settings (`stability_nroots`/`stability_tol` have struct defaults).
- The stability analysis is **mutually exclusive** with the excited-state run: when `stability != "off"` the task only performs the stability check (executed in `main_driver`, after SCF convergence and before the excited-state calculation).
- Results are printed per check and written to the `"stability"` field of `rest_results.json` (roots and stable flags for `internal`/`external`), machine-readable for regression testing.

## Implementation details

```text
stability(scf, mode)
    ├── resolve nroots/tol (struct defaults when [tddft] absent; "auto" → "full")
    ├── ROHF reference → error; DFT reference without grids → error (HF references need no grids)
    ├── internal:
    │   ├── (factor, xlet) = RHF: (4, 'S') / UHF: (2, 'R')
    │   ├── prepare_ao_data_with_spin(scf, Some("singlet") | None)
    │   └── hessian_roots(): batched Davidson on factor×(A+B)
    ├── external (RHF/RKS only):
    │   ├── prepare_ao_data_with_spin(scf, Some("triplet"))
    │   └── hessian_roots(): (factor, xlet) = (1, 'T')
    └── StabilityReport { stable_internal, roots_internal,
                          stable_external, roots_external }
```

- **Hessian-vector products**: the sum of `a_matvec_ao_batched` and `b_matvec_ao_batched`, scaled by the factor (`xlet` enters the matrix-vector products). The diagonal preconditioner is `factor × (ε_a − ε_i)`.
- **Davidson configuration**: `max_subspace = 60`, `max_iter = 200`, `add_dim = nroots + 2`, `tol = stability_tol` (with REST's convergence convention $\|r\| < \sqrt{\varepsilon}$, $10^{-8}$ corresponds to a residual of $10^{-4}$, matching PySCF's `STAB_TOL`).
- **HF references**: `prepare_ao_data_with_spin` returns `fxc_driver: None` for references without libxc components; the Hessian then contains only the RI J/K parts, so no DFT grid is required.
- **Data prepared per check**: internal (singlet kernel) and external (triplet kernel) each call `prepare_ao_data_with_spin` independently — the spin-adapted kernels are not shared.

## Limitations

- The real→complex ($^1(\mathrm{A}'-\mathrm{B}')$) check is not implemented;
- The UHF→GHF external check is not implemented (UHF/UKS only performs internal);
- ROHF references are not supported;
- Orbital following (rotating the orbitals and re-running the SCF after detecting an instability) is not implemented; the current mode is check-only.
