# SCF Wave-Function Stability Analysis

An SCF solution may be a **saddle point** of the energy surface rather than a minimum (the classic example being the unstable RHF solution for $\mathrm{O_2}$, $\mathrm{O_3}$, etc.). REST provides a wave-function stability analysis based on the (A+B) orbital Hessian (Seeger & Pople, JCP 66, 3045 (1977)), executed after SCF convergence: the solution is unstable when the lowest eigenvalue $\lambda_{\min} < -10^{-5}$, and the eigenvector of a negative eigenvalue spans the downhill orbital-rotation direction.

The checks supported by REST:

- **internal stability**: the RHF/RKS singlet channel or the UHF/UKS orbital Hessian;
- **external stability**: the RHF→UHF direction (triplet channel), restricted references only.

The analysis is **check-only**: it reports instabilities and their directions, but does not automatically rotate the orbitals and re-run the SCF.

## Input keywords

The following keywords are declared in the `[tddft]` block of the input deck:

| Keyword | Type | Default | Description |
|---|---|---|---|
| `stability` | String | `"off"` | Stability mode: `"internal"` (internal only), `"external"` (external only, RHF/RKS only), `"full"` (internal + external), `"auto"` (recommended; every check applicable to the reference). Any value other than `"off"` makes this task **stability-only — no excited-state calculation is performed** |
| `stability_nroots` | usize | 3 | Number of lowest stability-Hessian eigenvalues to solve for |
| `stability_tol` | f64 | 1.0e-8 | Davidson convergence threshold for the stability Hessian |

NOTE: the `check_stab` keyword in the top-level `[ctrl]` block accepts the same values as `stability` and serves as its fallback; when both are present, `[tddft] stability` wins. A stability-only deck needs no other excited-state settings (`stability_nroots`/`stability_tol` have defaults; the `[tddft]` block may contain the `stability` line alone).

## Example

```toml
[ctrl]
     xc =                        "b3lyp"
     basis_path =                "def2-tzvp"
     auxbas_path =               "def2-tzvp-rifit"
     eri_type =                  "ri-v"

[tddft]
     stability =                 "auto"
     stability_nroots =          3
```

The program runs the stability analysis after SCF convergence, printing the lowest eigenvalues and the verdict per check, e.g.:

```sh
RHF/RKS internal stability: lowest eigenvalues = [0.012, 0.045, 0.132]
RHF/RKS wavefunction is stable in the internal stability analysis
```

The results are also written to the `"stability"` field of `rest_results.json` (roots and stable flags for `internal`/`external`) for scripted post-processing and regression testing.

## Interpreting the results

| Reference | internal | external |
|--|--|--|
| RHF/RKS | $4(\mathbf{A}^{\mathrm{S}} + \mathbf{B}^{\mathrm{S}})$, singlet instability directions | $\mathbf{A}^{\mathrm{T}} + \mathbf{B}^{\mathrm{T}}$, RHF→UHF instability direction |
| UHF/UKS | $2(\mathbf{A} + \mathbf{B})$ (concatenated α/β rotation space) | not implemented (UHF→GHF), automatically skipped |

Criterion: the channel is stable iff $\lambda_{\min} \ge -10^{-5}$; a negative eigenvalue means the SCF solution is a saddle point in that channel — consider a different initial guess (`initial_guess`), mixer, or an unrestricted reference.

## Notes

- **Mutually exclusive with the excited-state run**: when `stability` is not `"off"`, the task only performs the stability check and does not compute excitation energies.
- **Numerical grids**: DFT references require the grids (prepared automatically by a regular SCF); HF references need no grids and automatically evaluate the RI J/K parts only.
- **Reference types**: ROHF references are not supported; the real→complex check and the UHF→GHF external check are not yet implemented.
- `"auto"` is recommended: it runs every check applicable to the current reference.

Developer details (theory factors, Davidson configuration, implementation flow): see [SCF wave-function stability analysis (contributor)](../contributor/tddft/tddft-stability.md).
