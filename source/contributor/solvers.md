# The `solvers` Module

> REST's unified iterative solvers, located in `src/solvers/`, providing a Davidson eigensolver and a Krylov linear-system solver. The module was introduced in [PR#183](https://gitee.com/restgroup/rest/pulls/183).

## Module layout

```
src/solvers/
├── mod.rs          # module declarations
├── davidson.rs     # Davidson solvers (davidson_solver / lr_davidson_solver)
├── krylov.rs       # Krylov solvers (krylov_tsr / krylov_vec)
└── unify.md        # solver design doc (architecture, features, test coverage)
```

---

## Davidson solver

### API

| Function | Purpose |
|---|---|
| `davidson_solver()` | TDA eigenvalue problem `A x = e x` (alias `tda_davidson_solver`) |
| `lr_davidson_solver()` | Symmetrized Casida equation `(A-B)^(1/2)(A+B)(A-B)^(1/2) z = omega^2 z` |
| `generate_initial_guess()` | builds initial trial vectors from the smallest elements of `diag` |

Both solvers take an `FnMut(&Vec<f64>) -> Vec<f64>` matvec closure and return `Vec<(f64, Vec<f64>)>` (eigenvalue, eigenvector).

### Convergence criterion

$$
\|\mathbf{r}\| < \sqrt{\mathrm{tol}} \quad \wedge \quad |\Delta e| < \mathrm{tol}
$$

where $\mathbf{r} = A\mathbf{x} - e\mathbf{x}$ is the true residual and $\Delta e$ is the eigenvalue change between successive iterations.

### `DavidsonConfig` defaults

| Field | Default | Description |
|---|---|---|
| `max_subspace` | `60` | maximum subspace dimension |
| `add_dim` | `2` | vectors added per iteration |
| `restart_dim` | `6` | subspace size kept after restart |
| `max_iter` | `100` | maximum number of iterations |
| `tol` | `1e-10` | convergence threshold (see criterion above) |
| `lindep` | `1e-14` | linear-dependence threshold; vectors with `\|\|v\|\|^2 < lindep` are dropped |
| `use_mgs` | `true` | Modified Gram-Schmidt orthogonalization (`false` uses CGS) |
| `divergence_restart` | `true` | detect divergence (`\|\|r\|\|/\|\|r_last\|\| > 3`), restore and restart |
| `track_states` | `false` | reorder eigenstates between iterations via overlap matrices |

### Features

- Configurable convergence threshold `tol`
- Preconditioner denominator clamp: `|denom|` floored at `1e-8` to avoid blow-up
- MGS / CGS orthogonalization toggle
- Divergence detection with automatic restart
- Eigenstate tracking for near-degenerate roots
- Linear-dependence filtering of trial vectors (`lindep`)

### Logging

- Output verbosity is controlled by the `log` framework, mapped from the `print_level` input keyword
- Each iteration prints two lines: `residues: [...]` / `|de|: [...]`, plus an `iter N: space=, converged=` summary

---

## Krylov solver

### API

| Function | Purpose |
|---|---|
| `krylov_tsr()` | Tsr/TsrView core with the full algorithm |
| `krylov_vec()` | `Vec<f64>` convenience wrapper (`krylov` is an alias) |

Solves `(I + A) x = b`. Returns `(Tsr, KrylovResult)`, where `KrylovResult` records `cycles`, `aop_calls`, `per_root_solves`, and `residual`.

### Convergence and recovery

The subspace-staleness signal is:

$$
\max(\|\mathbf{v}\|^2) < \max(\mathrm{lindep},\ \mathrm{tol}^2)
$$

When the subspace collapses (QR finds no new directions) but the true residual remains large, a per-RHS recovery mechanism kicks in:

1. **Per-RHS independent subspaces**: each RHS maintains its own Krylov basis (matvec remains batched)
2. **Recursive per-root solves**: RHS that still fail to converge are solved recursively via independent single-RHS `krylov_tsr` calls (fresh subspace, no shared coupling)
3. **Tolerance factor**: RHS with true residual `||r|| < max_residual_factor * tol` are accepted directly without a recursive solve

### `KrylovConfig` defaults

| Field | Default | Description |
|---|---|---|
| `tol` | `1e-9` | convergence tolerance |
| `max_cycle` | `50` | maximum number of cycles |
| `max_space` | `None` | hard-restart subspace cap (`None` = no restart) |
| `lindep` | `1e-15` | linear-dependence threshold |
| `max_residual_factor` | `1000.0` | tolerance factor; `\|\|r\|\| < factor * tol` is accepted |

---

## Callers and user-keyword mapping

| Solver | Caller | User keyword | Default |
|---|---|---|---|
| Davidson | BSE | `davidson_converge_threshold` | `1e-10` |
| Davidson | TD-DFT | `davidson_tol` | `1e-10` |
| Davidson | SCF stability (`addons.rs`) | — | `tol = 1e-5` |
| Krylov | `analdrv` | `cphf_tol` | `1e-9` |
| Krylov | `analdrv` | `cphf_lindep` | `1e-15` |
| Krylov | `analdrv` | `cphf_tol_inflation` | `1000.0` |
| Krylov | Hessian (`[ctrl.hessian]`) | `krylov_tol` | `1e-9` |
| Krylov | Hessian | `krylov_lindep` | `1e-15` |
| Krylov | Hessian | `krylov_tol_inflation` | `1000.0` |

> Note: Davidson's `tol` controls both convergence criteria — the true residual `||r|| < sqrt(tol)` and the eigenvalue change `|de| < tol`. Krylov's `tol` instead controls the subspace-staleness signal `max(||v||^2) < max(lindep, tol^2)`. The two semantics differ and should not be confused.

### Entry points

- BSE: `src/ri_bse/mod.rs` (four `davidson_solver` / `lr_davidson_solver` call sites)
- TD-DFT: `src/ri_tddft/tddft_solver.rs`
- SCF stability: `src/scf_io/addons.rs`
- CP-HF (Hessian): `src/ri_cphf/cphf_solver_pyscf.rs` → `src/hessian/rhf.rs`
- analdrv: `src/analdrv/krylov_block.rs` (thin wrapper around `krylov_tsr`)
