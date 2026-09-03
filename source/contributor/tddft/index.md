# TDDFT Module

The TDDFT module of REST (`src/ri_tddft/`) implements RI-accelerated linear-response time-dependent density functional theory. It reuses the Davidson solver, response-equation solvers, and the FEAST contour-integration solver from the BSE module (`ri_bse`), but builds the diagonal from KS orbital energies (rather than GW quasiparticle energies).

The module supports two kernel implementations, selected by the input keyword `tddft_mode`:

- **MO mode** (`tddft_mode = "mo"`, default): the three-center RI integrals are pre-transformed into the MO basis, and all matrix-vector products run entirely in the MO amplitude space (see [tddft-mo](tddft-mo.md));
- **AO mode** (`tddft_mode = "ao"`): the Davidson iteration still lives in the MO amplitude space, but each matrix-vector product builds AO transition densities; the Coulomb and exchange terms are evaluated through the density-driven interface of `ri_jk`, and the XC kernel is evaluated in batches on the numerical grid via `dft::numint_matmul` (see [tddft-ao](tddft-ao.md)).

Both modes share the same data structure `TDDFTData` (`src/ri_tddft/tddft.rs`): mode-specific members are stored as `Option`s, and the common scalar (`alpha_hybrid`) is stored directly.

```{toctree}
---
maxdepth: 1
---

tddft-mo
tddft-ao
```

## Theoretical Background

### The Casida equation

Within KS-DFT, the linear-response TDDFT eigenvalue problem can be written as the non-Hermitian Casida equation:

$$
\begin{pmatrix}
\mathbf{A} & \mathbf{B} \\
-\mathbf{B} & -\mathbf{A}
\end{pmatrix}
\begin{pmatrix}
\mathbf{X} \\
\mathbf{Y}
\end{pmatrix}
= \omega
\begin{pmatrix}
\mathbf{X} \\
\mathbf{Y}
\end{pmatrix}
$$

The matrix elements of the $\mathbf{A}$ and $\mathbf{B}$ blocks (in the MO basis) are:

$$
A_{ia,jb} = (\varepsilon_a - \varepsilon_i) \delta_{ij} \delta_{ab} + \kappa_c (ia|jb) - c_x (ij|ab) + f_{ia,jb}^{\text{xc}}
$$

$$
B_{ia,jb} = \kappa_c (ia|bj) - c_x (ib|aj) + f_{ia,jb}^{\text{xc}}
$$

where $\varepsilon_i, \varepsilon_a$ are the KS orbital energies (occupied and virtual), $(ia|jb)$ are the two-electron integrals in the MO basis, $f^{\text{xc}}$ is the exchange-correlation kernel, $c_x$ is the HF exchange mixing coefficient, and $\kappa_c$ is the Coulomb coupling factor:

- Pure functionals: $c_x = 0$
- Hybrid functionals: $c_x = \alpha_{\text{hybrid}}$
- Singlet: $\kappa_c = 2$
- Unpolarized ('R'): $\kappa_c = 1$
- Triplet: $\kappa_c = 0$

### Spin channels of the XC kernel

The XC kernel depends on the spin channel (CPL, 256, 454):

$$
f_s = f_{\uparrow\uparrow} + f_{\uparrow\downarrow} \qquad
f_t = f_{\uparrow\uparrow} - f_{\uparrow\downarrow}
$$

where $f_{\uparrow\uparrow}$ and $f_{\uparrow\downarrow}$ are the spin-resolved second functional derivatives. An unpolarized libxc evaluation only yields $f_u = \frac{1}{2}(f_{\uparrow\uparrow} + f_{\uparrow\downarrow})$, therefore:

- Singlet kernel $f_s = 2 f_u$: obtained directly from an unpolarized evaluation;
- Triplet kernel $f_t$: unobtainable from an unpolarized evaluation — a spin-polarized evaluation at $\rho_\uparrow = \rho_\downarrow = \rho/2$ is required, combined along the antisymmetric direction (the PySCF `nr_rks_fxc_st` recipe). REST currently supports triplets only in **AO mode**; requesting triplets in MO mode raises an explicit error.

### The Tamm-Dancoff approximation (TDA)

Neglecting the $\mathbf{B}$ block (i.e. setting $\mathbf{Y} = \mathbf{0}$) reduces the problem to a standard Hermitian eigenvalue problem:

$$
\mathbf{A} \mathbf{X} = \omega \mathbf{X}
$$

TDA usually gives excitation energies in reasonable agreement with experiment and is more efficient (no non-Hermitian matrix to handle).

### Full linear response (Full LR)

For the full Casida equation, a similarity transform yields a Hermitian form:

$$
(\mathbf{A} - \mathbf{B})^{1/2}(\mathbf{A} + \mathbf{B})(\mathbf{A} - \mathbf{B})^{1/2} \mathbf{Z} = \omega^2 \mathbf{Z}
$$

After solving, the X, Y vectors are recovered via the back-transformations $\mathbf{X} = (\mathbf{A} - \mathbf{B})^{1/2} \mathbf{Z} / \sqrt{\omega}$ and $\mathbf{Y} = \mathbf{X} - \mathbf{Z} / \sqrt{\omega}$.

### Response TDDFT

Frequency-domain response TDDFT solves the following four-component non-Hermitian linear system:

$$
\begin{pmatrix}
\mathbf{A} - \omega & \mathbf{B} & -\gamma & \mathbf{0} \\
\mathbf{B} & \mathbf{A} + \omega & \mathbf{0} & \gamma \\
\gamma & \mathbf{0} & \mathbf{A} - \omega & \mathbf{B} \\
\mathbf{0} & -\gamma & \mathbf{B} & \mathbf{A} + \omega
\end{pmatrix}
\begin{pmatrix}
\mathbf{r}^+ \\
\mathbf{r}^- \\
\mathbf{i}^+ \\
\mathbf{i}^-
\end{pmatrix}
=
\begin{pmatrix}
\boldsymbol{\mu}_z \\
\boldsymbol{\mu}_z \\
\mathbf{0} \\
\mathbf{0}
\end{pmatrix}
$$

where $\omega$ is the field frequency, $\gamma$ the lifetime broadening, and $\boldsymbol{\mu}_z$ the z-component transition dipole vector. The frequency-dependent polarizability $\alpha_{zz}(\omega)$ is computed after solving.

Response TDDFT (`response_tddft`) currently has an MO implementation only and is independent of `tddft_mode` — regardless of the keyword, the response solvers always use the MO machinery.

## Input keywords

The main TDDFT-related input keywords (defined in `src/ctrl_io/tddft_parameters.rs`):

| Keyword | Type | Default | Description |
|---|---|---|---|
| `tddft_method` | String | `"tda"` | `"tda"` (Tamm-Dancoff approximation) or `"lr"` (full linear response) |
| `tddft_mode` | String | `"mo"` | `"mo"` (MO-basis RI tensors) or `"ao"` (AO transition-density kernel) |
| `tddft_spin` | String | `"singlet"` | `"singlet"` or `"triplet"`; triplets currently supported in AO mode only |
| `nroots` | Integer | 1 | Number of excited states to solve for |
| `grid_batch` | Bool | `true` | AO mode only: evaluate the XC kernel in grid batches to keep the full AO-on-grid tensor out of memory; ignored in MO mode |
| `tddft_ao_rik_driver` | String | `"semitrans"` | AO mode only: exchange-K driver — `"semitrans"` (occupied-side semi-transformation, default), `"dm"` (exact batched density-driven), or `"lowrank"` (per-vector SVD) |
| `tddft_fxc_driver` | String | `"semitrans"` | AO mode only: fxc driver — `"semitrans"` (C_vir folded into the amplitudes; the vir side contracts against the raw AO on grid, so no psi_vir is ever formed) or `"mo"` (cached occ-side grid projections with a streamed vir side, the MO-mode fxc algorithm) or `"dm"` (assembled-density NIMatmul fallback) |
| `tddft_svd_tol` | Float | `1e-6` | AO mode only: relative singular-value threshold for low-rank K (keep $\sigma_i \ge \varepsilon \sigma_{\max}$) |
| `tddft_feast_solver` | Bool | `false` | FEAST contour-integration solver (MO mode only) |
| `response_tddft` | Bool | `false` | Enable response TDDFT (MO implementation) |
| `response_tddft_solver` | String | `"klopper"` | `"pople"` / `"gmres"` / `"klopper"` / `"dense"` |

## Workflow

### Eigenvalue TDDFT (`tddft_main`)

```
tddft_main(scf)
    │
    ├── Step 1: read TDDFT control parameters
    │   └── tddft_method (tda/lr), tddft_spin, tddft_mode, nroots, ...
    │
    ├── Step 2: determine orbital dimensions
    │   └── tddft_occupation_parameters() → (start_mo, occ_size, vir_size, dim)
    │       frozen core (< -2.0 Ha) and virtual truncation (tddft_cutoff_energy) handled here
    │
    ├── Step 3: prepare the shared TDDFTData
    │   ├── MO mode: prepare_mo_data() → fxc table + four MO-basis RI tensors
    │   └── AO mode: prepare_ao_data() → c_occ/c_vir + NIMatmul + raw fxc kernel table
    │
    ├── Step 4: build the diagonal preconditioner
    │   └── matvec::build_hdiag() → hdiag = ε_a - ε_i
    │
    ├── Step 5: generate the initial guess
    │   └── generate_initial_guess() (from solvers/davidson)
    │
    ├── Step 6: dispatch the solver
    │   ├── dim ≤ 15: dense diagonalization (build_a/build_b + LAPACK dsyev)
    │   ├── FEAST enabled (MO only): feast_solve_tddft_tda() or feast_solve_tddft_lr()
    │   └── default: iterative Davidson solver (batched interface)
    │       ├── TDA: davidson_solver_batched()
    │       └── Full LR: lr_davidson_solver_batched()
    │
    └── Step 7: print excitation energies and oscillator strengths
        └── transition_dipole_square(), normalize()
```

### Response TDDFT (`response_tddft`)

```
response_tddft(scf)
    │
    ├── prepare MO-mode data (prepare_mo_data, independent of tddft_mode)
    │
    ├── build the KS-energy diagonal vector + dipole vector
    │
    ├── build the 4-component matrix-vector closure
    │   (A, B, fxc, exchange all packed into H_4c * x)
    │
    ├── select the solver
    │   ├── "pople":   Pople-Krylov numerical trick
    │   ├── "gmres":   4-component GMRES (restart = 30)
    │   ├── "klopper": Klopper subspace solver (default)
    │   └── "dense":   dense LU factorization (LAPACK)
    │
    ├── compute the dynamic polarizability α_zz(ω)
    │
    └── (optional) export the polarization density to a spatial-grid file
```

## Code structure

The `src/ri_tddft/` directory contains:

| File | Responsibility |
|------|------|
| `mod.rs` | submodule declarations and public API exports (`tddft_main`, `response_tddft`) |
| `tddft.rs` | shared data structure `TDDFTData` with its two builders (`prepare_mo_data`, `prepare_ao_data`), plus the mode-dispatching dense small-system builders `build_a`/`build_b` |
| `tddft_solver.rs` | eigenvalue TDDFT driver (`tddft_main`): parameter parsing, data preparation, solver dispatch, output printing |
| `matvec.rs` | MO-mode matrix-vector products: the A and B blocks (`a_matvec`, `b_matvec`) |
| `matvec_ao.rs` | AO-mode matrix-vector products: transition-density construction, batched J/K via `ri_jk`, batched fxc via `numint_matmul`, batched and dense paths |
| `response.rs` | response TDDFT solvers (`response_tddft`): Pople, GMRES, Klopper, dense LU backends, and polarization-density export |
| `utils.rs` | utilities: `tddft_occupation_parameters` (orbital dimensions), `tddft_get_submatrix` (RI submatrix extraction), `compute_tddft_dipole_matrix` (transition dipoles) |
| `feast_solver.rs` | FEAST solver wrapper (MO mode only): adapts the TDDFT matrix-vector products to the generic `ri_bse::feast_solver::feast()` interface |

Input parameters are defined in `src/ctrl_io/tddft_parameters.rs`.

## Integration points

### Upstream dependencies

| Module | Purpose |
|----------|------|
| `scf_io::SCF` | core data: MO coefficients, KS orbital energies, RI integrals (`rimatr`), numerical grids, molecular information |
| `dft::num_int` | `FXCMatvecData` (MO-mode XC kernel data) and `fxc_matvec()` (XC kernel matrix-vector product) |
| `dft::numint_matmul` | AO-mode XC kernel: `NIMatmul` (grid AO cache, batched density construction and kernel contraction) |
| `ri_jk` | AO-mode J/K: `get_vj_ri_incore_nonsym`, `get_vk_ri_incore_dm`, `get_vk_ri_incore_dm_lowrank` (notation: the [ri documentation](../ri.md)) |
| `ri_bse` | Coulomb contribution (`coulomb_contribution`), response-equation solvers (Pople/GMRES/Klopper), dipole tools, FEAST algorithm |
| `solvers::davidson` | generic Davidson solvers: per-vector interface (`davidson_solver`, `lr_davidson_solver`) and batched interface (`davidson_solver_batched`, `lr_davidson_solver_batched`) |

### Downstream consumers

| Call site | Purpose |
|----------|------|
| `main_driver.rs` | invokes `tddft_main()` or `response_tddft()` after SCF convergence |
| `ri_cphf/cphf_solver_pyscf.rs` | uses `tddft_occupation_parameters()` for orbital dimensions |
| `dft/num_int.rs` | uses `tddft_occupation_parameters()` for orbital dimensions |
| `dft/response.rs` | uses `tddft_occupation_parameters()` for orbital dimensions |

### Relationship with BSE

The TDDFT module reuses substantial infrastructure from `ri_bse` (solvers, dipole tools, response solvers), with two key differences:

1. **Diagonal**: TDDFT uses KS orbital-energy differences $\varepsilon_a - \varepsilon_i$, whereas BSE uses GW quasiparticle-energy differences.
2. **Data flow**: TDDFT sits directly in the `main_driver` call chain (alongside PT2 and RPA), whereas BSE requires a preceding GW calculation.
