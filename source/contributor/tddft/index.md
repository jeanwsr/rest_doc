# TDDFT Module

The TDDFT module of REST (`src/ri_tddft/`) implements RI-accelerated linear-response time-dependent density functional theory. It reuses the generic Davidson solver (`solvers::davidson`) and the FEAST contour-integration solver (`solvers::feast`), together with the response-equation solvers and dipole tools from the BSE module (`ri_bse`), but builds the diagonal from KS orbital energies (rather than GW quasiparticle energies).

Capability overview:

- **Restricted (RHF/RKS) references**: `tddft_spin` selects the singlet (`"singlet"`), the triplet (`"triplet"`, AO mode only), or both (`"both"`, AO mode only; used for the PySOC export);
- **Unrestricted (UHF/UKS) references (UTDDFT)**: `spin_polarization = true`; the response has a single spin-coupled channel (the α/β excitation sectors are coupled through the spin-independent Coulomb kernel) and `tddft_spin` does not apply;
- **Range-separated hybrid (RSH) functionals**: the response exchange is decomposed as `c_LR·K_full + (c_SR − c_LR)·K_SR` (see "Range-separated hybrids" below);
- **SCF wave-function stability analysis**: the `stability`/`check_stab` keywords, built on the AO-mode (A+B) orbital Hessian (see [tddft-stability](tddft-stability.md));
- **Analytic excited-state gradients**: the `tddft_grad_state` keyword, an RI port of PySCF `grad/tdrks.py` (see [tddft-grad](tddft-grad.md)).

The module supports two kernel implementations, selected by the input keyword `tddft_mode`:

- **MO mode** (`tddft_mode = "mo"`, default): the three-center RI integrals are pre-transformed into the MO basis, and all matrix-vector products run entirely in the MO amplitude space (see [tddft-mo](tddft-mo.md));
- **AO mode** (`tddft_mode = "ao"`): the Davidson iteration still lives in the MO amplitude space, but each matrix-vector product builds AO transition densities; the Coulomb and exchange terms are evaluated through the density-driven interface of `ri_jk`, and the XC kernel is evaluated in batches on the numerical grid via `dft::numint_matmul` (see [tddft-ao](tddft-ao.md)).

Both modes share the same data structure `TDDFTData` (`src/ri_tddft/tddft.rs`): mode-specific members are stored as `Option`s; for an unrestricted reference the per-spin-sector members (`c_occ`/`c_vir`/`psi_occ`/`ri_terms` etc.) are stored as a `Vec` with one entry per sector, the number of sectors being given by `TDDFTData::n_sectors()` (1 for RHF, 2 for UHF). The exchange mixing coefficients are not stored in the data: they are derived per matrix-vector product from `scf.mol.xc_data`.

```{toctree}
---
maxdepth: 1
---

tddft-mo
tddft-ao
tddft-stability
tddft-grad
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
- Singlet (restricted reference): $\kappa_c = 2$
- Triplet (restricted reference): $\kappa_c = 0$

For an unrestricted reference the Coulomb term couples the $\alpha$/$\beta$ spin sectors with unit weight ($\kappa_c = 1$; see "Unrestricted references" below), so the restricted spin-adapted $\kappa_c$ factors above no longer apply. REST's matrix-vector product interface labels the channel by `xlet: char`: `'S'`/`'T'` correspond to the restricted singlet/triplet channels, while `'R'` is the generic non-spin-adapted marker used by the unrestricted and UHF-stability code paths.

### Spin channels of the XC kernel

The XC kernel depends on the spin channel (CPL, 256, 454):

$$
f_s = f_{\uparrow\uparrow} + f_{\uparrow\downarrow} \qquad
f_t = f_{\uparrow\uparrow} - f_{\uparrow\downarrow}
$$

where $f_{\uparrow\uparrow}$ and $f_{\uparrow\downarrow}$ are the spin-resolved second functional derivatives. An unpolarized libxc evaluation only yields $f_u = \frac{1}{2}(f_{\uparrow\uparrow} + f_{\uparrow\downarrow})$, therefore:

- Singlet kernel $f_s = 2 f_u$: obtained directly from an unpolarized evaluation;
- Triplet kernel $f_t$: unobtainable from an unpolarized evaluation — a spin-polarized evaluation at $\rho_\uparrow = \rho_\downarrow = \rho/2$ is required, combined along the antisymmetric direction (the PySCF `nr_rks_fxc_st` recipe). REST currently supports triplets only in **AO mode**; requesting triplets in MO mode raises an explicit error.
- `tddft_spin = "both"`: the eigenvalue problem is solved twice, singlet first and then triplet (the `run_spin` closure inside `tddft_main` prepares the `TDDFTData` per channel — the two channels **do not share** the spin-adapted kernel); both the output order and the JSON field order are singlet-then-triplet.

### Unrestricted references (UTDDFT)

The unrestricted (UHF/UKS) response has a **single** spin-coupled channel: the excitation space is the concatenation $[z_\alpha; z_\beta]$ of the α and β occupied→virtual rotation sectors, coupled through the spin-independent Coulomb kernel ($J[\sum_\tau z^\tau]$). There is no spin-adapted "factor 2 / factor 0" pair, so `tddft_spin` does not apply (an explicit value is rejected). The XC kernel is the spin-resolved kernel $f_{\sigma_1\sigma_2}[g,\alpha,\beta]$ (no singlet/triplet factors): stored as `fxc_u` (spin-resolved `FXCMatvecDataUnrestricted`) in MO mode, or as the spin-polarized `fxc_eff: [n_\mathrm{grid}, n_\mathrm{var}, 2, n_\mathrm{var}, 2]` in AO mode. The orbital windows (`tddft_cutoff_energy`, frozen core) are resolved independently on each spin channel (`tddft_occupation_parameters_u`); an empty sector (e.g. no β occupation) participates as a zero-dimensional sector. Amplitude post-processing (transition dipoles, oscillator strengths, leading-transition printing) follows the PySCF `uhf.py` conventions.

Both MO and AO modes support unrestricted references. The solver layering matches the restricted case, with separate (relaxed) dense thresholds: MO-U TDA $\dim \le 15$ (`dsyev`), MO-U Full LR $\dim \le 80$ (the explicit non-Hermitian $[\mathbf{A}\ \mathbf{B};-\mathbf{B}\ -\mathbf{A}]$ diagonalized via `dgeev`).

### Range-separated hybrid (RSH) functionals

For range-separated hybrids, the HF exchange operator is partitioned as $1/r_{12} = \mathrm{erfc}(\omega r_{12})/r_{12} + \mathrm{erf}(\omega r_{12})/r_{12}$ into short-range (SR) and long-range (LR) parts. REST's response-exchange convention (mirroring the ground-state Fock build in `scf_io`) is

$$
K^{\text{resp}} = c_{LR}\, K_{\text{full}} + (c_{SR} - c_{LR})\, K_{\text{SR}}
$$

where $K_{\text{full}}$ is the ordinary $1/r_{12}$ exchange, $K_{\text{SR}}$ the $\mathrm{erfc}(\omega r_{12})/r_{12}$ exchange, and $c_{LR}$, $c_{SR}$ the RSH parameters of the functional (`xc_data.rsh_params()` → $(\omega, c_{LR}, c_{SR})$). The coefficients are **not pre-stored** in `TDDFTData`: they are derived per matrix-vector product from `scf.mol.xc_data`. The SR three-center integral tensors are built only when $|c_{SR} - c_{LR}| > 10^{-12}$ (MO mode: the `oo_sr`/`vv_sr`/`ov_sr` triple of `RITensorTerms`; AO mode: `scf.rimatr_sr` must exist, i.e. it was built by the SCF with the same functional). For ordinary hybrids $c_{SR} - c_{LR} = 0$ and the formula reduces to $c_x K_{\text{full}}$.

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
| `tddft_method` | String | `"lr"` | `"tda"` (Tamm-Dancoff approximation) or `"lr"` (full linear response) |
| `tddft_mode` | String | `"mo"` | `"mo"` (MO-basis RI tensors) or `"ao"` (AO transition-density kernel) |
| `tddft_spin` | String | `"singlet"` | `"singlet"` / `"triplet"` / `"both"`; triplet and both are AO-mode-only and restricted-reference-only (explicit values are rejected for unrestricted references) |
| `nroots` | Integer | 6 | Number of excited states to solve for |
| `davidson_tol` | Float | `1e-10` | Davidson convergence threshold: $\|r\| < \sqrt{\varepsilon}$ and $|\Delta E| < \varepsilon$ |
| `davidson_max_iter` | Integer | 50 | Maximum Davidson iterations |
| `davidson_max_subspace` | Integer | 60 | Maximum subspace dimension multiplier (full LR needs a large subspace to converge before the first restart; do not shrink it casually) |
| `grid_batch` | Bool | `true` | AO mode only: evaluate the XC kernel in grid batches to keep the full AO-on-grid tensor out of memory; ignored in MO mode |
| `tddft_ao_rik_driver` | String | `"semitrans"` | AO mode only: exchange-K driver — `"semitrans"` (occupied-side semi-transformation, default), `"dm"` (exact batched density-driven), or `"lowrank"` (per-vector SVD) |
| `tddft_fxc_driver` | String | `"semitrans"` | AO mode only: fxc driver — `"semitrans"` (C_vir folded into the amplitudes; the vir side contracts against the raw AO on grid, so no psi_vir is ever formed) or `"mo"` (cached occ-side grid projections with a streamed vir side, the MO-mode fxc algorithm) or `"dm"` (assembled-density NIMatmul fallback; unknown values warn and fall back here) |
| `tddft_svd_tol` | Float | `1e-6` | AO mode only: relative singular-value threshold for low-rank K (keep $\sigma_i \ge \varepsilon \sigma_{\max}$) |
| `stability` | String | `"off"` | SCF stability analysis (`"internal"` / `"external"` / `"full"` / `"auto"`); any value other than `"off"` makes the task stability-only (see [tddft-stability](tddft-stability.md)) |
| `stability_nroots` | Integer | 3 | Number of lowest Hessian eigenvalues to solve for |
| `stability_tol` | Float | `1e-8` | Davidson convergence threshold for the stability Hessian |
| `tddft_grad_state` | Integer | 0 | Analytic excited-state gradient (1-based, 0 disables; restricted references only, see [tddft-grad](tddft-grad.md)) |
| `pysoc` | Bool | `false` | Export the PySOC JSON (`rest_pysoc_export.json`); requires a restricted reference with `tddft_spin = "both"` |
| `tddft_cutoff_energy` | Float | `1e6` | Virtual orbital energy cutoff (Hartree); resolved independently per spin channel for unrestricted references |
| `tddft_feast_solver` | Bool | `false` | FEAST contour-integration solver (MO mode and restricted references only) |
| `response_tddft` | Bool | `false` | Enable response TDDFT (MO implementation, restricted references) |
| `response_tddft_solver` | String | `"klopper"` | `"pople"` / `"gmres"` / `"klopper"` / `"dense"` |

NOTE: the top-level `[ctrl]` keyword `check_stab` accepts the same values as `[tddft] stability`; when both are given the latter wins.

## Workflow

### Eigenvalue TDDFT (`tddft_main`)

```
tddft_main(scf)
    │
    ├── Step 1: read TDDFT control parameters + reference-type gating
    │   ├── tddft_method (tda/lr), tddft_spin, tddft_mode, nroots, ...
    │   ├── ROHF reference → error (use spin_polarization = true instead)
    │   ├── unrestricted reference → tddft_spin must be absent; pysoc needs restricted + "both"
    │   └── restricted reference → resolve tddft_spin (singlet/triplet/both) → xlet
    │       triplet/both in MO mode → error (AO mode only)
    │
    ├── Step 2: resolve the orbital sectors
    │   ├── restricted: tddft_occupation_parameters() → (start_mo, occ_size, vir_size, dim)
    │   └── unrestricted: tddft_occupation_parameters_u() → [α sector, β sector]
    │       frozen core (< -2.0 Ha) and virtual truncation (tddft_cutoff_energy) handled here
    │
    ├── Step 3: initial guess + diagonal preconditioner
    │   └── build_hdiag() + generate_initial_guess() (from solvers/davidson)
    │
    ├── Step 4: solve per spin channel (run_spin closure; data prepared per channel)
    │   ├── data preparation (per channel, spin-adapted kernels not shared):
    │   │   ├── MO mode: prepare_mo_data() → fxc/fxc_u table + per-sector MO-basis RI bundles (ri_terms)
    │   │   └── AO mode: prepare_ao_data_with_spin(scf, Some(spin)) → c_occ/c_vir + NIMatmul + raw kernel table
    │   │
    │   └── layered solver dispatch:
    │       ├── FEAST (restricted + MO only): feast_solve_tddft_tda/lr()
    │       ├── dense diagonalization (small systems):
    │       │   ├── restricted dim ≤ 15: build_a/build_b → dsyev (LR via the (A−B) symmetrized reduction, TDA fallback)
    │       │   ├── MO-U TDA dim ≤ 15: dsyev
    │       │   └── MO-U LR dim ≤ 80: explicit [A B; -B -A] + dgeev
    │       └── iterative Davidson:
    │           ├── AO mode: tda/lr_davidson_solver_batched (batched closures)
    │           ├── MO restricted: tda/lr_davidson_solver (per-vector closures)
    │           └── MO unrestricted: per-vector Davidson over the concatenated [α;β] amplitudes
    │
    ├── Step 5: print excitation energies and oscillator strengths per channel
    │   ├── transition_dipole_square(), normalize() (the *_u variants for unrestricted, PySCF uhf.py conventions)
    │   └── tddft_spin = "both": singlet first, then triplet; pysoc = true exports rest_pysoc_export.json
    │
    └── collected into the "tddft" field of rest_results.json (energies, oscillator_strength)
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

### SCF stability analysis (`stability::stability`)

```
main_driver (after SCF convergence, before the excited-state run)
    │
    ├── resolve the stability mode: [tddft] stability wins, [ctrl] check_stab as fallback
    │
    ├── stability::stability(scf, mode)
    │   ├── internal: RHF/RKS → H = 4(A^S+B^S); UHF/UKS → H = 2(A+B) (concatenated [α;β])
    │   ├── external (RHF/RKS only): H = (A^T+B^T) (the RHF→UHF direction)
    │   └── batched Davidson for the lowest nroots eigenvalues → unstable if λ_min < -1e-5
    │
    └── print the results + the "stability" field of rest_results.json (internal/external roots)
```

See [tddft-stability](tddft-stability.md) for details.

### Analytic excited-state gradient (`tddft_grad::TddftGradEngine`)

When `tddft_grad_state > 0`, gradient tasks add the excited-state response gradient onto the ground-state gradient inside `main_driver::eval_force`. Module details are yet to be written; see [tddft-grad](tddft-grad.md).

## Code structure

The `src/ri_tddft/` directory contains:

| File | Responsibility |
|------|------|
| `mod.rs` | submodule declarations and public API exports (`tddft_main`, `response_tddft`, `TddftGradEngine`) |
| `tddft.rs` | shared data structure `TDDFTData`, the `TDDFTMode`/`FxcDriver` enums, the two builders (`prepare_mo_data`, `prepare_ao_data_with_spin`), plus the mode-dispatching dense small-system builders `build_a`/`build_b` |
| `tddft_solver.rs` | eigenvalue TDDFT driver (`tddft_main`): reference-type gating, sector resolution, per-spin-channel data preparation and solver dispatch (FEAST/dense/Davidson), output printing and JSON summary |
| `matvec.rs` | MO-mode matrix-vector products: the A and B blocks (`a_matvec`, `b_matvec`, with the RSH two-coefficient exchange), per-sector `RITensorTerms` |
| `matvec_ao.rs` | AO-mode matrix-vector products: transition-density construction, batched J/K via `ri_jk`, batched fxc via `numint_matmul`, batched and dense paths |
| `response.rs` | response TDDFT solvers (`response_tddft`): Pople, GMRES, Klopper, dense LU backends, and polarization-density export |
| `stability.rs` | SCF wave-function stability analysis (`stability`, `StabilityReport`): the AO-mode (A+B) orbital Hessian + batched Davidson |
| `tddft_grad.rs` | analytic excited-state gradient (`TddftGradEngine`): an RI port of PySCF `grad/tdrks.py`/`tdrhf.py` |
| `utils.rs` | utilities: `tddft_occupation_parameters`/`_u`/`tddft_sector_params` (orbital sectors), `tddft_get_submatrix` (RI submatrix extraction), `compute_tddft_dipole_matrix` (transition dipoles) |
| `feast_solver.rs` | FEAST solver wrapper (MO mode, restricted references): adapts the TDDFT matrix-vector products to the generic `solvers::feast` interface |

Input parameters are defined in `src/ctrl_io/tddft_parameters.rs`.

## Integration points

### Upstream dependencies

| Module | Purpose |
|----------|------|
| `scf_io::SCF` | core data: MO coefficients, KS orbital energies, RI integrals (`rimatr`, plus `rimatr_sr` for RSH), numerical grids, molecular information |
| `dft::num_int` | `FXCMatvecData`/`FXCMatvecDataUnrestricted` (MO-mode XC kernel data) and `fxc_matvec()` (XC kernel matrix-vector product) |
| `dft::numint_matmul` | AO-mode XC kernel: `NIMatmul` (grid AO cache, batched density construction and kernel contraction), `eval_vxc_fxc_from_rho` (raw kernel table) |
| `dft::xceff` | libxc evaluation wrappers (`libxc_eval_eff`, `determine_den_type`); the triplet spin-polarized kernel is composed on top of these |
| `ri_jk` | AO-mode J/K: `get_vj_ri_incore_nonsym`, `get_vk_ri_incore_dm`, `get_vk_ri_incore_dm_lowrank`, `get_vk_ri_incore_coeff_pair` (notation: the [ri documentation](../ri.md)) |
| `ri_bse` | Coulomb contribution (`coulomb_contribution`), response-equation solvers (Pople/GMRES/Klopper, adapted from `ri_bse::response`), dipole tools (`dipoles::normalize` etc.), `pysoc_export` |
| `solvers::davidson` | generic Davidson solvers: per-vector interface (`davidson_solver`, `lr_davidson_solver`) and batched interface (`davidson_solver_batched`, `lr_davidson_solver_batched`) |
| `solvers::feast` | generic FEAST contour-integration algorithm |
| `ri_cphf` | CPHF Z-vector solve for the gradient module (`CPHFSolverPySCF`) |
| `ri_gw::gw_grad` | raw RI tensors and atom-derivative blocks reused by the gradient module (`RawRiTensors`, `AtomDerivBlocks`) |

### Relationship with BSE

The TDDFT module reuses infrastructure from `ri_bse` (response solvers, dipole tools), with two key differences:

1. **Diagonal**: TDDFT uses KS orbital-energy differences $\varepsilon_a - \varepsilon_i$, whereas BSE uses GW quasiparticle-energy differences.
2. **Data flow**: TDDFT sits directly in the `main_driver` call chain (alongside PT2 and RPA), whereas BSE requires a preceding GW calculation.
