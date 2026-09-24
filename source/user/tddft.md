# TDDFT Excited-State Calculations

Time-dependent density functional theory (TDDFT) is the workhorse for computing molecular excitation energies and properties. REST supports RI (Resolution of Identity)-accelerated linear-response TDDFT, including both the TDA (Tamm-Dancoff Approximation) and the full linear response (Full LR) schemes, as well as frequency-domain response TDDFT.

In addition, REST provides:

- **Unrestricted TDDFT (UTDDFT)**: excited-state calculations on an unrestricted (UHF/UKS) reference with `spin_polarization = true`;
- **TDDFT for range-separated hybrid (RSH) functionals** ($\omega$B97X, CAM-B3LYP, etc.);
- **SCF wave-function stability analysis**: detects whether the SCF solution is a saddle point (RHF/RKS and UHF/UKS references);
- **Analytic excited-state gradients**: the analytic nuclear gradient of a selected excited state (`tddft_grad_state`), usable directly in geometry optimizations and similar tasks;
- **PySOC export**: exports the singlet/triplet amplitudes to a PySOC-readable JSON file for subsequent spin-orbit coupling calculations.

## Calculation modes

REST's TDDFT module provides two calculation modes (the input keyword `tddft_mode`):

1. **Eigenvalue TDDFT**: solves the Casida equation and directly yields excitation energies $\omega_n$ and oscillator strengths. Both TDA and Full LR are supported.
2. **Response TDDFT**: solves the frequency-space linear system to obtain the frequency-dependent dynamic polarizability $\alpha_{zz}(\omega)$. Several iterative solvers are available (MO implementation and restricted references only).

Both modes build the diagonal from the KS orbital energies and accelerate the Coulomb and exchange matrix-vector products through RI. The exchange-correlation (fxc) kernel is evaluated via libxc, supporting LDA, GGA, hybrid, and range-separated hybrid functionals.

### Kernel implementation mode `tddft_mode`

| Value | Description |
|---|---|
| `"mo"` (default) | The three-center RI integrals are pre-transformed into the MO basis; all matrix-vector products run in the MO amplitude space. The mature, extensively validated path |
| `"ao"` | The Davidson iteration stays in the MO amplitude space, but each matrix-vector product builds AO transition densities; the Coulomb/exchange terms go through `ri_jk` and the XC kernel through batched grid evaluation. Significantly lower memory usage; **triplets and `tddft_spin = "both"` are supported in this mode only**, as are HF references for stability analysis |

General advice: use the default `"mo"` mode for routine singlet valence excitations; use `"ao"` for large systems (where the MO-basis RI tensors become a memory bottleneck), for triplet excitations, and for stability analysis on HF references.

### Unrestricted references (UTDDFT)

Setting `spin_polarization = true` in `[ctrl]` (UHF/UKS reference) runs an unrestricted TDDFT. The unrestricted response has a single spin-coupled channel (the α and β excitation sectors are coupled by the spin-independent Coulomb kernel), so the `tddft_spin` keyword does not apply — an explicit `tddft_spin` in an unrestricted run raises an error. Both the MO and AO kernel modes support unrestricted references.

## Input keywords

All TDDFT-related keywords are declared in the `[tddft]` block of the input deck.

### Basic settings

| Keyword | Type | Default | Description |
|---|---|---|---|
| `tddft_method` | String | `"lr"` | TDDFT method. `"tda"` for the Tamm-Dancoff approximation, `"lr"` for full linear response |
| `tddft_mode` | String | `"mo"` | Kernel implementation mode: `"mo"` (MO-basis RI tensors) or `"ao"` (AO transition-density kernel) |
| `tddft_spin` | String | `"singlet"` | Spin channel (restricted references only). `"singlet"` singlet excitations, `"triplet"` triplets (AO mode only), `"both"` singlet + triplet (AO mode only). An explicit value in an unrestricted run (`spin_polarization = true`) raises an error |
| `nroots` | usize | 6 | Number of excited states to compute |
| `tddft_use_optimized_fxc` | bool | true | Use the rayon-parallel fxc matrix-vector kernel |
| `tddft_cutoff_energy` | f64 | 1.0e6 | Virtual orbital energy cutoff (Hartree). Virtuals above this energy are excluded from the excitation space; resolved independently per spin channel for unrestricted references |

### Davidson solver settings

| Keyword | Type | Default | Description |
|---|---|---|---|
| `davidson_tol` | f64 | 1.0e-10 | Davidson convergence threshold: residual norm $\|r\| < \sqrt{\varepsilon}$ and energy change between iterations $|\Delta E| < \varepsilon$ |
| `davidson_max_iter` | usize | 50 | Maximum Davidson iterations |
| `davidson_max_subspace` | usize | 60 | Maximum subspace dimension multiplier (subspace dimension = `nroots × davidson_max_subspace`). Full LR needs a large enough subspace to converge before the first restart; avoid shrinking it |

### AO-mode advanced options

The following keywords take effect only with `tddft_mode = "ao"` (ignored in MO mode):

| Keyword | Type | Default | Description |
|---|---|---|---|
| `grid_batch` | bool | true | Evaluate the XC kernel in grid batches, keeping the full AO-on-grid tensor out of memory. Trades a small time overhead for roughly 40% lower peak memory |
| `tddft_ao_rik_driver` | String | `"semitrans"` | Exchange-K scheme: `"semitrans"` occupied-side semi-transformation contraction (exact, default); `"dm"` exact batched density-driven; `"lowrank"` per-vector SVD low-rank approximation (threshold `tddft_svd_tol`, generally not recommended) |
| `tddft_fxc_driver` | String | `"semitrans"` | XC-kernel scheme: `"semitrans"` folds $C_{vir}$ into the amplitudes and contracts the virtual side against the raw AO on grid (default); `"mo"` caches the occupied-side grid projections with a streamed virtual side (the MO-mode fxc algorithm); `"dm"` the assembled-transition-density NIMatmul fallback |
| `tddft_svd_tol` | f64 | 1.0e-6 | Relative singular-value threshold for the `"lowrank"` K driver |

### SCF stability analysis

| Keyword | Type | Default | Description |
|---|---|---|---|
| `stability` | String | `"off"` | SCF wave-function stability analysis: `"internal"` (RHF/RKS singlet or UHF/UKS orbital Hessian), `"external"` (RHF→UHF triplet channel check, RHF/RKS only), `"full"` (internal + external), `"auto"` (recommended; every check applicable to the reference). Any value other than `"off"` makes this task **stability-only — no excited-state calculation is performed** |
| `stability_nroots` | usize | 3 | Number of lowest stability-Hessian eigenvalues to solve for |
| `stability_tol` | f64 | 1.0e-8 | Davidson convergence threshold for the stability Hessian |

NOTE: the `check_stab` keyword in the top-level `[ctrl]` block accepts the same values as `stability`. When both are present, `[tddft] stability` wins. A stability-only deck needs no other excited-state settings.

### Analytic excited-state gradient

| Keyword | Type | Default | Description |
|---|---|---|---|
| `tddft_grad_state` | usize | 0 | Compute the analytic nuclear gradient for the N-th excited state (1-based); 0 disables. Restricted references only; the TDDFT solve must run first in the same task (`tddft_grad_state` must not exceed `nroots`) |

With `tddft_grad_state` set, the total gradient of gradient tasks (`jobtype = force`) and geometry optimizations automatically includes the response contribution of the selected state. With `tddft_spin = "both"` the gradient uses the singlet channel.

### PySOC export

| Keyword | Type | Default | Description |
|---|---|---|---|
| `pysoc` | bool | false | Export the TDDFT results to `rest_pysoc_export.json` (PySOC format) for subsequent spin-orbit coupling calculations. Requires a restricted reference with `tddft_spin = "both"` (spin-orbit coupling needs both the singlet and triplet amplitudes) |

### Response TDDFT settings

| Keyword | Type | Default | Description |
|---|---|---|---|
| `response_tddft` | bool | false | Enable the response TDDFT calculation (MO implementation, restricted references) |
| `response_tddft_solver` | String | `"klopper"` | Response-equation solver: one of `"pople"`, `"gmres"`, `"klopper"`, `"dense"` |
| `response_tddft_tol` | f64 | 1.0e-6 | Iterative solver convergence threshold |
| `response_tddft_max_iter` | usize | 200 | Maximum iterations of the iterative solver |
| `external_field_freq` | f64 | 0.5 | External field frequency $\omega$ (Hartree) |
| `lifetime_gamma` | f64 | 0.001 | Lifetime broadening $\gamma$ (Hartree) |

### Response TDDFT spatial grid sampling

When `response_tddft = true`, the following keywords generate the spatial sampling grid for the polarization-density export:

| Keyword | Type | Default | Description |
|---|---|---|---|
| `response_tddft_x_start` | f64 | 0.0 | X start coordinate (Bohr) |
| `response_tddft_x_end` | f64 | 1.0 | X end coordinate (Bohr) |
| `response_tddft_x_points` | usize | 2 | Number of X sampling points |
| `response_tddft_y_start` | f64 | 0.0 | Y start coordinate (Bohr) |
| `response_tddft_y_end` | f64 | 1.0 | Y end coordinate (Bohr) |
| `response_tddft_y_points` | usize | 2 | Number of Y sampling points |
| `response_tddft_z_start` | f64 | 0.0 | Z start coordinate (Bohr) |
| `response_tddft_z_end` | f64 | 1.0 | Z end coordinate (Bohr) |
| `response_tddft_z_points` | usize | 2 | Number of Z sampling points |

NOTE: the sampling points are generated with X as the outer, Y the middle, and Z the inner loop. Directions with ≤ 1 points are not sampled.

### FEAST solver settings

The following keywords take effect when `tddft_feast_solver = true` (the contour-integration solver is supported for MO mode and restricted references only):

| Keyword | Type | Default | Description |
|---|---|---|---|
| `tddft_feast_solver` | bool | false | Use the FEAST contour-integration solver |
| `tddft_feast_eigenrange_min` | f64 | 0.0 | Lower bound of the energy window (Hartree) |
| `tddft_feast_eigenrange_max` | f64 | 0.5 | Upper bound of the energy window (Hartree) |
| `tddft_feast_m_expected` | usize | 20 | Expected number of eigenvalues inside the window |
| `tddft_feast_max_iter` | usize | 30 | Maximum FEAST outer iterations |
| `tddft_feast_tol` | f64 | 1.0e-8 | FEAST convergence threshold |
| `tddft_feast_gmres_restart` | usize | 200 | GMRES restart dimension |
| `tddft_feast_gmres_max_iter` | usize | 500 | Maximum GMRES inner iterations |
| `tddft_feast_cg_max_iter` | usize | 100 | Maximum CG inner iterations |
| `tddft_feast_cg_tol` | f64 | 1.0e-8 | CG convergence threshold |
| `tddft_feast_init_guess_type` | String | `"random"` | Initial guess type |
| `tddft_feast_gaussian_width_factor` | f64 | 0.5 | Gaussian width factor |

## Examples

### Example 1: TDA singlet excitations of H2

```toml
[ctrl]
     print_level =               2
     num_threads =               4
     xc =                        "pbe0"
     basis_path =                "def2-tzvp"
     auxbas_path =               "def2-tzvp-rifit"
     eri_type =                  "ri-v"
     use_ri_symm =               true
     charge =                    0.0
     spin =                      1.0
     spin_polarization =         false
     mixer =                     "diis"
     max_scf_cycle =             100
     scf_acc_rho =               1.0e-8
     initial_guess =             "sad"

[tddft]
     tddft_method =              "tda"
     nroots =                    3

[geom]
     name = "H2"
     unit = "angstrom"
     position = '''
        H  0.0  0.0  0.0
        H  0.0  0.0  0.74
     '''
```

The program first converges the SCF, then calls the TDDFT module to compute the lowest 3 excited-state energies and oscillator strengths. The output prints, per state, the excitation energy (eV), oscillator strength, transition-dipole square, etc., and the energies/oscillator strengths are summarized into the `"tddft"` field of `rest_results.json`.

### Example 2: triplet excitations (AO mode)

```toml
[tddft]
     tddft_method =              "tda"
     tddft_mode =                "ao"
     tddft_spin =                "triplet"
     nroots =                    6
```

### Example 3: SCF stability analysis

```toml
[tddft]
     stability =                 "auto"
     stability_nroots =          3

[ctrl]
     check_stab =                "auto"
```

NOTE: either `stability` or `check_stab` being other than `"off"` triggers the stability analysis (the former wins); the task only performs the stability check and does not compute excited states. The lowest eigenvalue $\lambda_{\min} < -10^{-5}$ reports an instability; the results are both printed and written to the `"stability"` field of `rest_results.json`.

### Example 4: analytic excited-state gradient

```toml
[tddft]
     tddft_method =              "tda"
     nroots =                    3
     tddft_grad_state =          1
```

In a gradient task (`jobtype = force`), the program solves the TDDFT first and then adds the response contribution of the first excited state to the ground-state gradient, printing the total gradient. The same gradient feeds geometry optimizations.

### Example 5: unrestricted TDDFT

```toml
[ctrl]
     spin_polarization =         true
     spin =                      1.0

[tddft]
     tddft_method =              "tda"
     nroots =                    6
```

The unrestricted output covers both spin channels (e.g. `#3a->#5b` means an excitation from α orbital 3 to β orbital 5), with post-processing conventions matching PySCF's UHF-TDDFT.

## Notes

- **RI acceleration prerequisite**: the TDDFT module relies on RI acceleration. The deck must provide `auxbas_path` and set `eri_type = "ri-v"`. The short-range exchange integrals of RSH functionals are built automatically during the SCF; no extra settings are needed.
- **Solver selection**: for small systems (restricted excitation space dim ≤ 15; unrestricted MO mode TDA ≤ 15, Full LR ≤ 80) the program automatically uses dense diagonalization (full LR diagonalizes the non-Hermitian $[\mathbf{A}\ \mathbf{B}; -\mathbf{B}\ -\mathbf{A}]$ directly); for medium and larger systems the Davidson iterative solver is the default (batched interface in AO mode); the FEAST solver is available for specific energy windows (restricted references + MO mode only).
- **Singlet/triplet excitations**: controlled by `tddft_spin` for restricted references; triplets and `"both"` require `tddft_mode = "ao"`. Unrestricted references have a single spin-coupled channel and do not accept `tddft_spin`.
- **Stability analysis**: mutually exclusive with the excited-state run; a DFT reference requires the numerical grids (HF references work without grids and automatically evaluate the RI J/K parts only); ROHF references are not supported; the real→complex and UHF→GHF external checks are not yet implemented.
- **Analytic gradient**: restricted references only; the TDDFT solve must run first in the same task; the gradient is analytic within the RI approximation.
- **Response TDDFT**: the response-mode linear system is a 4-component non-Hermitian system with 4 times the dimension of the eigenvalue mode. The Klopper subspace solver is the recommended choice.
