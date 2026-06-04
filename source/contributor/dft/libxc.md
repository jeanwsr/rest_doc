# libxc Interface

REST uses the external crate [`libxc-rs`](https://github.com/RESTGroup/libxc-rs) ([docs.rs](https://docs.rs/libxc/latest/libxc/)) as a safe Rust wrapper around the C libxc library. It is dynamically loaded at runtime via the Cargo feature `dynamic_loading`; no static linking is needed.

## Architecture

```
┌─────────────────────┐
│  DFA4REST (mod.rs)   │  Core DFT data structure, calls libxc via DFA
│  libxc_itrf.rs       │  Numerical integration interface (high-order derivatives)
└────────┬────────────┘
         │ xc_func_init / lda/gga/mgga_exc_vxc / eval_libxc_func_new
┌────────▼────────────┐
│  libxc_helper.rs     │  Thin wrapper over libxc-rs
└────────┬────────────┘
         │ LibXCFunctional / compute_lda/gga/mgga
┌────────▼────────────┐
│  libxc (external)    │  Safe Rust FFI bindings to C libxc
└─────────────────────┘
```

## Key Types

### `LibXCFunctional` (external crate `libxc`)

Constructed via `LibXCFunctional::from_number(id, spin)` or REST's `xc_func_init(code, spin)`.

Implements `Drop` — resources are automatically freed when the object goes out of scope. No manual `xc_func_end()` needed.

| Method | Returns | Description |
|---|---|---|
| `family()` | `LibXCFamily` | Functional family (`LDA`/`GGA`/`MGGA`/`HybGGA`/`HybMGGA`/`HybLDA`) |
| `hyb_exx_coef()` | `Option<f64>` | Exact-exchange coefficient, `None` for pure functionals |
| `cam_coef()` | `Option<(f64,f64,f64)>` | CAM parameters `(omega, alpha, beta)` |
| `is_hyb_cam()` | `bool` | Whether it is a range-separated hybrid |
| `needs_tau()` | `bool` | Whether kinetic energy density is needed (MGGA) |
| `needs_laplacian()` | `bool` | Whether Laplacian is needed |
| `flags()` | `LibXCFlags` | Feature flags (e.g. `VV10`) |
| `references()` | `&[LibXCReference]` | Reference list |
| `describe()` | `String` | Human-readable description |

### `LibXCFamily` enum changes

| Old (hand-written) | New (external crate) |
|---|---|
| `HybridGGA` | `HybGGA` |
| `HybridMGGA` | `HybMGGA` |
| (none) | `HybLDA` |

### `DFA4REST` (`src/dft/mod.rs`)

| Method | Description |
|---|---|
| `init_libxc(code)` | Initialize `LibXCFunctional` from func ID |
| `init_libxc_and_set_param(code)` | Initialize with parameters (e.g. omega) |
| `get_hybrid_libxc(...)` | Compute total HFX coefficient from SCF functionals |
| `get_rsh_libxc(...)` | Get RSH parameters from SCF functionals |
| `xc_exc_vxc(...)` | Dispatch to `lda/gga/mgga_exc_vxc` by family |

## `libxc_helper.rs` Public API

| Function | Description |
|---|---|
| `xc_func_init(func_id, spin)` | Initialize functional from libxc id; spin=1 unpolarized / 2 polarized |
| `xc_code_fdqc(name)` | Map functional name (e.g. `"pbe"`) to `[corr, exch, mixed]` code array |
| `lda_exc_vxc(func, rho)` | LDA energy + potential |
| `gga_exc_vxc(func, rho, sigma)` | GGA energy + potential |
| `mgga_exc_vxc(func, rho, sigma, lapl, tau)` | MGGA energy + potential |
| `eval_libxc_func_new(func, spin, deriv, ...)` | General evaluator writing into pre-allocated buffer |

## `libxc_itrf.rs` Integration Interface

`eval_xc_eff(func_ids, factors, xc_type, spin, rho_array, np, deriv)` is the main entry point for `num_int.rs` and `response.rs`. It accepts multiple libxc functional IDs with coefficients, computes XC energy, first-order potential, and high-order derivatives, then converts from libxc's raw output format to REST's internal representation.

## API Migration Reference

| Old (`XcFuncType`) | New (`libxc-rs`) |
|---|---|
| `XcFuncType::xc_func_init(code, spin)` | `xc_func_init(code, spin)` (compatible) |
| `xc_func.xc_func_end()` | (removed, RAII) |
| `xc_func.xc_func_family` (field) | `xc_func.family()` (method) |
| `xc_func.xc_hyb_exx_coeff()` | `xc_func.hyb_exx_coef()` returns `Option` |
| `xc_func.xc_hyb_cam_coef()` | `xc_func.cam_coef()` returns `Option` |
| `xc_func.is_rsh()` | `xc_func.is_hyb_cam()` |
| `xc_func.use_laplacian()` | `xc_func.needs_laplacian()` |
| `xc_func.use_kinetic_density()` | `xc_func.needs_tau()` |
| `xc_func.get_libxc_family()` | `xc_func.family()` |
| `xc_func.get_libxc_references()` | `xc_func.references()` |
| `xc_func.xc_func_info_printout()` | `println!("{}", xc_func.describe())` |
| `unsafe { xc_version(...) }` | `libxc::util::libxc_version()` |
| `unsafe { xc_functional_get_name(id) }` | `libxc::util::libxc_functional_get_name(id)` |

## Dependency

```toml
libxc = { version = "0.1", features = ["api-v7_0", "dynamic_loading"] }
```
