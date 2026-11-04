# libxc Interface

REST uses the external crate [`libxc-rs`](https://github.com/RESTGroup/libxc-rs) ([docs.rs](https://docs.rs/libxc/latest/libxc/)) as a safe Rust wrapper around the C libxc library. It is dynamically loaded at runtime via the Cargo feature `dynamic_loading`; no static linking is needed.

## Architecture

The call hierarchy, top-down:

1. `DFAdef` (`parse_xc`, the primary DFA definition, see [ADR-0003](../adr/adr-0003-dfa-definitions.md); the old `DFA4REST` lives in `legacy.rs`, legacy route only) and `libxc_itrf.rs` (numerical integration interface providing high-order derivatives for num_int/response), calling via `xc_func_init` / `lda/gga/mgga_exc_vxc` / `eval_libxc_func_new`
2. `libxc_helper.rs`: a thin wrapper hiding the libxc-rs details
3. the external `libxc` crate: safe Rust bindings FFI-calling the C libxc

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

### `DFAdef` (`src/dft/parse_xc/parse.rs`, primary definition)

The primary definition of a DFA (principles in [ADR-0003](../adr/adr-0003-dfa-definitions.md); the full API table is on [DFA Definition API](dfa-definitions.md)): the SCF part and the non-SCF (post-SCF) part are both lists of `DFAComponent`s; whole-functional parameters are recommended to be read through accessors.

| Method | Description |
|---|---|
| `scf_func_list(spin)` / `nscf_func_list(spin)` | Build the `(factor, LibXCFunctional)` evaluation list (with `_omega` set), consumed by the exc/vxc free functions |
| `is_dfa_scf()` / `is_hybrid()` / `is_rsh()` | SCF functional predicates |
| `family_pos()` / `fifth_rung_component()` / `is_rpa()` | Post-SCF (fifth-rung) family and component |
| `rsh_params_scf()` | RSH parameters `(omega, c_LR, c_SR)` of the SCF functional (computed at parse time, stored in `dfa_rsh_scf_lrsr`) |
| `omega()` / `update_pt2_params(os, ss)` | RSH omega; PT2 `(c_os, c_ss)` factor override |
| `is_nlc()` / `vv10_params()` | Non-local correlation (VV10) predicate and `(b, C)` parameters |
| `x_part(spin)` / `c_part(spin)` | Pure exchange / correlation parts (decomposed by the libxc functional kind) |

### `DFA4REST` (`src/dft/legacy.rs`, legacy)

The old code-array definition, kept only on the legacy/DeepLearning routes (its parse result is converted into a DFAdef via `DFAdef::from_dfa4rest`).

| Method | Description |
|---|---|
| `init_libxc(code)` | Initialize `LibXCFunctional` from func ID |
| `init_libxc_and_set_param(code)` | Initialize with parameters (e.g. omega) |
| `get_hybrid_libxc(...)` | Compute total HFX coefficient from SCF functionals |
| `get_rsh_libxc(...)` | Get RSH parameters from SCF functionals |

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

`eval_xc_eff(funcs, xc_type, spin, rho_array, np, deriv)` is the main entry point for `num_int.rs` and `response.rs`. It accepts a pre-built list of `(factor, LibXCFunctional)` pairs (libxc functionals are initialized once at the definition boundary rather than re-initialized per grid block), computes XC energy, first-order potential, and high-order derivatives, then converts from libxc's raw output format to REST's internal representation.

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
