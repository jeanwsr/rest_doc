# DFA Definition API (DFAdef / DFAComponent)

The DFA-definition API centers on `DFAdef` / `DFAComponent` (design principles in [ADR-0003](../adr/adr-0003-dfa-definitions.md)): everything is a component, and the definition is decoupled from the computation. The parse syntax is documented in [parse_xc](parse_xc.md); the libxc interface in [the libxc interface](libxc.md).

## DFAdef

```rust
pub struct DFAdef {
    pub xc_scf: Option<Vec<DFAComponent>>,        // component list of the SCF part
    pub xc_nscf: Option<Vec<DFAComponent>>,       // component list of the final-energy (post-SCF) part
    pub dfa_hybrid_scf: f64,                      // total HFX coefficient of the SCF functional
    pub dfa_rsh_scf: (Option<f64>, f64, f64),     // SCF RSH raw parameters (omega, c_LR, c_SR - c_LR)
    pub dfa_rsh_scf_lrsr: Option<(f64, f64, f64)>,// SCF RSH parameters (omega, c_LR, c_SR), computed at parse time
    pub dfa_hybrid_nscf: Option<f64>,             // final-energy HFX coefficient (two-step functionals only)
    pub spin_channel: usize,
    pub reference: Vec<String>,
}
```

When both `xc_scf` and `xc_nscf` are present it is a two-step (double-hybrid etc.) functional; single-step functionals have only the former as `Some`.

### Constructors (parsing / definition-construction phase)

| Method | Description |
|---|---|
| `DFAdef::new()` | an empty definition |
| `parse_and_derive(xc, spin, print)` | the parse entry: xc string -> `DFAdef` |
| `from_dfa4rest(&DFA4REST)` | the legacy/DeepLearning adapter (code arrays -> components; the only conversion direction on the production path) |
| `to_dfa4rest()` | convert into the old DFA4REST representation (full): used solely by the round-trip parity tests; no production consumer |
| `to_xc_data()` | SCF-side-only DFA4REST: test-only, planned for removal |
| `update_pt2_params(os, ss)` | override the PT2 `(c_os, c_ss)` factors (ctrl keywords; also called during the DL-DFA iteration) |
| `get_hybrid_scf(spin)` / `get_hybrid_nscf(spin)` | derive the total HFX coefficient of the SCF / final-energy functional from the components (written into the stored fields at parse time) |
| `get_rsh_scf(spin)` / `get_rsh_scf_lrsr(spin)` | derive the RSH parameters from the components (raw `(omega, c_LR, c_SR - c_LR)` / consumer convention `(omega, c_LR, c_SR)`; written into the stored fields at parse time) |
| `check_sanity()` | definition self-check |

### Accessors (consumer reads)

SCF/grad/hessian/analdrv consumers read the definition through the following methods, without modifying it.

| Method | Returns | Description |
|---|---|---|
| `is_dfa_scf()` | `bool` | whether the SCF functional evaluates XC on a grid (has libxc components) |
| `is_hybrid()` / `is_rsh()` | `bool` | hybrid / range-separated predicates |
| `rsh_params_scf()` | `Option<(f64, f64, f64)>` | SCF RSH parameters `(omega, c_LR, c_SR)` (stored at parse time; `omega()` reads the single parameter) |
| `family_pos()` | `Option<DFAFamily>` | post-SCF family (PT2/RPA/SCSRPA/SBGE2), backed by `FIFTH_RUNG_COMPONENT_TYPES` |
| `fifth_rung_component()` / `is_fifth_dfa()` / `is_rpa()` | - | the fifth-rung component and its predicates |
| `is_nlc()` / `vv10_params()` | - | non-local correlation (VV10) predicate and `(b, C)` parameters |
| `x_part(spin)` / `c_part(spin)` | `Option<DFAdef>` | pure exchange / correlation part (decomposed by the libxc functional kind; combined functionals yield `None`) |
| `has_nscf()` / `use_density_gradient()` / `use_kinetic_density()` / `use_laplacian()` | `bool` | definition-shape predicates |
| `summary()` / `formatted_output()` | - | printout |

### Definition -> evaluation list (boundary construction)

The boundary between definition and computation: the `(factor, LibXCFunctional)` list is built here (libxc functionals initialized once, with `_omega` set) and then consumed by the exc/vxc free functions.

| Method | Description |
|---|---|
| `scf_func_list(spin)` | the SCF functional's evaluation list |
| `nscf_func_list(spin)` | the DFT part of the final-energy functional (`None` for pure-MP2 families) |
| `func_list_from_components(comps, spin)` | the door for an arbitrary component list (used by the response/grad interfaces) |

## DFAComponent

```rust
pub struct DFAComponent {
    pub factor: f64,
    pub func: String,
    pub func_full_name: String,
    pub id: usize,
    pub param_positional: Vec<f64>,
    pub param_keyword: HashMap<String, Value>,
    pub component_type: ComponentType,
    pub reference: Vec<String>,
}

pub enum ComponentType { HF, RSHF, PT2, RPA, SCSRPA, SBGE2, Disp, Libxc, Unknown }
```

### Constructors / normalization (parse phase)

| Method | Description |
|---|---|
| `DFAComponent::new(factor, func)` / `from_libxc(id, factor)` | construct a component |
| `canonicalize()` / `normalize_rsh()` / `normalize_pt2_param()` | parse-time normalization (RSH convention, PT2 parameter shape) |

### Accessors (consumer reads)

| Method | Description |
|---|---|
| `to_libxc_functional(spin)` | initialize this component's `LibXCFunctional` (`None` for non-Libxc components) |
| `get_hybrid(spin)` / `get_rsh(spin)` | per-component HFX / RSH parameters |
| `is_rsh()` / `is_nlc()` | component predicates |
| `get_dftd_params()` | Disp component parameters (dftd3-rs / dftd4-rs toml) |
| `formatted_output()` | component printout |

The non-libxc component types, their input names, aliases and parameters are listed in [parse_xc](parse_xc.md).

## Computation API (free functions)

Definition and computation are decoupled: the computation functions live under `dft` (`mod.rs` etc.) and receive `(factor, LibXCFunctional)` lists - never a definition object or `(factor, id)`.

| Function | Description |
|---|---|
| `xc_exc_vxc(&funcs, grids, spin, dm, mo, occ, print_level)` | XC energy + potential |
| `xc_exc_vxc_slots(&funcs, range_grids, ...)` / `_dm_only` | blocked (screened) variants |
| `xc_exc(&funcs, grids, ...)` | XC energy only |
| `post_xc_exc(&names, ...)` / `post_tabulated_exc(&funcs, ...)` | post-SCF functional energies |
| `DFAdef::scf_func_list(spin)` / `nscf_func_list(spin)` / `func_list_from_components` | definition -> evaluation-list builders |
| `xc_part_of_name(name, spin, kind)` / `DFAdef::x_part` / `c_part` | name/definition -> pure X or pure C part |

### XCData (the derivative-time bundle)

`XCData` (`src/dft/num_int.rs`) is the gradient/Hessian-side evaluation bundle: it holds the func list + `xc_type` + the density sources (`dm` or `mo_coeffs`/`occ`). Production chain: `mol.dfadef` -> `scf_func_list(spin)` / `func_list_from_components` -> `XCData` built via `NumInt::gen_xc_data` (`grad/rks.rs`, `grad/uks.rs`) or `xc_skeleton_dao` (`grad/rdh.rs`, `grad/udh.rs`) -> `get_vxc_rayon_new` / `eval_xc_eff`.

Note: `DFAdef::to_xc_data()` (the old DFA4REST conversion, planned for removal) is unrelated to `XCData` - same name, different things.

## Correspondence with the old DFA4REST API

`DFA4REST` (now at `src/dft/legacy.rs`) was the old code-array representation; the legacy/DeepLearning routes convert their parse result into a DFAdef via `DFAdef::from_dfa4rest`. Status markers: → equivalent replacement; ≈ replaced but not identical; ✕ removed; L legacy-only; N new.

| Old (DFA4REST) | New (DFAdef / free fn) | Status |
| --- | --- | --- |
| `dfa_compnt_scf: Vec<usize>` + `dfa_paramr_scf: Vec<f64>` | `xc_scf: Option<Vec<DFAComponent>>` | ≈ |
| `dfa_compnt_pos` + `dfa_paramr_pos` | `xc_nscf: Option<Vec<DFAComponent>>` | ≈ (pos arrays held only Libxc codes; xc_nscf holds all non-SCF components) |
| `dfa_paramr_adv: Option<Vec<f64>>` (os/ss) | `fifth_rung_component().param_positional` | ≈ |
| `dfa_hybrid_scf` / `dfa_hybrid_nscf` (old `dfa_hybrid_pos`) | same-name fields + `get_hybrid_scf/get_hybrid_nscf` | → |
| `dfa_rsh_scf` (consumer triple) | `dfa_rsh_scf` (raw) + `dfa_rsh_scf_lrsr` (consumer) + `rsh_params_scf()` | ≈ |
| `dfa_family_pos` (field) | `family_pos()` (accessor, `FIFTH_RUNG_COMPONENT_TYPES`) | ≈ |
| `nlc_vv10: Option<(b, c)>` | `is_nlc()` / `vv10_params()` | → |
| `is_dfa_scf()/is_hybrid()/is_rsh()/is_fifth_dfa()/is_rpa()` | same names | → |
| `use_density_gradient()/use_kinetic_density()` | same names | → |
| `omega()` | `omega()` (field read) | → |
| `rsh_params()` (re-derived) | `rsh_params_scf()` (parse-time field) | → |
| `rsh_alpha()` / `use_eri()` | - | ✕ |
| `xc_exc_vxc(&self, ...)` and other exc/vxc methods | `dft::xc_exc_vxc(&funcs, ...)` free fns | ≈ (funcs currency) |
| `libxc_func_list()` | `DFAdef::scf_func_list(spin)` | → |
| `get_hybrid_libxc` / `get_rsh_libxc` | `DFAComponent::get_hybrid` / `get_rsh` | ≈ (parity-tested) |
| `libxc_code_fdqc` / `xc_func_init_fdqc` | `xc_part_of_name` / `x_part` / `c_part` | ≈ |
| `summary(print_level)` / `describe()` / `formatted_output()` | `summary()` / `formatted_output()` | ≈ (component view) |
| `new(xc, spin, print)` (legacy parser) | `parse_and_derive(xc, spin, print)` | ≈ |
| `new_nonstandard(...)` | - | ✕ (NonStandard input removed) |
| - | `x_part` / `c_part` / `xc_part_of_name` / `from_libxc` | N |

## See also

- [ADR-0003: unify DFA definitions into DFAdef](../adr/adr-0003-dfa-definitions.md)
- [parse_xc: parse syntax](parse_xc.md)
- [the libxc interface](libxc.md)
