# ADR-0003: DFA Definitions

- **Status**: In use
- **Date**: 2026-11-03
- **Scope**: The DFA definition (`src/dft/parse_xc/`) of crate `rest` and all DFA consumers (`scf_io`, `grad`, `hessian`, `analdrv`, `ri_tddft`, `ri_gw`, `ri_pt2`, etc.). Legacy code still contains `DFA4REST` (moved to `src/dft/legacy.rs`) and `Molecule.xc_data` (now an `Option`), kept for debugging; both will be retired together with the legacy parser.
- **Author**: Shirong Wang

## Decision

- **Principle 1: the definition and the computation are decoupled.** All definition work is completed before entering the computation. Computation functions do not do definition-related work; for example, they receive `(factor, LibXCFunctional)` rather than `(factor, id)`.
- **Principle 2: everything is a component.** Every part of a DFA is defined as a component, including D3/NLC/HF.
  - Component parameters are not added as fields of the whole; they are obtained from the component's (`DFAComponent`) API/fields. (That said, for compatibility with the old API and for implementation convenience, `DFAdef` does keep a few whole-functional fields, such as the overall hybrid/RSH parameters.)
  - Whole-functional parameters are recommended to be read through accessor functions, e.g. `rsh_params_scf()`.
- `DFAdef` (`src/dft/parse_xc/`) is the primary definition of a DFA.
  - `DFA4REST` moves to `src/dft/legacy.rs`; the legacy and DeepLearning routes convert their parse result into a DFAdef via `DFAdef::from_dfa4rest`;
  - computation functions such as `xc_exc_vxc` become free functions.

The core pattern of this decision (the full field/method tables are on the [DFA Definition API](../dft/dfa-definitions.md) page):

```rust
let dfa = parse_and_derive("xyg3", 1, 0);

// component access: iterate the SCF / final-energy component lists
for comp in dfa.xc_scf.as_ref().unwrap() {
    // component fields: factor, func, id, param_positional, param_keyword, component_type, reference
    println!("{} {} {}", comp.factor, comp.func, comp.id);
}

// the fifth-rung (post-SCF) component
let pt2 = dfa.fifth_rung_component().unwrap();
assert_eq!(pt2.component_type, ComponentType::PT2);
assert_eq!(pt2.param_positional, [0.3211, 0.3211]); // (c_os, c_ss)

// functional-level parameters: read through accessors, not re-derived from components
let (omega, c_lr, c_sr) = dfa.rsh_params_scf().unwrap();
let hyb = dfa.dfa_hybrid_scf;
```

## 1. DFAComponent and ComponentType

Every part of a DFA is a component, `DFAComponent`:

| field | description |
|---|---|
| `factor` | coefficient |
| `func` / `func_full_name` | component name / full libxc name |
| `id` | libxc number (Libxc components) |
| `param_positional` / `param_keyword` | parameters (e.g. PT2 (c_os, c_ss), RSHF omega) |
| `component_type` | component type |
| `reference` | references |

ComponentType values: Libxc (a libxc functional), HF (exact exchange), RSHF (range-separated HF: SR_HF/LR_HF/RSH), PT2 (MP2-type correlation), RPA / SCSRPA / SBGE2 (fifth-rung RPA families), Disp (empirical dispersion D3/D4/VV10), Unknown.

The full field/method tables are on the [DFA Definition API](../dft/dfa-definitions.md) page.

## 2. Background

The original purpose was to make `parse_xc` the default parser, avoiding the trouble of manually selecting a parser. Since `DFA4REST` and `DFAdef` are both structs describing a DFA definition (the former is not only a definition - it also carried computation functions - but its own fields are definition only), the DFA-definition-related APIs are unified into `DFAdef`. `Molecule` now holds both `dfadef` (always present) and `xc_data: Option<DFA4REST>` (`Some` only on the legacy/DeepLearning routes); all DFA consumers read the definition through the accessors of `mol.dfadef` (`is_dfa_scf`, `is_hybrid`, `is_rsh`, `family_pos`, `fifth_rung_component`, `rsh_params_scf`, `is_nlc`, `vv10_params`, ...).

## 3. Alternatives

The old API has not been deleted: all `DFA4REST` functions and the `Molecule.xc_data` field still exist, although no consumer remains on the normal call path. Temporarily restoring the old API for debugging is therefore feasible; after such debugging, new implementations are expected to be written within the new framework.

## 4. Notes

Specific API details may still change; this ADR only constrains the design principles above. The old/new API correspondence is documented in [the libxc interface](../dft/libxc.md) and [parse_xc](../dft/parse_xc.md).
