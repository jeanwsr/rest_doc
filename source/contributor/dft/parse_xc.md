# parse_xc

To handle arbitrary functional inputs, a new xc_parser was written following PySCF's `parse_xc` and the extensions from ajz34/dh. In REST it is used as follows:

```
xc = "B3LYP"
xc_parser = "parse_xc"
```

The entire input string is **case-insensitive**.

## Level 1 Structure

The xc string can be split by `,` into two regions, or kept as a single region. More than two regions are not allowed. (Commas inside `()` are ignored.)

Examples: `0.5*PBE + 0.5*B88, PBE` or `B3LYP`.

With two regions, the left is the X (exchange) region and the right is the C (correlation) region. With one region, it is called the XC region.

## Level 2 Structure

Each region can be written as a linear combination of functional components. Each term has an optional coefficient before the component name, connected by `*`. If no coefficient is given, it defaults to `1.0`.

Example: `PBE - 0.1*B88 + 23.45*PBEsol`.

## Level 3 Structure

Rules for each functional component.

### Name

Resolved by the following priority:

1. **WHITELIST_NONLIBXC**: e.g. HF, MP2, SCSRPA, RPA, SBGE2, DFTD3, DFTD4, etc.
2. **Predefined CODES**: e.g. `LDA -> 1`, `B3LYPG -> 402`. These directly give the libxc name and id.
3. **ALIAS**: recursively expanded. e.g. `PBE -> PBE,PBE`. When there are two regions, X/C region aliases cannot expand into functionals of a different type. For example `"BLYP,"` is invalid but `BLYP` (without comma) is valid:

```
Parsing xc: blyp
SCF components:
  type: Libxc, factor: 1, func: B88, id: 106, func_full_name: GGA_X_B88,
    reference: ["A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)"]
  type: Libxc, factor: 1, func: LYP, id: 131, func_full_name: GGA_C_LYP,
    reference: ["C. Lee, W. Yang, and R. G. Parr.,  Phys. Rev. B 37, 785 (1988)", "B. Miehlich, A. Savin, H. Stoll, and H. Preuss.,  Chem. Phys. Lett. 157, 200 (1989)"]
  Total hybrid: 0
```

4. **Other regular libxc functionals**:
   - Can be a libxc full name, e.g. `GGA_X_PBE`. In the X region, C-type names are not allowed, and vice versa. The XC region accepts both X and C types.
   - Can have a prefix `X_`, `C_`, or `XC_`. Again, prefix must match the region.
   - Can have no prefix, e.g. `PBE`. In this case, the appropriate prefix (`X_`, `C_`, or `XC_`) is inferred from the region.

For cases 3 and 4, all possible full names are searched. If multiple matches are found, an error is raised.

Examples:

```
Parsing xc: .2*hf + 0.08*lda + 0.72*b88, 0.81*lyp + 0.19*vwn3
Detected REST_DATA_DIR environment variable: /home/wsr/rest_test/xc
Loaded JSON file: "/home/wsr/rest_test/xc/xc_user.json"
SCF components:
  type: HF, factor: 0.2, func: HF,
    reference: Not available yet
  type: Libxc, factor: 0.08, func: LDA, id: 1, func_full_name: LDA_X,
    reference: ["P. A. M. Dirac.,  Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)", "F. Bloch.,  Z. Phys. 57, 545 (1929)"]
  type: Libxc, factor: 0.72, func: B88, id: 106, func_full_name: GGA_X_B88,
    reference: ["A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)"]
  type: Libxc, factor: 0.81, func: LYP, id: 131, func_full_name: GGA_C_LYP,
    reference: ["C. Lee, W. Yang, and R. G. Parr.,  Phys. Rev. B 37, 785 (1988)", "B. Miehlich, A. Savin, H. Stoll, and H. Preuss.,  Chem. Phys. Lett. 157, 200 (1989)"]
  type: Libxc, factor: 0.19, func: VWN3, id: 8, func_full_name: LDA_C_VWN_RPA,
    reference: ["S. H. Vosko, L. Wilk, and M. Nusair.,  Can. J. Phys. 58, 1200 (1980)"]
  Total hybrid: 0.2
```

Invalid case:

```
Parsing xc: 0.5*pbe + 0.5*b88, x_pbe

thread 'main' (1822079) panicked at rest/src/dft/parse_xc/parse.rs:130:17:
Error: functional X_PBE has illegal prefix for type C
```

### NONLIBXC Components

| Name | ComponentType | Alias | Parameters |
| --- | --- | --- | --- |
| HF | HF | | |
| SR_HF | RSHF | | 1 positional (`omega`) |
| LR_HF | RSHF | | 1 positional (`omega`) |
| RSH | RSHF | | 3 positional (`omega`, `alpha`, `beta`) |
| MP2 | PT2 | | `os=, ss=` or 2 positional |
| MP2_OS | PT2 | | 1 positional |
| MP2_SS | PT2 | | 1 positional |
| RPA | RPA | | |
| SCSRPA | SCSRPA | | like MP2 |
| SBGE2 | SBGE2 | | like MP2 |
| DFTD3 | Disp | D3 | toml-param of dftd3-rs |
| DFTD4 | Disp | D4 | toml-param of dftd4-rs |
| VV10 | Disp | | |

Note: VV10 is currently only parsed; its non-local correlation calculation is not yet implemented.

### Parameters

Each component can have parameters inside `()`. Parameters must be **all positional** or **all keyword** (mixed usage is not allowed). Positional: `(0.3, 0.4)`. Keyword: `(a=0.3, b=0.4)`.

Currently, parameters for libxc components are parsed but have no functional effect.

MP2, SCSRPA, etc. can specify parameters as described in the table above.

D3/D4 also support custom parameters; see [dftd3-rs](https://github.com/RESTGroup/dftd3-rs#example-custom-parameters-by-toml).

### Dispersion

Dispersion suffix shorthand is supported. A suffix appended to a functional name is automatically expanded into a Disp component:

| Suffix | Expansion |
| --- | --- |
| `-D3ZERO` | `D3(version=zero, xc=...)` |
| `-D3BJ` | `D3(version=bj, xc=...)` |
| `-D4` | `D4(version=bj, xc=...)` |
| `-V` / `-VV10` | `VV10(xc=...)` |
| `-RVV10` | `rVV10(xc=...)` |

For example, `B3LYP-D3BJ` is equivalent to `B3LYP + D3(version=bj, xc=b3lyp)` (though this equivalence does not always hold; see special cases below). When `xc` contains the above suffixes, it cannot coexist with the `empirical_dispersion` setting.

Certain special dispersion functionals also support direct alias expansion:

| Name | Expansion |
| --- | --- |
| `CF22D` | `CF22D + D3(version=zero, xc=cf22d)` |
| `WB97X-D3BJ` | `wb97x_v + D3(version=bj, xc=wb97x)` |
| `WB97X-D3` | `wb97x_d3 + D3(version=zero, xc=wb97x)` |
| `WB97X-V` | `wb97x_v + VV10(xc=wb97x_v)` |
| `WB97M-D3BJ` | `wb97m_v + D3(version=bj, xc=wb97m)` |
| `WB97M-V` | `wb97m_v + VV10(xc=wb97m_v)` |
| `SCAN-VV10` | `SCAN,SCAN_VV10 + VV10(xc=scan_vv10)` |
| `SCAN-RVV10` | `SCAN,SCAN_RVV10 + VV10(xc=scan_rvv10)` |
| `REVSCAN-VV10` | `REVSCAN,REVSCAN_VV10 + VV10(xc=revscan_vv10)` |

## Multi-Step Functionals

Built-in JSON files (in a format similar to ajz34/dh) support several multi-step functionals. Multi-step parsing has the highest priority: the parser first checks if the input is a multi-step functional, and then proceeds with either multi-step or single-step parsing accordingly.

Example:

```
Parsing xc: xyg3
Detected REST_DATA_DIR environment variable: /home/wsr/rest_test/xc
Loaded JSON file: "/home/wsr/rest_test/xc/xc_user.json"
Detected multi-step functional xyg3, which is parsed to:
Step for SCF         : B3LYPg
Step for final energy: 0.8033*HF - 0.0140*LDA + 0.2107*B88, 0.6789*LYP + 0.3211*MP2
SCF components:
  type: Libxc, factor: 1, func: B3LYPG, id: 402, func_full_name: HYB_GGA_XC_B3LYP,
    reference: ["P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and M. J. Frisch.,  J. Phys. Chem. 98, 11623 (1994)"]
  Total hybrid: 0.2
Final energy components:
  type: HF, factor: 0.8033, func: HF,
    reference: Not available yet
  type: Libxc, factor: -0.014, func: LDA, id: 1, func_full_name: LDA_X,
    reference: ["P. A. M. Dirac.,  Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)", "F. Bloch.,  Z. Phys. 57, 545 (1929)"]
  type: Libxc, factor: 0.2107, func: B88, id: 106, func_full_name: GGA_X_B88,
    reference: ["A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)"]
  type: Libxc, factor: 0.6789, func: LYP, id: 131, func_full_name: GGA_C_LYP,
    reference: ["C. Lee, W. Yang, and R. G. Parr.,  Phys. Rev. B 37, 785 (1988)", "B. Miehlich, A. Savin, H. Stoll, and H. Preuss.,  Chem. Phys. Lett. 157, 200 (1989)"]
  type: PT2, factor: 1, func: MP2,
    param_positional: [0.3211, 0.3211], reference: Not available yet
  Total hybrid: 0.8033
References:
  10.1073/pnas.0901093106
```

The list of supported double-hybrid functionals can be found in `src/dft/parse_xc/*json`.

Users can define custom multi-step functionals by setting the `REST_DATA_DIR` environment variable and creating `xc_*.json` files in that directory:

```json
{
    "B2GP-PLYP-D4": {
        "code": "0.65*HF + 0.35*B88, 0.64*LYP + 0.36*MP2 + D4(xc=b2gpplyp)"
    },
    "R-XYG3": {
        "code_scf": "B3LYPg",
        "code": "0.8033*HF - 0.0140*LDA + 0.2107*B88, 0.6789*LYP + 0.3211*SBGE2"
    }
}
```

Support for other multi-step functional input formats (e.g. `String`, `Vec<String>`) is under consideration.

## Result

The parsed result is a `DFAdef` struct.

```
pub struct DFAdef {
    pub xc_scf: Option<Vec<DFAComponent>>,
    pub xc_nscf: Option<Vec<DFAComponent>>,
    pub reference: Vec<String>,
    pub spin_channel: usize,
    pub dfa_hybrid_scf: f64,                    // SCF HFX coefficient
    pub dfa_rsh_scf: (Option<f64>, f64, f64),   // SCF RSH parameters (omega, alpha, beta)
    pub dfa_hybrid_nscf: Option<f64>,           // final-energy HFX coefficient (two-step only)
}
```

When both `xc_scf` and `xc_nscf` are present, it is a two-step functional; single-step functionals have only the former as `Some`.

The `DFAComponent` struct is defined as:

```rust
pub struct DFAComponent {
    pub factor: f64,
    pub func: String,
    pub func_full_name: String,
    pub id: usize,
    pub param_positional: Vec<f64>,
    pub param_keyword: HashMap<String, Value>,   // Value can be f64 or String
    pub component_type: ComponentType,
    pub reference: Vec<String>,
}

pub enum ComponentType {
    HF,
    RSHF,      // SR_HF, LR_HF, RSH
    PT2,
    RPA,
    SCSRPA,
    SBGE2,
    Disp,
    Libxc,
    Unknown,
}
```
