# DFA 定义 API（DFAdef / DFAComponent）

DFA 定义的 API 以 `DFAdef` / `DFAComponent` 为中心（设计原则见 [ADR-0003](../adr/adr-0003-dfa-definitions.md)）：一切皆组分，定义与计算解耦。解析语法见 [parse_xc](parse_xc.md)；libxc 接口见 [libxc 接口](libxc.md)。

## DFAdef

```rust
pub struct DFAdef {
    pub xc_scf: Option<Vec<DFAComponent>>,        // SCF 部分的组分列表
    pub xc_nscf: Option<Vec<DFAComponent>>,       // 最终能量（post-SCF）部分的组分列表
    pub dfa_hybrid_scf: f64,                      // SCF 泛函的总 HFX 系数
    pub dfa_rsh_scf: (Option<f64>, f64, f64),     // SCF 泛函 RSH 原始参数 (omega, c_LR, c_SR - c_LR)
    pub dfa_rsh_scf_lrsr: Option<(f64, f64, f64)>,// SCF 泛函 RSH 参数 (omega, c_LR, c_SR)，parse 时计算
    pub dfa_hybrid_nscf: Option<f64>,             // 最终能量泛函的 HFX 系数（仅两步泛函）
    pub spin_channel: usize,
    pub reference: Vec<String>,
}
```

同时有 `xc_scf` 与 `xc_nscf` 时为两步（双杂化等）泛函；一步泛函只有前者非 `None`。

### 构造（解析 / 定义构造期）

| 方法 | 说明 |
|---|---|
| `DFAdef::new()` | 空定义 |
| `parse_and_derive(xc, spin, print)` | 解析入口：xc 字符串 -> `DFAdef` |
| `from_dfa4rest(&DFA4REST)` | legacy/DeepLearning 路线的适配器（code-array -> 组分；生产路径的唯一转换方向） |
| `to_dfa4rest()` | 转为旧版 DFA4REST 表示（全量）：仅 round-trip parity 测试使用，生产路径无 consumer |
| `to_xc_data()` | 仅 SCF 侧的 DFA4REST：仅测试使用，计划移除 |
| `update_pt2_params(os, ss)` | 覆写 PT2 `(c_os, c_ss)` 因子（ctrl 关键词；DL-DFA 迭代中亦会调用） |
| `get_hybrid_scf(spin)` / `get_hybrid_nscf(spin)` | 从组分推导 SCF / 最终能量泛函的总 HFX 系数（parse 时写入存储 field） |
| `get_rsh_scf(spin)` / `get_rsh_scf_lrsr(spin)` | 从组分推导 RSH 参数（原始 `(omega, c_LR, c_SR - c_LR)` / 消费者约定 `(omega, c_LR, c_SR)`；parse 时写入存储 field） |
| `check_sanity()` | 定义自检 |

### accessor（消费者读取）

SCF/grad/hessian/analdrv 等消费者通过以下方法读取定义，不做任何修改。

| 方法 | 返回 | 说明 |
|---|---|---|
| `is_dfa_scf()` | `bool` | SCF 泛函是否含 libxc 组分（在格点上计算 XC） |
| `is_hybrid()` / `is_rsh()` | `bool` | 是否杂化 / 是否 range-separated |
| `rsh_params_scf()` | `Option<(f64, f64, f64)>` | SCF 泛函 RSH 参数 `(omega, c_LR, c_SR)`（parse 时存储；`omega()` 读取单个参数） |
| `family_pos()` | `Option<DFAFamily>` | post-SCF 家族（PT2/RPA/SCSRPA/SBGE2），由 `FIFTH_RUNG_COMPONENT_TYPES` 支撑 |
| `fifth_rung_component()` / `is_fifth_dfa()` / `is_rpa()` | — | 第五阶组分与其谓词 |
| `is_nlc()` / `vv10_params()` | — | 非局域相关（VV10）判定与 `(b, C)` 参数 |
| `x_part(spin)` / `c_part(spin)` | `Option<DFAdef>` | 纯交换 / 纯相关部分（按 libxc functional kind 分解；组合泛函返回 `None`） |
| `has_nscf()` / `use_density_gradient()` / `use_kinetic_density()` / `use_laplacian()` | `bool` | 定义形状谓词 |
| `summary()` / `formatted_output()` | — | 输出 |

### 定义 -> 计算列表（边界构造）

定义与计算的边界：`(factor, LibXCFunctional)` 列表在此构造（libxc 泛函初始化一次，含 `_omega` 参数设置），之后供 exc/vxc 等 free function 使用。

| 方法 | 说明 |
|---|---|
| `scf_func_list(spin)` | SCF 泛函的计算列表 |
| `nscf_func_list(spin)` | 最终能量泛函的 DFT 部分列表（纯 MP2 家族为 `None`） |
| `func_list_from_components(comps, spin)` | 任意组分列表的 door（response/grad 接口使用） |

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

### 构造 / 归一化（解析期）

| 方法 | 说明 |
|---|---|
| `DFAComponent::new(factor, func)` / `from_libxc(id, factor)` | 构造组分 |
| `canonicalize()` / `normalize_rsh()` / `normalize_pt2_param()` | 解析期归一化（RSH 约定、PT2 参数形状） |

### accessor（消费者读取）

| 方法 | 说明 |
|---|---|
| `to_libxc_functional(spin)` | 初始化该组分的 `LibXCFunctional`（非 Libxc 组分返回 `None`） |
| `get_hybrid(spin)` / `get_rsh(spin)` | 单组分的 HFX / RSH 参数 |
| `is_rsh()` / `is_nlc()` | 组分谓词 |
| `get_dftd_params()` | Disp 组分的参数（dftd3-rs / dftd4-rs toml） |
| `formatted_output()` | 组分信息输出 |

非 libxc 组分类型与输入名称的对应（`ComponentType`、alias、参数）见 [parse_xc](parse_xc.md)。

## 计算 API（free function）

定义与计算解耦：计算函数位于 `dft`（`mod.rs` 等）下，接收 `(factor, LibXCFunctional)` 列表，而不接收定义对象或 `(factor, id)`。

| 函数 | 说明 |
|---|---|
| `xc_exc_vxc(&funcs, grids, spin, dm, mo, occ, print_level)` | XC 能量 + 势 |
| `xc_exc_vxc_slots(&funcs, range_grids, ...)` / `_dm_only` | 分块（screening）版本 |
| `xc_exc(&funcs, grids, ...)` | 仅 XC 能量 |
| `post_xc_exc(&names, ...)` / `post_tabulated_exc(&funcs, ...)` | post-SCF 泛函的能量 |
| `DFAdef::scf_func_list(spin)` / `nscf_func_list(spin)` / `func_list_from_components` | 定义 -> 计算列表的构造入口 |
| `xc_part_of_name(name, spin, kind)` / `DFAdef::x_part` / `c_part` | 名称/定义 -> 纯 X 或纯 C 部分 |

### XCData（导数求值包）

`XCData`（`src/dft/num_int.rs`）是梯度/Hessian 侧的求值包：持有 func 列表 + `xc_type` + 密度来源（`dm` 或 `mo_coeffs`/`occ`）。产生链路：`mol.dfadef` -> `scf_func_list(spin)` / `func_list_from_components` -> 经 `NumInt::gen_xc_data`（`grad/rks.rs`、`grad/uks.rs`）或 `xc_skeleton_dao`（`grad/rdh.rs`、`grad/udh.rs`）构造 `XCData` -> `get_vxc_rayon_new` / `eval_xc_eff`。

注意：`DFAdef::to_xc_data()`（旧版 DFA4REST 转换，计划移除）与 `XCData` 无关——同名不同物。

## 与旧版 DFA4REST API 的对照

`DFA4REST`（现位于 `src/dft/legacy.rs`）是旧版 code-array 表示；legacy/DeepLearning 路线解析后经 `DFAdef::from_dfa4rest` 转为 DFAdef。状态标记：→ 等价替换；≈ 替换但语义不完全相同；✕ 移除；L 仅存于 legacy；N 新增。

| 旧（DFA4REST） | 新（DFAdef / free fn） | 状态 |
| --- | --- | --- |
| `dfa_compnt_scf: Vec<usize>` + `dfa_paramr_scf: Vec<f64>` | `xc_scf: Option<Vec<DFAComponent>>` | ≈ |
| `dfa_compnt_pos` + `dfa_paramr_pos` | `xc_nscf: Option<Vec<DFAComponent>>` | ≈（pos 数组只含 Libxc code；xc_nscf 含全部非 SCF 组分） |
| `dfa_paramr_adv: Option<Vec<f64>>`（os/ss） | `fifth_rung_component().param_positional` | ≈ |
| `dfa_hybrid_scf` / `dfa_hybrid_nscf`(旧 `dfa_hybrid_pos`) | 同名 field + `get_hybrid_scf/get_hybrid_nscf` | → |
| `dfa_rsh_scf`（消费者约定三元组） | `dfa_rsh_scf`（原始）+ `dfa_rsh_scf_lrsr`（消费者约定）+ `rsh_params_scf()` | ≈ |
| `dfa_family_pos`（field） | `family_pos()`（accessor，`FIFTH_RUNG_COMPONENT_TYPES`） | ≈ |
| `nlc_vv10: Option<(b, c)>` | `is_nlc()` / `vv10_params()` | → |
| `is_dfa_scf()/is_hybrid()/is_rsh()/is_fifth_dfa()/is_rpa()` | 同名方法 | → |
| `use_density_gradient()/use_kinetic_density()` | 同名方法 | → |
| `omega()` | `omega()`（field 读取） | → |
| `rsh_params()`（重复推导） | `rsh_params_scf()`（parse 时存储） | → |
| `rsh_alpha()` / `use_eri()` | - | ✕ |
| `xc_exc_vxc(&self, ...)` 等 exc/vxc 方法 | `dft::xc_exc_vxc(&funcs, ...)` 等 free fn | ≈（funcs 货币） |
| `libxc_func_list()` | `DFAdef::scf_func_list(spin)` | → |
| `get_hybrid_libxc` / `get_rsh_libxc` | `DFAComponent::get_hybrid` / `get_rsh` | ≈（parity-tested） |
| `libxc_code_fdqc` / `xc_func_init_fdqc` | `xc_part_of_name` / `x_part` / `c_part` | ≈ |
| `summary(print_level)` / `describe()` / `formatted_output()` | `summary()` / `formatted_output()` | ≈（组分视角） |
| `new(xc, spin, print)`（legacy parser） | `parse_and_derive(xc, spin, print)` | ≈ |
| `new_nonstandard(...)` | - | ✕（NonStandard 输入已移除） |
| - | `x_part` / `c_part` / `xc_part_of_name` / `from_libxc` | N |

## 相关

- [ADR-0003：DFA 定义统一到 DFAdef](../adr/adr-0003-dfa-definitions.md)
- [parse_xc：解析语法](parse_xc.md)
- [libxc 接口](libxc.md)
