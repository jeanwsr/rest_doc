# ADR-0003：DFA 定义

- **状态**：使用中
- **日期**：2026-11-03
- **适用范围**：crate `rest` 的 DFA 定义（`src/dft/parse_xc/`）与全部 DFA 消费者（`scf_io`、`grad`、`hessian`、`analdrv`、`ri_tddft`、`ri_gw`、`ri_pt2` 等）。存量代码中 `DFA4REST`（已移至 `src/dft/legacy.rs`）与 `Molecule.xc_data`（现为 `Option`）尚存，供调试使用，将随 legacy parser 退役。
- **作者**：王石嵘

## 决定

- **原则一：定义与计算解耦**。完成全部的定义工作后，再转入计算。计算函数中不做定义相关的事情。例如，计算函数接收 `(factor, LibXCFunctional)`，而不是 `(factor, id)`。
- **原则二：一切皆组分**。DFA 的每个部分都要被定义为组分，包括 D3/NLC/HF。
  - 组分的参数不添加为整体的 field，而是从组分（`DFAComponent`）的 api/field 获取。（不过考虑到和旧版 api 的一致性，以及本身实现的方便，`DFAdef` 还是有一些整体的 field，例如整体的 hybrid/rsh 参数。）
  - DFA 整体的参数推荐从 accessor fn 获取，例如 `rsh_params_scf()`。
- `DFAdef`（`src/dft/parse_xc/`）是 DFA 的主定义。
  - `DFA4REST` 移入 `src/dft/legacy.rs`，legacy 与 DeepLearning 路线解析后经 `DFAdef::from_dfa4rest` 转为 DFAdef；
  - `xc_exc_vxc` 之类的计算函数成为 free function。

决定的核心模式（完整的 field/方法表见 [DFA 定义 API](../dft/dfa-definitions.md)）：

```rust
let dfa = parse_and_derive("xyg3", 1, 0);

// 组分访问：遍历 SCF / 最终能量部分的组分列表
for comp in dfa.xc_scf.as_ref().unwrap() {
    // 组分的 field：factor, func, id, param_positional, param_keyword, component_type, reference
    println!("{} {} {}", comp.factor, comp.func, comp.id);
}

// 第五阶（post-SCF）组分
let pt2 = dfa.fifth_rung_component().unwrap();
assert_eq!(pt2.component_type, ComponentType::PT2);
assert_eq!(pt2.param_positional, [0.3211, 0.3211]); // (c_os, c_ss)

// 泛函级参数：从 accessor 获取，而不是从组分重新推导
let (omega, c_lr, c_sr) = dfa.rsh_params_scf().unwrap();
let hyb = dfa.dfa_hybrid_scf;
```

## 1. DFAComponent 与 ComponentType

DFA 的每个部分是一个组分 `DFAComponent`：

| field | 说明 |
|---|---|
| `factor` | 系数 |
| `func` / `func_full_name` | 组分名 / libxc 全名 |
| `id` | libxc 编号（Libxc 组分） |
| `param_positional` / `param_keyword` | 参数（如 PT2 的 `(c_os, c_ss)`、RSHF 的 `omega`） |
| `component_type` | 组分类型 |
| `reference` | 文献 |

`ComponentType` 取值：`Libxc`（libxc 泛函）、`HF`（精确交换）、`RSHF`（range-separated HF：`SR_HF`/`LR_HF`/`RSH`）、`PT2`（MP2 型相关）、`RPA` / `SCSRPA` / `SBGE2`（第五阶 RPA 家族）、`Disp`（经验色散 D3/D4/VV10）、`Unknown`。

完整的 field/方法表见 [DFA 定义 API](../dft/dfa-definitions.md)。

## 2. 背景

本决定的最初目的是将 `parse_xc` 设为默认 parser，以避免各种需要手动设置 parser 的麻烦。由于 `DFA4REST` 与 `DFAdef` 实际上都是描述 DFA 定义的结构体（从功能上来说前者并不只是定义，也包括计算函数，但其本身的 field 只有定义），因此将 DFA 定义相关的 API 统一到 `DFAdef` 中。`Molecule` 现在同时持有 `dfadef`（必有）与 `xc_data: Option<DFA4REST>`（仅 legacy/DeepLearning 路线为 `Some`）；DFA 的全部消费者通过 `mol.dfadef` 的 accessor（`is_dfa_scf`、`is_hybrid`、`is_rsh`、`family_pos`、`fifth_rung_component`、`rsh_params_scf`、`is_nlc`、`vv10_params` 等）读取定义。

## 3. 备选方案

旧版 API 并未删除：`DFA4REST` 的全部函数与 `Molecule.xc_data` field 仍然存在，只是正常调用路径上已无 consumer。因此在调试时暂时恢复旧版 API 是可行的；但恢复旧版 API 调试后，希望实现尽量写在新版 API 的框架下。

## 4. 备注

具体的 API 细节仍然可以变动；本 ADR 只约束上述设计原则。新/旧 API 的对照见 [libxc 接口](../dft/libxc.md) 与 [parse_xc](../dft/parse_xc.md)。
