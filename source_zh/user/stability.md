# SCF 波函数稳定性分析

SCF 收敛解可能只是能量表面的**鞍点**而非极小点。REST 提供 SCF 波函数稳定性分析（Seeger & Pople, JCP 66, 3045 (1977)）：在 SCF 收敛后，检查能量对占据→虚轨道旋转的二阶导数（轨道 Hessian）是否正定。最低本征值 $\lambda_{\min} < -10^{-5}$ 即报告不稳定。

REST 支持的检查：

- **internal (内稳定性)**：RHF/RKS 单重通道或 UHF/UKS 轨道 Hessian；
- **external (外稳定性)**：RHF→UHF 方向 (三重通道) 检查，仅 RHF/RKS 参考。

分析为**只查 (check-only)** 模式：输出最低若干本征值与稳定/不稳定结论，但不自动旋转轨道并重跑 SCF。

## 输入关键字

以下关键字在输入卡的 `[tddft]` 区块中声明：

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `stability` | String | `"off"` | 稳定性分析模式：`"internal"` (仅内稳定性)、`"external"` (仅外稳定性，仅 RHF/RKS)、`"full"` (internal + external)、`"auto"` (推荐，自动选择适用于当前参考类型的全部检查)。设为 `"off"` 以外的值时，本任务**只做稳定性分析、不做激发态计算** |
| `stability_nroots` | usize | 3 | 稳定性 Hessian 求解的最低本征值数目 |
| `stability_tol` | f64 | 1.0e-8 | 稳定性 Hessian 的 Davidson 收敛阈值 |

NOTE: 顶层 `[ctrl]` 区块中的 `check_stab` 关键字接受与 `stability` 相同的取值，作为其回退；两者同时存在时 `[tddft] stability` 优先。仅需稳定性分析时，输入卡无需准备激发态计算的其他设置（`stability_nroots`/`stability_tol` 有默认值，`[tddft]` 区块甚至可以只含 `stability` 一行）。

## 用法示例

```toml
[ctrl]
     xc =                        "b3lyp"
     basis_path =                "def2-tzvp"
     auxbas_path =               "def2-tzvp-rifit"
     eri_type =                  "ri-v"

[tddft]
     stability =                 "auto"
     stability_nroots =          3
```

程序在 SCF 收敛后执行稳定性分析，逐检查打印最低本征值与结论，例如：

```sh
RHF/RKS internal stability: lowest eigenvalues = [0.012, 0.045, 0.132]
RHF/RKS wavefunction is stable in the internal stability analysis
```

结果同时写入 `rest_results.json` 的 `"stability"` 字段 (`internal`/`external` 各自的 roots 与 stable 标志)，供脚本化后处理与回归测试比对。

## 结果解读

| 参考 | internal | external |
|--|--|--|
| RHF/RKS | $4(\mathbf{A}^{\mathrm{S}} + \mathbf{B}^{\mathrm{S}})$，检测单重通道失稳 | $\mathbf{A}^{\mathrm{T}} + \mathbf{B}^{\mathrm{T}}$，检测 RHF→UHF (三重通道) 失稳 |
| UHF/UKS | $2(\mathbf{A} + \mathbf{B})$ (拼接 α/β 旋转空间) | 未实现 (UHF→GHF)，自动跳过 |

判据：$\lambda_{\min} \ge -10^{-5}$ 则该通道稳定；出现负本征值说明 SCF 解是该通道下的鞍点，建议更换初始猜测 (`initial_guess`)、混合器或改用非限制参考重新计算。

## 注意事项

- **与激发态计算互斥**：`stability` 非 `"off"` 时本任务只做稳定性检查，不计算激发能。
- **参考类型**：ROHF 参考不支持；实→复稳定性与 UHF→GHF 外稳定性尚未实现。
- 推荐使用 `"auto"`：自动执行适用于当前参考的全部检查。

开发者细节 (理论因子、Davidson 配置、实现流程) 见 [SCF 波函数稳定性分析 (contributor)](../contributor/tddft/tddft-stability.md)。
