# SCF 波函数稳定性分析

函数路径：`ri_tddft::stability::stability`（`src/ri_tddft/stability.rs`）

## 概述

SCF 收敛解可能只是能量表面的**鞍点**而非极小点。稳定性分析在收敛解附近对能量作二阶展开，检查轨道旋转空间的 Hessian 是否正定：存在负本征值即说明存在能量下降的轨道旋转方向，该 SCF 解不稳定。

REST 的稳定性分析完全构建在 AO 模式 TDDFT 的 (A/B) 矩阵-矢量积机制之上（见 [tddft-ao](tddft-ao.md)）：稳定性算符就是 **(A+B) 轨道 Hessian**，本征值用批量 Davidson 求解器取最低若干个。分析为**只查 (check-only)** 模式：输出最低若干本征值与稳定/不稳定结论，不输出本征矢量，也不自动旋转轨道并重跑 SCF (orbital following 未实现)。

## 理论背景

稳定性算符是占据→虚轨道旋转下的能量二阶导数 (Seeger & Pople, JCP 66, 3045 (1977))，包含 B 型 (去激发/对称化密度) 响应，通过 TDDFT 的 A、B 算符在同一参考上表达：

$$
\mathbf{H}_{\mathrm{stab}} = f \cdot (\mathbf{A}^{x} + \mathbf{B}^{x})
$$

REST 采用的通道与因子 (与 PySCF `rhf_internal`/`uhf_internal`/`rhf_external` 一致)：

| 检查 | Hessian | 通道 `$x$` | 因子 `$f$` |
|--|--|--|--|
| RHF/RKS internal | $4(\mathbf{A}^{\mathrm{S}} + \mathbf{B}^{\mathrm{S}})$ | `'S'` (单重) | 4 |
| RHF/RKS external (RHF→UHF) | $\mathbf{A}^{\mathrm{T}} + \mathbf{B}^{\mathrm{T}}$ | `'T'` (三重) | 1 |
| UHF/UKS internal | $2(\mathbf{A} + \mathbf{B})$ (拼接 $[\alpha;\beta]$ 旋转空间) | `'R'` | 2 |

**稳定性判据**：$\lambda_{\min} \ge -10^{-5}$ (`STABILITY_THRESHOLD`，PySCF 阈值) 则稳定；否则该通道不稳定（存在能量下降的轨道旋转方向）。

参考文献：Seeger & Pople, JCP 66, 3045 (1977)；Bauernschmitt & Ahlrichs, JCP 104, 9047 (1996)；PySCF `scf/stability.py`。

## 使用方式

```toml
[tddft]
     stability =                 "auto"    # off | internal | external | full | auto
     stability_nroots =          3         # 求解的最低本征值数目
     stability_tol =             1.0e-8    # Davidson 收敛阈值
```

- `stability` 取值：

| 取值 | 行为 |
|--|--|
| `"off"` (默认) | 不做稳定性分析 |
| `"internal"` | RHF/RKS 单重或 UHF/UKS 轨道 Hessian |
| `"external"` | RHF→UHF (三重通道) 检查，仅 RHF/RKS；UHF/UKS 下打印提示并跳过 (UHF→GHF 未实现) |
| `"full"` | internal + external |
| `"auto"` (推荐) | 自动选择适用于当前参考类型的全部检查；当前解析为 `"full"` |

- 顶层 `[ctrl]` 区块的 `check_stab` 接受相同取值，作为 `stability = "off"` 时的回退；两者同时给出时 `[tddft] stability` 优先。仅做稳定性分析的输入卡无需 `[tddft]` 之外的激发态设置 (`stability_nroots`/`stability_tol` 有结构默认值)。
- 稳定性分析与激发态计算**互斥**：`stability != "off"` 时任务只做稳定性检查 (在 `main_driver` 中、SCF 收敛后、激发态计算之前执行)。
- 结果逐检查打印本征值与结论，并写入 `rest_results.json` 的 `"stability"` 字段 (`internal`/`external` 各自的 roots 与 stable 标志)，供回归测试机器可读地比对。

## 实现细节

```text
stability(scf, mode)
    ├── 解析 nroots/tol（[tddft] 缺省时取结构默认；"auto" → "full"）
    ├── ROHF 参考 → 报错；DFT 参考无格点 → 报错（HF 参考无需格点）
    ├── internal:
    │   ├── (factor, xlet) = RHF: (4, 'S') / UHF: (2, 'R')
    │   ├── prepare_ao_data_with_spin(scf, Some("singlet") | None)
    │   └── hessian_roots(): batched Davidson on factor×(A+B)
    ├── external (仅 RHF/RKS):
    │   ├── prepare_ao_data_with_spin(scf, Some("triplet"))
    │   └── hessian_roots(): (factor, xlet) = (1, 'T')
    └── StabilityReport { stable_internal, roots_internal,
                          stable_external, roots_external }
```

- **Hessian-矢量积**：`a_matvec_ao_batched` + `b_matvec_ao_batched` 之和乘因子 (`xlet` 进入矩阵-矢量积)。对角预条件器取 `factor × (ε_a − ε_i)`。
- **Davidson 配置**：`max_subspace = 60`、`max_iter = 200`、`add_dim = nroots + 2`、`tol = stability_tol` (REST 收敛约定 $\|r\| < \sqrt{\varepsilon}$，故 $10^{-8}$ 对应残差 $10^{-4}$，与 PySCF `STAB_TOL` 一致)。
- **HF 参考**：`prepare_ao_data_with_spin` 对无 libxc 组分的参考返回 `fxc_driver: None`，Hessian 只含 RI J/K 部分，因此无需 DFT 格点即可运行。
- **每检查独立准备数据**：internal (单重核) 与 external (三重核) 各自调用 `prepare_ao_data_with_spin`，不共享自旋适配核。

## 限制与未实现

- 实→复稳定性 ($^1(\mathrm{A}'-\mathrm{B}')$ 检查) 未实现；
- UHF→GHF 外稳定性未实现 (UHF/UKS 只做 internal)；
- ROHF 参考不支持；
- 轨道跟随 (检测不稳定后旋转轨道并重跑 SCF) 未实现，当前为只查模式。
