# TDDFT 激发态解析梯度

模块路径：`ri_tddft::tddft_grad`（`src/ri_tddft/tddft_grad.rs`，入口 `TddftGradEngine`）

## 概述

`TddftGradEngine` 计算单个激发态激发能的解析核梯度。实现是 PySCF `pyscf/grad/tdrks.py`（`_contract_xc_kernel`、`grad_elec`）与 `pyscf/grad/tdrhf.py`（`grad_elec`、CPHF Z-vector、Pulay 项）到 REST 的 RI 表示的忠实移植。梯度组装为

$$
\mathbf{F}_{\mathrm{TDDFT}} = \mathbf{F}_{\mathrm{GS}} + \bigl[\mathbf{F}_{\mathrm{resp}}(x, y) - \mathbf{F}_{\mathrm{resp}}(\mathbf{0}, \mathbf{0})\bigr]
$$

其中 $\mathbf{F}_{\mathrm{GS}}$ 是 REST 已验证的基态 RKS 梯度，方括号内是振幅依赖的响应部分（`assemble(x,y) - assemble(0,0)`，后者抵消了过渡密度中的基态成分）。

## 理论背景

激发能为 $\omega = \sum_{ia} (x+ y)_{ia} (A_{ia,jb} + B_{ia,jb}) (x - y)_{jb} / 2$ 形式的 Rayleigh 商（TDA 时 $y = 0$），其梯度由三部分贡献：

1. **电子响应梯度**：过渡密度 ($\mathbf{P} = C (x+y) C^{\mathrm{T}}$ 与差分密度) 的 0 阶 J/K 势、导数 J/K 积分表、XC 核 (fxc) 的两遍格点收缩（第二遍需要 Z-vector）；
2. **CPHF Z-vector**：弛豫的密度响应通过解 Z-vector 方程获得（`solve_zvector`，复用 `ri_cphf::CPHFSolverPySCF` 的 CPHF 机制）；
3. **Pulay 项与核梯度**：RI 辅助基与格点的位置依赖带来的附加项（`aux_forces_atom`、生成函数导数、重叠导数）。

PySCF 中两次 `get_jk` 调用的分工在 RI 近似下均被复现：

- 0 阶 `get_jk(mol, dm, hermi=0)`：普通 J/K AO 矩阵（`j_potential`/`k_potential`，用于 `veff0` 系列）；
- 导数 `get_jk`：积分导数表 $v^{J}_{t,\mu\nu} = -\sum_{ls} \frac{\partial}{\partial R_\mu} (\mu\nu|ls) D_{ls}$（`vj_bra_atom`/`vk_bra_atom` + `AtomDerivBlocks`，用于 `veff1`）；辅助基与度量的响应由双线性形式从完整势导数中加入（`RawRiTensors::d_i_from_blocks` + `d_j_atom`，复用 `ri_gw::gw_grad` 的基础设施）。

## 使用方式

```toml
[tddft]
     tddft_method =              "tda"    # 或 "lr"
     nroots =                    3
     tddft_grad_state =          1        # 1-based；0 关闭
```

- 仅支持**限制性参考**（`spin_polarization = true` 时 panic）；
- `tddft_grad_state = N` 表示对第 N 个激发态求梯度，须不超过 `nroots`，且同一任务中必须先完成 TDDFT 求解；
- `main_driver::eval_force` 在基态梯度之上叠加所选态的响应贡献后输出总梯度；几何优化/MD 等调用同一梯度接口，因此可直接对激发态做结构优化；
- 通道选择：`tddft_spin = "triplet"` 时取三重态通道，其余情况（含 `"both"`）取单重态通道；
- RSH 泛函：交换贡献带短程修正 (引擎内 `hyb` = $c_{LR}$，`hyb_sr` = $c_{SR} - c_{LR}$)。

## 关键数据结构

| 成员/结构 | 意义 |
|--|--|
| `TddftGradEngine` | 单激发态梯度引擎：`scf`、`state`、`singlet`、`tda`、振幅 `(x, y)`、RI 原始张量 `raw`、预求逆的 RI 度量 `jinv`、交换系数 (`hyb`, `hyb_sr`)、fxc 缓存 `fxc_cache`、libcint 句柄 `cint` |
| `RawRiTensors` / `AtomDerivBlocks` | 复用自 `ri_gw::gw_grad`：三中心积分的原子导数分块 |
| `FxcHessianCache` / `VindWorkspace` | 复用自 `dft::response`：fxc 二阶核的格点缓存与求值工作区 |
| `AOMat` / `Density` / `RankFactor` | 列优先 AO 矩阵与显式秩因子分解 $D = \sum l\, r^{\mathrm{T}}$（交换收缩需要因子结构） |
| `AoBlocks` | 整格点 AO 表缓存（可选；见下文内存调优） |

## 计算流程

```text
TddftGradEngine::new(scf, state, singlet, tda, x, y)
    ├── 构建原始 RI 张量 (build_raw_ri_tensors) + 度量逆 jinv (一次)
    └── response_gradient()
        └── assemble(x, y) - assemble(0, 0)
            ├── 过渡密度: P = x+y / x-y 组合 → dvv, doo, P_pl
            ├── 0 阶势: vj / vk (kf ≠ 0 时才算 K) → veff0 系列
            ├── XC 核第一遍: contract_xc_kernel (fxc0)
            ├── Z-vector: solve_zvector(wvo) (CPHF, ri_cphf)
            ├── XC 核第二遍: contract_xc_kernel (fxcz1, 需 Z-vector)
            ├── 导数势: vj_bra_atom / vk_bra_atom (AtomDerivBlocks)
            ├── 生成函数导数 hcore、重叠导数、aux_forces_atom (Pulay)
            └── 组装 [3, natm]
```

`contract_xc_kernel` 每次梯度调用两遍（第二遍的 `fxcz1` 需要 Z-vector），两遍走同一格点、同一 `ao_deriv`。`eval_ao_batch` 约占响应成本的 30%，第二遍是纯重复计算：是否缓存取决于内存预算（见下）。

## 性能与内存调优 (环境变量)

| 环境变量 | 作用 |
|--|--|
| `REST_TDDFT_GRAD_AOCACHE_MB` | 覆盖整格点 AO 表缓存的内存预算（仅当 AO 张量能装下时才缓存，否则逐批重算） |
| `REST_TDDFT_GRAD_NO_AOCACHE` | 禁用 AO 表缓存 |
| `REST_TDDFT_GRAD_XCBLK_MB` | 覆盖 XC 格点分块大小（缺省按 `max_memory` 预算） |
| `REST_TDDFT_GRAD_AO_LEGACY=1` | 恢复重构前的 `eval_ao_batch` 求值路径（仅回归交叉验证用；libcint 路径快 5-6 倍） |
| `REST_TDDFT_GRAD_TIME` / `REST_TDDFT_GRAD_MEM` | stderr 分步计时 / RSS 峰值内存追踪 |

## 限制

- 仅限制性 (RHF/RKS) 参考；非限制梯度未实现；
- 当前为激发能梯度；激发态之间的态交叉 (conical intersection) 场景未做特殊处理；
- 梯度为 RI 近似下的解析梯度（辅助基响应已含）。
