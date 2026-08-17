# 求解器模块 `solvers`

> REST 的统一迭代求解器模块，位于 `src/solvers/`，提供 Davidson 本征值求解器与 Krylov 线性方程组求解器。该模块在 [PR#183](https://gitee.com/restgroup/rest/pulls/183) 中引入。

## 模块结构

```
src/solvers/
├── mod.rs          # 模块声明
├── davidson.rs     # Davidson 求解器 (davidson_solver / lr_davidson_solver)
├── krylov.rs       # Krylov 求解器 (krylov_tsr / krylov_vec)
└── unify.md        # 求解器设计文档 (架构、特性、测试覆盖)
```

---

## Davidson 求解器

### API

| 函数 | 用途 |
|---|---|
| `davidson_solver()` | TDA 本征值问题 `A x = e x`（别名 `tda_davidson_solver`） |
| `lr_davidson_solver()` | 对称化 Casida 方程 `(A-B)^(1/2)(A+B)(A-B)^(1/2) z = omega^2 z` |
| `generate_initial_guess()` | 由 `diag` 最小元生成初始猜测向量 |

两个求解器均接受 `FnMut(&Vec<f64>) -> Vec<f64>` 形式的矩阵向量积闭包，返回 `Vec<(f64, Vec<f64>)>`（本征值、本征向量）。

### 收敛判据

$$
\|\mathbf{r}\| < \sqrt{\mathrm{tol}} \quad \wedge \quad |\Delta e| < \mathrm{tol}
$$

其中 $\mathbf{r} = A\mathbf{x} - e\mathbf{x}$ 为真实残差，$\Delta e$ 为本征值相邻两次迭代的变化。

### `DavidsonConfig` 默认值

| 字段 | 默认值 | 说明 |
|---|---|---|
| `max_subspace` | `60` | 子空间最大维度 |
| `add_dim` | `2` | 每轮迭代新增的向量数 |
| `restart_dim` | `6` | 重启后保留的子空间维度 |
| `max_iter` | `100` | 最大迭代次数 |
| `tol` | `1e-10` | 收敛阈值（见上方收敛判据） |
| `lindep` | `1e-14` | 线性相关阈值，`\|\|v\|\|^2 < lindep` 的向量被丢弃 |
| `use_mgs` | `true` | Modified Gram-Schmidt 正交化（`false` 时用 CGS） |
| `divergence_restart` | `true` | 检测发散（`\|\|r\|\|/\|\|r_last\|\| > 3`）后恢复并重启 |
| `track_states` | `false` | 通过重叠矩阵在迭代间重排本征态顺序 |

### 特性

- 可配置收敛阈值 `tol`
- 预处理器分母夹取：`|denom|` 下限为 `1e-8`，避免除零爆炸
- MGS / CGS 正交化可切换
- 发散检测与自动重启
- 本征态跟踪（近简并根）
- 线性相关向量过滤（`lindep`）

### 日志输出

- 输出等级由 `log` 框架控制，由输入卡 `print_level` 关键词映射（见 [logger.md](logger.md)）
- 每轮迭代输出两行：`residues: [...]` / `|de|: [...]`，以及 `iter N: space=, converged=` 汇总

---

## Krylov 求解器

### API

| 函数 | 用途 |
|---|---|
| `krylov_tsr()` | Tsr/TsrView 核心，完整算法 |
| `krylov_vec()` | `Vec<f64>` 便捷包装（`krylov` 为其别名） |

`solve (I + A) x = b`。返回 `(Tsr, KrylovResult)`，其中 `KrylovResult` 记录 `cycles`、`aop_calls`、`per_root_solves`、`residual` 等统计信息。

### 收敛判据

子空间收敛信号为：

$$
\max(\|\mathbf{v}\|^2) < \max(\mathrm{lindep},\ \mathrm{tol}^2)
$$

当子空间坍缩（QR 分解后无新方向）但真残差仍较大时，触发 per-RHS 恢复机制：

1. **per-RHS 独立子空间**：每个 RHS 维护独立的 Krylov 基（matvec 仍批处理）
2. **递归 per-root 求解**：仍未收敛的 RHS 以独立的单 RHS `krylov_tsr` 递归求解（新子空间、无共享耦合）
3. **容忍阈值**：真残差 `||r|| < max_residual_factor * tol` 的 RHS 直接接受，不触发递归求解

### `KrylovConfig` 默认值

| 字段 | 默认值 | 说明 |
|---|---|---|
| `tol` | `1e-9` | 收敛阈值 |
| `max_cycle` | `50` | 最大迭代轮数 |
| `max_space` | `None` | 硬重启子空间上限（`None` 不重启） |
| `lindep` | `1e-15` | 线性相关阈值 |
| `max_residual_factor` | `1000.0` | 容忍系数，`\|\|r\|\| < factor * tol` 视为可接受 |

---

## 调用方与用户关键字映射

| 求解器 | 调用方 | 用户关键字 | 默认值 |
|---|---|---|---|
| Davidson | BSE | `davidson_converge_threshold` | `1e-10` |
| Davidson | TD-DFT | `davidson_tol` | `1e-10` |
| Davidson | SCF 稳定性 (`addons.rs`) | — | `tol = 1e-5` |
| Krylov | `analdrv` | `cphf_tol` | `1e-9` |
| Krylov | `analdrv` | `cphf_lindep` | `1e-15` |
| Krylov | `analdrv` | `cphf_tol_inflation` | `1000.0` |
| Krylov | Hessian (`[ctrl.hessian]`) | `krylov_tol` | `1e-9` |
| Krylov | Hessian | `krylov_lindep` | `1e-15` |
| Krylov | Hessian | `krylov_tol_inflation` | `1000.0` |

> 注：Davidson 的 `tol` 同时控制两个收敛判据——真残差 `||r|| < sqrt(tol)` 与本征值变化 `|de| < tol`。Krylov 的 `tol` 则控制子空间收敛信号阈值 `max(||v||^2) < max(lindep, tol^2)`。两者语义不同，请勿混淆。

### 代码入口

- BSE：`src/ri_bse/mod.rs`（4 处 `davidson_solver` / `lr_davidson_solver` 调用）
- TD-DFT：`src/ri_tddft/tddft_solver.rs`
- SCF 稳定性：`src/scf_io/addons.rs`
- CP-HF（Hessian）：`src/ri_cphf/cphf_solver_pyscf.rs` → `src/hessian/rhf.rs`
- analdrv：`src/analdrv/krylov_block.rs`（`krylov_tsr` 的薄包装）
