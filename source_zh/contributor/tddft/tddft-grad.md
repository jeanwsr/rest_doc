# TDDFT 激发态解析梯度

`TddftGradEngine`（`src/ri_tddft/tddft_grad.rs`）计算单个激发态激发能的解析核梯度，由输入关键字 `tddft_grad_state` 触发；用户侧用法见 [TDDFT 激发态计算](../../user/tddft.md)。

## 实现要点

- **对标 PySCF**：梯度装配流程是对 `pyscf/grad/tdrks.py`（`_contract_xc_kernel`、`grad_elec`）与 `pyscf/grad/tdrhf.py`（`grad_elec`、CPHF Z-vector、Pulay 项）的移植，改写为 REST 的 RI 表示。
- **总体分解**：`de_TDDFT = de_GS + response`。`de_GS` 复用 REST 已验证的基态 RKS 梯度；`response = assemble(x, y) - assemble(0, 0)` 是依赖于激发振幅的部分。
- **两类 J/K 调用**（与 PySCF 一致）：零阶 J/K 为普通 AO 矩阵（对应 `veff0doo`/`veff0mop`/`veff0mom`）；导数表 `vj[t,μ,ν] = -Σ_ls ∂(μν|ls)/∂R D_ls`、`vk[t,μ,ν] = -Σ_ls ∂(μl|νs)/∂R D_ls`（对应 `veff1`）。两者均在 RI 近似下复现。
- **辅助基/度规响应**：由完整势导数 `RawRiTensors::d_i_from_blocks` + `d_j_atom` 按双线性形式逐项加入。

> TODO: 本页尚待补充更详细的公式推导与数值验证说明。
