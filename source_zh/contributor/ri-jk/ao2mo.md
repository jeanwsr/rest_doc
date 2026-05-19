# 双指标原子-分子轨道转换


我们将处理下述类型的双指标原子-分子轨道转换：

$$
g_{i a, P} = \sum_{\mu \nu} C_{\mu i} C_{\nu a} g_{\mu \nu, P}
$$

这里的程序主要是为 post-SCF RI-incore 算法作准备的；但 $g_{\mu \nu, P}$ 未必一定要是三维张量；指标 $P$ 也同样可以指代 $\kappa \lambda$ 等双指标。


## 算法说明


上述转换的输入是三角存储的 $g_{\mathrm{tp}(\mu\nu), P}$，需要先将其展开为 $g_{\mu\nu, P}$，然后再与轨道系数作两次矩阵乘法。由于 $i$ 与 $a$ 属于不同的轨道空间，输出 $g_{ia, P}$ 关于 $(i, a)$ 没有对称性。

在实现上，我们支持多组 (set) 的 bra 与 ket 轨道系数；每组可以有不同的 $n_i$ 与 $n_a$。程序会根据 $n_i$ 与 $n_a$ 的大小关系，选择先缩并较小维度的指标，以减少临时存储的大小。

目前提供两种输出指标顺序的版本：
- **notrans**：输出指标顺序为 $(i, a, P)$，即 $P$ 在最后；适用于需要按分子轨道 $(i, a)$ 对访问的情形。post-SCF 方法经常采用这种模式。
- **trans**：输出指标顺序为 $(P, i, a)$，即 $P$ 在最前；适用于需要按辅助基 $P$ 访问的情形。

此外，函数支持输入与输出的数据浮点类型不同 (通过 `cast` 函数转换，一般是类似于 `|x| T::from_f64(x).unwrap()` 等简单闭包转换)。需要留意当输入输出类型相同时，会使用 `unsafe` 的类型转换来避免额外的内存分配与类型转换开销，无论用户输入了哪种函数。


## 具体算法与函数接口


**函数 `get_ao2mo_s2ij_to_s1_notrans_with_output`**


函数路径：`ri_jk::pure_ao2mo::get_ao2mo_s2ij_to_s1_notrans_with_output`


输出指标顺序为 $(i, a, P)$。对辅助基指标 $P$ 作并行循环。

当 $n_i < n_a$ 时 (先缩并 bra)：

$$
\begin{aligned}
g_{\mu\nu, P} &\bowtie g_{\mathrm{tp}(\mu\nu), P}
&& \text{(eq.1)} \\
S_{\nu i, P}^\mathbb{A} &= \sum_\mu C_{\mu i}^\mathbb{A} \, g_{\mu\nu, P}
&& \text{(eq.2)} \\
g_{ia, P}^\mathbb{A} &= \sum_\nu S_{\nu i, P}^\mathbb{A} \, C_{\nu a}^\mathbb{A}
&& \text{(eq.3)}
\end{aligned}
$$

当 $n_a \le n_i$ 时 (先缩并 ket)：

$$
\begin{aligned}
g_{\mu\nu, P} &\bowtie g_{\mathrm{tp}(\mu\nu), P}
&& \text{(eq.1)} \\
S_{\mu a, P}^\mathbb{A} &= \sum_\nu g_{\mu\nu, P} \, C_{\nu a}^\mathbb{A}
&& \text{(eq.2)} \\
g_{ia, P}^\mathbb{A} &= \sum_\mu C_{\mu i}^\mathbb{A} \, S_{\mu a, P}^\mathbb{A}
&& \text{(eq.3)}
\end{aligned}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `j3c` | $g_{\mathrm{tp}(\mu\nu), P}$ | $(\mathrm{tp}(\mu\nu), \ldots)$ | $(n_\mathrm{tp(basis)}, \ldots)$ | 至少 1 维；$P$ 可为多指标 |
| `bra` | $C_{\mu i}^\mathbb{A}$ | $(\mu, i)$ | $(n_\mathrm{basis}, n_i^\mathbb{A})$ | 各 set 的 $n_i$ 可不同 |
| `ket` | $C_{\nu a}^\mathbb{A}$ | $(\nu, a)$ | $(n_\mathrm{basis}, n_a^\mathbb{A})$ | 各 set 的 $n_a$ 可不同 |
| `output` | $g_{ia, P}^\mathbb{A}$ | $(i, a, \ldots)$ | $(n_i^\mathbb{A}, n_a^\mathbb{A}, \ldots)$ | 其余维度与 `j3c` 一致 |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| thread | (eq.1) | $g_{\mu\nu, P}$ | $(\mu, \nu)$ | $(n_\mathrm{basis}, n_\mathrm{basis})$ | 对每个 $P$ 并行 |
| thread | (eq.2) | $S$ | $(n_\mathrm{basis}, n_s)$ | $(n_\mathrm{basis}, n_s)$ | $n_s = \min(n_i, n_a)$ |
| thread | (eq.3) | $S_\mathrm{out}$ (类型不同时) | $(n_i, n_a)$ | $(n_i, n_a)$ | 仅当输入输出类型不同时 |


**函数 `get_ao2mo_s2ij_to_s1_trans_with_output`**


函数路径：`ri_jk::pure_ao2mo::get_ao2mo_s2ij_to_s1_trans_with_output`


输出指标顺序为 $(P, i, a)$，即 $P$ 在最前。为尽可能使用 BLAS3 函数以提升性能，我们对 $P$ 指标进行分批处理。默认分批大小为 $\max(\lceil 0.2 \, n_P \rceil, 16)$。也因此，用户最好要在使用该函数前，保证有额外 20% 的内存可用来存储中间结果；或者直接指定一个较小的 `nbatch` 来降低内存需求。

当 $n_i < n_a$ 时 (先缩并 bra)：

$$
\begin{aligned}
g_{\mu\nu, \underline{P}} &\bowtie g_{\mathrm{tp}(\mu\nu), \underline{P}}
&& \text{(eq.1)} \\
S_{\nu \underline{P}, i}^\mathbb{A} &= \sum_\mu g_{\mu\nu, \underline{P}} \, C_{\mu i}^\mathbb{A}
&& \text{(eq.2)} \\
g_{\underline{P}, ia}^\mathbb{A} &= \sum_\nu S_{\nu \underline{P}, i}^\mathbb{A} \, C_{\nu a}^\mathbb{A}
&& \text{(eq.3)}
\end{aligned}
$$

当 $n_a \le n_i$ 时 (先缩并 ket)：

$$
\begin{aligned}
g_{\mu\nu, \underline{P}} &\bowtie g_{\mathrm{tp}(\mu\nu), \underline{P}}
&& \text{(eq.1)} \\
S_{\mu \underline{P}, a}^\mathbb{A} &= \sum_\nu g_{\mu\nu, \underline{P}} \, C_{\nu a}^\mathbb{A}
&& \text{(eq.2)} \\
g_{\underline{P}, ia}^\mathbb{A} &= \sum_\mu C_{\mu i}^\mathbb{A} \, S_{\mu \underline{P}, a}^\mathbb{A}
&& \text{(eq.3)}
\end{aligned}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `j3c` | $g_{\mathrm{tp}(\mu\nu), P}$ | $(\mathrm{tp}(\mu\nu), \ldots)$ | $(n_\mathrm{tp(basis)}, \ldots)$ | 至少 1 维；$P$ 可为多指标 |
| `bra` | $C_{\mu i}^\mathbb{A}$ | $(\mu, i)$ | $(n_\mathrm{basis}, n_i^\mathbb{A})$ | 各 set 的 $n_i$ 可不同 |
| `ket` | $C_{\nu a}^\mathbb{A}$ | $(\nu, a)$ | $(n_\mathrm{basis}, n_a^\mathbb{A})$ | 各 set 的 $n_a$ 可不同 |
| `nbatch` | | | | 默认 $\max(\lceil 0.2 \, n_P \rceil, 16)$ |
| `output` | $g_{P, ia}^\mathbb{A}$ | $(\ldots, i, a)$ | $(\ldots, n_i^\mathbb{A}, n_a^\mathbb{A})$ | 前置维度与 `j3c` 一致 |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| batched | (eq.1) | $g_{\mu\nu, \underline{P}}$ | $(\mu, \nu, \underline{P})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{batch})$ | |
| batched | (eq.2) | $S$ | $(n_\mathrm{basis}, n_\mathrm{batch}, n_s)$ | $(n_\mathrm{basis}, n_\mathrm{batch}, n_s)$ | $n_s = \max(n_i, n_a)$ |
| thread | (eq.3) | $S_\mathrm{out}$ (类型不同时) | $(n_\mathrm{batch}, n_s)$ | $(n_\mathrm{batch}, n_s)$ | 仅当输入输出类型不同时 |
