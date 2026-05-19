# RI 的 2c-2e ERI 分解


## 3c-2e 分解后的 ERI `cderi`

RI 算法经常需要使用下述四指标算子：

$$
g_{\mu \nu, \kappa \lambda} \approx g_{\mu \nu, P} (\mathbf{J}^{-1})_{PQ} g_{\kappa \lambda, Q}
$$

从以下两个角度出发：

- 三张量计算比较复杂，张量越少越好
- 整个表达式对称性很高

我们可能更希望将上述表达式作分解：

$$
\boxed{
Y_{\mu \nu, P} = \sum_Q (\mathbf{J}^{-1/2})_{PQ} g_{\mu \nu, Q}
}
$$

该计算过程是 $O(n_\text{basis}^2 n_\text{aux}^2)$ 即 $O(N^4)$ 标度算法；但相比于 4c-2e ERI，其计算量更小 (与矩阵乘法相同)。这样之后，四角标算子可以简化为两个张量的乘法：

$$
g_{\mu \nu, \kappa \lambda} \approx \sum_P Y_{\mu \nu, P} Y_{\kappa \lambda, P}
$$

该分解在程序中通常记为 `cderi`；其名称是 Cholesky Decomposed ERI 的缩写，但也可能有其他分解策略。这点我们后面会提到。

并不是所有算法都一定要采用上述分解。具体来说：
- SCF 计算通常需要分解，特别是需要 RI-K 的情景。在自洽场迭代过程中 (除开分解步) 之外，RI-K 的计算量是 $O(n_\text{occ} n_\text{basis}^2 n_\text{aux})$；但如果不分解，每次迭代都要增加至少 $O(n_\text{occ} n_\text{basis} n_\text{aux}^2)$ 的计算量，而这在基组没有大到 4-5 zeta 时计算量要比 $O(n_\text{basis}^2 n_\text{aux}^2)$ 更大，而且程序实现起来也不方便。
- MP2 计算通常习惯作分解。分解可以简化程序实现；是否作分解对 MP2 的计算量没有影响，且分解这一步的 $O(N^4)$ 标度要小于 MP2 的 $O(N^5)$ 标度。
- RPA 方法 / LT 算法可以不作分解。RPA 算法标度是 $O(n_\text{occ} n_\text{vir} n_\text{aux}^2)$ 的，是比分解过程相当或更小的。

## 2c-2e ERI 分解：Cholesky 分解方法

需要留意，这里的 $(\mathbf{J}^{-1/2})_{PQ}$ 不一定必须是严格的矩阵平方根；我们也可以使用其他分解方法，例如 Cholesky 分解，来得到一个三角矩阵 $L_{PR}$。请留意，上下三角的处理方式是不同的：

- 下三角 (row-major)：$J_{PQ} = \sum_R L_{PR} L_{QR}$ 或 $\mathbf{J} = \mathbf{L} \mathbf{L}^\dagger$
- 上三角 (col-major)：$J_{PQ} = \sum_R L_{RP} L_{RQ}$ 或 $\mathbf{J} = \mathbf{L}^\dagger \mathbf{L}$

总之，$(\mathbf{J}^{-1/2})_{PQ}$ 的目的仅仅是将矩阵 $\mathbf{J}$ 分解到两边，使得 $g_{\mu \nu, \kappa \lambda} \approx \sum_P Y_{\mu \nu, P} Y_{\kappa \lambda, P}$ 表达式成立；至于怎样分解，做法并不是唯一的。

在使用 Cholesky 分解方法时，必须要注意，Cholesky 分解的矩阵并不是对称的，因此转置与否是非常重要的问题。这与狭义的 $(\mathbf{J}^{-1/2})_{PQ}$ 可以任意转置是不同的。

我们的程序取 col-major，因此 $g_{\mu \nu, P}$ 的指标顺序是 $(\mu, \nu, P)$ 即辅助基指标 $P$ 是内存最不连续的。那么 Cholesky 分解的做法是

$$
Y_{\mu \nu, P} = \sum_Q g_{\mu \nu, Q} (\mathbf{L}^{-1})_{PQ}
$$

注意三角矩阵求逆后乘法，在 BLAS 中有 DTRSM (在 RSTSR 中有 solve_triangular) 实现。这种情形是 `side = Left`, `trans = NoTrans`。

采用 Cholesky 分解方案的好处是
- Cholesky 分解与 DTRSM 的计算量都非常小；
- DTRSM 可以 inplace 计算，不需要额外内存。

Cholesky 分解方案的坏处是
- 这要求 2c-2e ERI $J_{PQ}$ 具有比较好的正定性 (本征值不能太小)；否则在处理矩阵求解的 $(\mathbf{L}^{-1})_{PQ}$ 时可能会有数值不稳定的问题，从而导致失败。
- 也因此，Cholesky 最好不要作为唯一的分解方案。

同时兼顾 Cholesky 分解好处与坏处的方案，可能会是一种比较特殊的 pivotted Cholesky 与矩形三角矩阵求解，但目前我们暂时不考虑用这种方案。

## 2c-2e ERI 分解：本征值分解方法

本征值分解方法就是依矩阵线性代数的严格定义，求取 $(\mathbf{J}^{-1/2})_{PQ}$。对于厄米对称的矩阵，由于本征值的存在，这个定义是容易实现的：记本征向量 $\mathbf{U}$ 与本征值 $\mathbf{\Lambda}$，则

$$
\mathbf{J}^{-1/2} = \mathbf{U} \mathbf{\Lambda}^{-1/2} \mathbf{U}^\dagger
$$

但对于正定性不强的情景，$\mathbf{\Lambda}$ 中可能会有比较小的本征值，导致 $\mathbf{\Lambda}^{-1/2}$ 中有比较大的元素，从而导致数值不稳定的问题。我们本意上求 $\mathbf{J}^{-1/2}$ 的目的不是得到它本身，而是为了求得近似的 $g_{\mu \nu, \kappa \lambda}$。因此，通过设置一个阈值来丢弃掉比较小的本征值，是合理且不会严重损失精度的做法。


## 相关的函数


**函数 `get_solved_j3c`**


函数路径：`ri_jk::pure_decompose::get_solved_j3c`


该函数将输入 $\mathbf{g}$ 并返回分解后的结果 $\mathbf{Y}$。

也需要指出，输入 $\mathbf{g}$ 的维度不需要是 2-dim 或 3-dim，它可以是任意维度的。对于该张量，
- 我们**要求** $\mathbf{g}$ 的最后一个维度 (最不连续的维度) 是待分解的辅助基维度。
- 我们**非常建议** $\mathbf{g}$ 的其他维度是 col-major 连续的，因为我们会使用 reshape 函数来将 $\mathbf{g}$ 视作一个二维矩阵来进行分解计算；如果其他维度不连续，reshape 将产生额外内存开销。

具体的计算过程上，

- 对于 Cholesky 上三角情形，$\mathbf{Y} = \mathbf{g} \mathbf{L}^{-1} = \big( (\mathbf{L}^\mathrm{T})^{-1} \mathbf{g}^\mathrm{T} \big)^\mathrm{T}$；
- 对于 Cholesky 下三角情形，$\mathbf{Y} = \mathbf{g} (\mathbf{L}^\mathrm{T})^{-1} = \big( (\mathbf{L}^{-1}) \mathbf{g}^\mathrm{T} \big)^\mathrm{T}$；
- 对于本征值分解情形，$\mathbf{Y} = \mathbf{g} \mathbf{J}^{-1/2}$。

Cholesky 分解之所以要作复杂的转置，是因为 RSTSR 中的 `rt::linalg::solve_triangular((a, b, uplo))` 被求解的问题是 $\mathbf{A}^{-1} \mathbf{b}$，左右顺序无法调整；而我们需要求解的问题是 $\mathbf{b} \mathbf{A}^{-1}$，因此需要转置来调整左右顺序。

