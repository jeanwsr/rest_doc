# nimatmul DFT 格点积分与导数

nimatmul (numerical integration by matrix multiplication) 是基于矩阵乘法的 DFT 格点积分模块，不使用 einsum 或其他张量缩并工具。REST 中对应的实现位于 `rest/src/dft/numint_matmul` 与 `rest/src/dft/gen_grids/becke_partitioning_deriv`。

这组文档说明 nimatmul 的 API 设计概念、程序接口、设计决定，以及求解 DFT Hessian (能量二阶与 Fock 矩阵一阶 Skeleton 导数) 的公式约定与实现策略，包括格点偏移 (grid-shift) 导数的贡献。前一部分文档移植自 rstsr-showcase-dft-grids 原型仓库，Hessian 相关文档移植自 rstsr-showcase-hessian 原型仓库，并均已依 REST 中的正式实现作校正。

:::{note}
除特别注明外，各文档均只处理闭壳层问题；开壳层需要对自旋分量作相应展开 (参见各文档开头的小节)。张量维度采用 column-major 约定；`becke-deriv` 一文依其程序实现风格采用 row-major 约定，见该文开头的说明。
:::

:::{note}
`concept`、`basic-api`、`adr` 三篇文档由 AI (Claude Code + glm-5.3) 从 rstsr-showcase-dft-grids 的设计文档 (cn-dft-api-concepts) 改写并拆分而来，公式与接口均依 REST 当前实现作过事实核查；原文的记号已统一为 `def` 一文的约定。
:::

```{toctree}
:maxdepth: 1
def
concept
basic-api
adr
skeleton2
vmat1
becke-deriv
becke-grid-shift
```
