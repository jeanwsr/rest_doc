# RI 对 4c-2e ERI 近似原理简介


## RI 算法原理

RI 算法 (这里仅讨论与实现 RI-V) 的策略是将计算量较大的 4c-2e ERI，拆解成 3c-2e ERI 与 2c-2e ERI 的缩并计算。我们需要利用到
- 3c-2e ERI $g_{\mu \nu, P} = (\mu \nu | P)$
- 2c-2e ERI $J_{PQ} = (P | Q)$

对 4c-2e ERI 的近似是

$$
\boxed{
(\mu \nu | \kappa \lambda) \approx \sum_{PQ} g_{\mu \nu, P} (\mathbf{J}^{-1})_{PQ} g_{\kappa \lambda, Q}
}
$$

原理性的内容到这里就结束了。其他都是具体的算法与实现细节。


## RI-JK 问题

电子积分的计算目的总是实现 Fock 矩阵计算 (在这里是与密度矩阵 $D_{\mu \nu}$ 或轨道系数 $C_{\mu i}$ 的缩并)，因此仅仅讨论拆分电子积分是不行的，还要具体地考虑张量缩并。

RI 算法可以应用于 SCF 与 post-SCF 的计算。对于 SCF，RI-JK 得到 Fock 矩阵的整体过程非常简单：

- J 积分：(我们这里约定，辅助基指标的 $J_{PQ}$ 是电子积分，原子轨道指标的 $J_{\mu \nu} [\mathbf{D}]$ 是 Fock 矩阵的库伦贡献)

    $$
    J_{\mu \nu}[\mathbf{D}] = \sum_{PQ} \sum_{\kappa \lambda} g_{\mu \nu, P} (\mathbf{J}^{-1})_{PQ} g_{\kappa \lambda, Q} D_{\kappa \lambda}
    $$

- K 积分：

    $$
    K_{\mu \nu}[\mathbf{D}] = \sum_{PQ} \sum_{\kappa \lambda} g_{\mu \kappa, P} (\mathbf{J}^{-1})_{PQ} g_{\nu \lambda, Q} D_{\kappa \lambda}
    $$

## RI 算法的一些重要特征

- RI 的 3c-2e ERI 一般需要全部存储于内存 (内存标度 $O(n_\text{basis}^2 n_\text{aux})$)；否则，电子积分的计算代价将会比较大。可以不严格地认为 RI 相对于 4c-2e ERI 是一个内存换时间的算法。
- RI 是强耦合的问题，或者说 RI 算法依某一个指标，作不增加计算量的分批，经常是难以做到的。
- 对于 MP2/RPA 级别的 post-SCF 方法，可以将内存用量减少到 $O(n_\text{occ} n_\text{virt} n_\text{aux})$；这仍然是三次标度的，但在基组较大、或冻结轨道较多时，内存占用会大幅降低。
