# LibXC/XCFun 记号约定变换

:::{note}
**该文档有 AI 参与翻译**

中文版由 AI (Claude Code + glm-5.3-flash) 自英文原版翻译；代码、表格与函数名保持原文。
:::

:::{note}
本文移植自原型仓库 [rstsr-showcase-dft-grids](https://github.com/RESTGroup/rstsr-showcase-dft-grids)；同一文档也以 rustdoc 形式内嵌于 `rest/src/dft/xceff/xc_deriv.rs` 的函数 `libxc_to_xcfun_mapping` 上 (代码内的副本不含本页末尾的 Python 参考实现)。

在 REST 中，该变换在 `libxc_eval_eff` 内部施加于 LibXC 原始输出，先于 $\sigma \to \rho_r$ 的展开；周边 API 见 [nimatmul/basic-api](nimatmul/basic-api.md)。
:::

## 由具体例子说明

LibXC 与 XCFun 在处理自旋极化情形时的约定有所不同。

我们以 GGA 自旋极化的二阶导数 (fxc) 为例进行说明。

首先指出，对能量部分 (zk) 与一阶导数部分 (vxc)，LibXC 与 XCFun 的约定相同。其顺序为

| 索引 | 记号 |
|--|--|
| 0 | zk |
| 1 | r_u |
| 2 | r_d |
| 3 | s_uu |
| 4 | s_ud |
| 5 | s_dd |

对上表的记号：

- Rho 是自旋极化的两分量变量：`r_u` 表示 $\rho^\uparrow$，`r_d` 表示 $\rho^\downarrow$。
- Sigma 是自旋极化的三分量变量：`s_uu` 表示 $\sigma^{\uparrow \uparrow}$，`s_ud` 表示 $\sigma^{\uparrow \downarrow}$，`s_dd` 表示 $\sigma^{\downarrow \downarrow}$。

对 fxc 部分，涉及两个密度变量。LibXC 与 XCFun 的根本区别在于

- LibXC 先按密度类型分类 (rho 在 sigma 之前)，再按自旋排序；
- XCFun 先按自旋分量分类 (依 r_u, r_d, s_uu, s_ud, s_dd 的顺序)。

对未来的 DFT 计算而言，XCFun 风格更为有利。但 LibXC 支持的泛函更多、API 设计更好、使用也更广泛。因此需要一个索引变换映射。对本例而言，

| 索引 | 记号<br>LibXC | 记号<br>XCFun | Map |
|-:|--|--|-:|
|  6 | r_u  / r_u  | r_u  / r_u  |  6 |
|  7 | r_u  / r_d  | r_u  / r_d  |  7 |
|  8 | r_d  / r_d  | r_u  / s_uu |  9 |
|  9 | r_u  / s_uu | r_u  / s_ud | 10 |
| 10 | r_u  / s_ud | r_u  / s_dd | 11 |
| 11 | r_u  / s_dd | r_d  / r_d  |  8 |
| 12 | r_d  / s_uu | r_d  / s_uu | 12 |
| 13 | r_d  / s_ud | r_d  / s_ud | 13 |
| 14 | r_d  / s_dd | r_d  / s_dd | 14 |
| 15 | s_uu / s_uu | s_uu / s_uu | 15 |
| 16 | s_uu / s_ud | s_uu / s_ud | 16 |
| 17 | s_uu / s_dd | s_uu / s_dd | 17 |
| 18 | s_ud / s_ud | s_ud / s_ud | 18 |
| 19 | s_ud / s_dd | s_ud / s_dd | 19 |
| 20 | s_dd / s_dd | s_dd / s_dd | 20 |

该映射读作 `out[k] = in[map[k]]`：XCFun 第 $k$ 个槽位上的元素取 LibXC 索引 `map[k]` 处的值。

我们可以看到

- LibXC 先分为 `r/r` (6--8)、`r/s` (9--14)、`s/s` (15--20) 三类；每类内部按自旋排序。
- XCFun 先分为 `r_u` (6--10)、`r_d` (11--14)、`s_uu` (15--17)、`s_ud` (18--19)、`s_dd` (20) 五类；每类内部，第二个变量与第一类中的排序方式一致。

## 例子的一般化

- **更高阶导数**：我们可能遇到更高阶的导数 (通常至四阶，但可能更高)。
- **更多种类的密度**：我们可能使用更多密度类型。其优先级为 RHO > SIGMA > TAU > LAPL。
  - TAU：`t_u`、`t_d`
  - LAPL：`l_u`、`l_d`
- 我们假定 RHO (LDA) 只输入 RHO；SIGMA (GGA) 输入 RHO 与 SIGMA；TAU (部分 meta-GGA) 输入 RHO、SIGMA、TAU；LAPL (部分 meta-GGA) 输入 RHO、SIGMA、TAU、LAPL 全部 (尽管部分 LAPL 型 meta-GGA 实际并不输入 tau，为简单起见我们仍要求 TAU 可用)。

## REST 中的实现

该变换实现于 `rest/src/dft/xceff/xc_deriv.rs`：

- `libxc_to_xcfun_mapping_parts(den_type, deriv)` 保存硬编码的映射表，目前覆盖 SIGMA 与 TAU 至四阶导数。依赖 LAPL 的泛函、以及超过四阶的导数尚未支持。
- `libxc_to_xcfun_mapping(den_type, spin, deriv)` 在 `deriv <= 1`、自旋非极化、或密度类型为 RHO 时返回 `None` (恒等变换，不重排)；这些恰是两种约定一致的场合。否则返回 `0..=deriv` 各阶映射的拼接，与 LibXC 输出的扁平布局相对应。
- `libxc_transform_xcfun_indices` 通过 `index_select` 将映射施加于 LibXC 输出的最后一维 (索引维)。它在 `libxc_eval_eff` 中被调用，先于 $\sigma \to \rho_r$ 的展开 (`transform_xc_inner`)。

## 自动生成索引映射的代码

```python
from itertools import combinations_with_replacement
from math import comb


def libxc_to_xcfun_indices_map(den_type: str, deriv: int) -> list[int]:
    """Spin-Polarized Indices Map from LibXC to XCFun.

    Parameters
    ----------
    den_type : str
        Density Type. Supports `rho`, `sigma`, `tau`, `lapl`.
    deriv : int
        Derivative level.

    Example
    -------
    >>> libxc_to_xcfun_indices_map("sigma", 2)
    [6, 7, 9, 10, 11, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    """
    # Each variable: (type_priority, spin_index)
    # RHO has 2 spin components, SIGMA has 3, TAU has 2, LAPL has 2
    group_specs = [
        ("rho", 0, 2),
        ("sigma", 1, 3),
        ("tau", 2, 2),
        ("lapl", 3, 2),
    ]

    type_map = {
        "rho": ["rho"],
        "sigma": ["rho", "sigma"],
        "tau": ["rho", "sigma", "tau"],
        "lapl": ["rho", "sigma", "tau", "lapl"],
    }

    if den_type not in type_map:
        raise ValueError(f"Unknown den_type: {den_type}")

    active_groups = set(type_map[den_type])

    # Build variable list: each variable is (type_priority, spin_index)
    variables = []
    for group_name, priority, n_spin in group_specs:
        if group_name in active_groups:
            for spin in range(n_spin):
                variables.append((priority, spin))

    d = len(variables)

    # Generate all non-decreasing multi-indices of length deriv
    # combinations_with_replacement yields them in lexicographic order = XCFun order
    xcfun_order = list(combinations_with_replacement(range(d), deriv))

    # LibXC order: sort by density type signature first, then by variable indices
    def libxc_key(mi):
        return tuple(variables[i][0] for i in mi) + mi

    libxc_order = sorted(xcfun_order, key=libxc_key)

    # Build reverse lookup: multi-index -> LibXC position
    libxc_pos = {mi: pos for pos, mi in enumerate(libxc_order)}

    # Base offset = sum of outputs for all previous derivative levels
    base_offset = sum(comb(d + i - 1, i) for i in range(deriv))

    return [base_offset + libxc_pos[mi] for mi in xcfun_order]
```
