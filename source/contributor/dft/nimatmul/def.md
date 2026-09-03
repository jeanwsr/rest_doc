# DFT Grid Integration Conventions

The main purpose of this document is to describe the DFT grid integration conventions used in the program. This includes formula notation, common program variables, shape information, and so on.

This document only describes the conventions of nimatmul, a simple DFT grid integration module.

## 1. DFT Energy

This section explains the core goal of DFT grid integration in the quickest way. Even though nimatmul is a very naive DFT grid integration implementation, and the formulas and the program implementation involve a large amount of tedious details, all of the work amounts to nothing more than the computation of the DFT energy and its derivatives.

DFT grid integration is the process of numerically integrating over space under the functional kernel $f(\rho, \gamma, \tau, \cdots)$:

$$
E^\text{xc} = \int f \cdot \rho (\boldsymbol{r}) \, \mathrm{d} \boldsymbol{r}
$$

Typical DFT programs are grid-discretized. Grid discretization is an approximation, but we always state it with an equals sign.

$$
E^\text{xc} = \sum_g w_g f_g \rho_g
$$

where
- $g$ is the grid index,
- $\rho_g$ is the electron density $\rho(\boldsymbol{r}_g)$ at grid point $g$;
- $f_g$ is the functional kernel $f(\rho_g, \gamma_g, \tau_g, \cdots)$ at grid point $g$;
- $w_g$ is the weight at grid point $g$.

:::{note}
Additional convention: We do not discuss Laplacian-type functionals for now. Therefore, in the end our basic DFT variables only include $\rho$, $\gamma$, $\tau$, and not $\nabla^2 \rho$.

Additional remark: Many programs use sigma to denote the GGA $\gamma = \nabla \rho \cdot \nabla \rho$. Since we use that notation as the index of the spin components, we do not use sigma to denote the GGA gradient density.
:::

## 2. Basic DFT Variables $\xi$

For convenience in later program implementation and formula derivation, we write the basic DFT variables in the vector form $\boldsymbol{\xi}$[^note-1]:

[^note-1]: This notation is my personal habit. The notation is taken from Su, N. Q.; Zhang, I. Y.; Xu, X. *J. Comput. Chem.* **2013**, *34* (20), 1759–1774. doi: [10.1002/jcc.23312](https://doi.org/10.1002/jcc.23312). For the $\xi_\chi$ of this document, the corresponding notation in that paper is $\zeta_\eta$. The reason for changing $\zeta$ to $\xi$ is that $\zeta$ is often used in DFT to denote the spin polarization $\zeta = (\rho^\alpha - \rho^\beta) / \rho$, and is also used in basis set functions to denote the exponentially decaying parameter $\zeta$. The reason for changing $\eta$ to $\chi$ is simply that we may later use `x` for the index of this variable in einsum notation.

$$
\boldsymbol{\xi} = (\rho, \rho^x, \rho^y, \rho^z, \tau)
$$

The component superscript of the basic DFT variables $\boldsymbol{\xi}$ is denoted $\chi$. Note that $\boldsymbol{\xi}$ is not the same for LDA, GGA, and mGGA. LDA has only the first component ($\chi$ has length 1); GGA has only the first four components ($\chi$ has length 4); mGGA contains all five components ($\chi$ has length 5).

$$
\xi_{g \mu \nu}^\chi = \begin{cases}
\phi_{g \mu} \phi_{g \nu} & \chi = \rho \\
\phi_{g \mu}^r \phi_{g \nu} + \phi_{g \mu} \phi_{g \nu}^r & \chi = \rho^r \\
\sum_r \frac{1}{2} \phi_{g \mu}^r \phi_{g \nu}^r & \chi = \tau
\end{cases} 
$$

Here, we use $r \in \{x, y, z\}$ to denote the spatial dimension index. The spatial dimension is also often denoted by $t, s$. In the future, we will try to use the index $r$ for the spatial dimension index in density expressions involving GGA or mGGA, while using the indices $t, s$ in gradient computations, in order to distinguish the two.

It should be pointed out that GGA is traditionally expressed with $\gamma = \nabla \rho \cdot \nabla \rho$ as the variable; in fact, the DFT functional engines LibXC and XCFun are designed this way in their programs. However, since $\gamma$ is a second-order quantity in the density matrix, it brings a considerable degree of complexity to program implementation and formula derivation. In future work, we will always convert $\gamma$ to the form of $\rho_x, \rho_y, \rho_z$ for processing.

Introducing the density matrix $D_{\mu \nu}$, we can obtain the spatial grid representation of the basic DFT variables:

$$
\xi_g^\chi = \sum_{\mu \nu} D_{\mu \nu} \xi_{g \mu \nu}^\chi
$$

Note that, in quite a few situations, using the orbital coefficients $C_{\mu i}$ in place of the density matrix $D_{\mu \nu}$ is often the more efficient approach in the program. Here we only give the definition and will not elaborate.

## 3. DFT Functional Kernel Derivatives

We define the functional kernel derivatives

$$
\begin{aligned}
f^\chi &= \frac{\partial (f \rho)}{\partial \xi^{\chi}} \\
f^{\chi \chi'} &= \frac{\partial^2 (f \rho)}{\partial \xi^{\chi} \partial \xi^{\chi'}} \\
f^{\chi \chi' \chi''} &= \frac{\partial^3 (f \rho)}{\partial \xi^{\chi} \partial \xi^{\chi'} \partial \xi^{\chi''}}
\end{aligned}
$$

Note that the quantity being differentiated here is not the functional kernel $f$ itself, but its product with the density, $f \rho$.

## 4. Program Conventions

**Shape**

| Program index | Formula index | Shape | Meaning |
|--|--|--|--|
| `u, v` | subscript $\mu, \nu$ | `nao` $n_\mathrm{basis}$ | basis function |
| `g` | subscript $g$ | `ngrids` $n_\mathrm{grids}$ | DFT grid |
| `t, s, r` | superscript $t, s, r \in \{x, y, z\}$ | 3 | spatial components |
| `x, y` | superscript $\chi, \chi'$ | `nvar` $n_\mathrm{var}$ | basic DFT variable components |
| `A, B` | superscript $A, B$ | `natm` $n_\mathrm{atom}$ | atoms |
| `A, B` | superscript $\mathbb{A}, \mathbb{B}$ | `nprop` $n_\mathrm{prop}$ | arbitrary property |
| `i, j` | subscript $i, j$ | `nocc` $n_\mathrm{occ}$ | occupied orbitals |
| | superscript $*$ | `ncomp` | atomic-orbital derivative components |
| `σ, ς` | superscript $\sigma, \varsigma \in \{ \alpha, \beta \}$ | 2 | spin components |

**Tensors (closed-shell)**

| Variable | Formula | Shape (PySCF/mixed-major) | Shape (REST/col-major) |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(*, g, \mu)$ <br> `[ncomp, ngrids, nao]` | $(g, \mu, *)$ <br> `[ngrids, nao, ncomp]` |
| `rho` | $\xi_{g}^{\chi}$ | $(\chi, g)$ <br> `[nvar, ngrids]` | $(g, \chi)$ <br> `[ngrids, nvar]` |
| `vxc` | $f_g^\chi$ | $(\chi, g)$ <br> `[nvar, ngrids]` | $(g, \chi)$ <br> `[ngrids, nvar]` |
| `fxc` | $f_g^{\chi \chi'}$ | $(\chi, \chi', g)$ <br> `[nvar, nvar, ngrids]` | $(g, \chi, \chi')$ <br> `[ngrids, nvar, nvar]` |

The star $*$ appearing here generically refers to the atomic-orbital derivative components; its exact number depends on the highest order of derivative involved (counting the zeroth-order component, there are 4 in total up to first-order derivatives, 10 up to second order, and 20 up to third order). It is not the complex conjugate.

Note that, although the NumPy used by PySCF follows the row-major convention, for $\phi_{g \mu}^{*}$ there, even though its apparent shape is $(*, g, \mu)$, its memory contiguity has $g$, i.e. the grid index, in the innermost position, followed by $\mu$, i.e. the basis function, and finally the atomic-orbital derivative components. This memory-contiguity order is in fact the same as in the column-major REST. For the other variables, REST (column-major) and PySCF/NumPy (row-major) have opposite shapes but the same memory contiguity.

**Tensors (open-shell)**

| Variable | Formula | Shape (PySCF/mixed-major) | Shape (REST/col-major) |
|--|--|--|--|
| `rho` | $\xi_{g}^{\sigma \chi}$ | $(\sigma, \chi, g)$ <br> `[2, nvar, ngrids]` | $(g, \chi, \sigma)$ <br> `[ngrids, nvar, 2]` |
| `vxc` | $f_g^{\sigma \chi}$ | $(\sigma, \chi, g)$ <br> `[2, nvar, ngrids]` | $(g, \chi, \sigma)$ <br> `[ngrids, nvar, 2]` |
| `fxc` | $f_g^{\sigma \chi \sigma' \chi'}$ | $(\sigma, \chi, \sigma', \chi', g)$ <br> `[2, nvar, 2, nvar, ngrids]` | $(g, \chi, \sigma, \chi', \sigma')$ <br> `[ngrids, nvar, 2, nvar, 2]` |

Note that the `fxc` tensor possesses the $(\sigma, \sigma')$ and $(\chi, \chi')$ symmetries. For convenience in programming, we give up the symmetry and use the full-tensor representation.
