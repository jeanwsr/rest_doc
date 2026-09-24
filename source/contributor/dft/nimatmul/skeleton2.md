# Second-Order Skeleton Derivatives of the DFT Energy

This document discusses implementation strategies for the second-order skeleton derivatives of DFT. We use only matrix multiplication strategies, and no other optimization tools.

:::{note}
**This document handles only the closed-shell case**

For open shells, part of the computational procedure needs to be reprocessed in terms of spin components. We will discuss this issue in other documents.
:::

## 1. Second-Order Skeleton Derivatives: Overview

We split the DFT second-order skeleton derivatives into three parts:
- `fxc`, the part related to $f^{\chi \chi'}$;
- `vxc_diag`, the diagonal part (double derivatives on a single atom);
- `vxc_off`, the off-diagonal part (one derivative on each of two atoms).

Let us elaborate a little here. First, recall the computational expression for the DFT energy:

$$
E^\text{xc} = \int f \rho \, \mathrm{d} \boldsymbol{r} = \sum_g w_g f_g \rho_g = \sum_g w_g (f \rho)_g
$$

The skeleton derivative is defined as the derivative obtained by varying the nuclear coordinates while the density matrix (orbital coefficients) remains unchanged. Skeleton derivatives can be obtained by numerical differencing, with the orbital coefficients fixed and the nuclear coordinates varied.

For a DFT task, its first-order skeleton derivative is (for now, the partial-derivative notation denotes skeleton derivatives; we only need to remember not to differentiate the orbital coefficients):

$$
\frac{\partial E^\text{xc}}{\partial \mathbb{A}} = \int \frac{\partial (f \rho)}{\partial \mathbb{A}} \, \mathrm{d} \boldsymbol{r} = \int \sum_\chi \frac{\partial (f \rho)}{\partial \xi^\chi} \frac{\partial \xi^\chi}{\partial \mathbb{A}} \, \mathrm{d} \boldsymbol{r}
$$

Written in the form of grid integration:

$$
\frac{\partial E^\text{xc}}{\partial \mathbb{A}} = \sum_{g \chi} w_g f_g^\chi \frac{\partial \xi_g^\chi}{\partial \mathbb{A}}
$$

Differentiating the formula above once more, using the product rule:

$$
\begin{aligned}
\frac{\partial^2 E^\text{xc}}{\partial \mathbb{A} \partial \mathbb{B}}
&= \sum_{g \chi} w_g \left( \frac{\partial f_g^\chi}{\partial \mathbb{B}} \frac{\partial \xi_g^\chi}{\partial \mathbb{A}} + f_g^\chi \frac{\partial^2 \xi_g^\chi}{\partial \mathbb{A} \partial \mathbb{B}} \right) \\
&= \sum_{g \chi \chi'} w_g f_g^{\chi \chi'} \frac{\partial \xi_g^{\chi'}}{\partial \mathbb{B}} \frac{\partial \xi_g^\chi}{\partial \mathbb{A}} + \sum_{g \chi} w_g f_g^\chi \frac{\partial^2 \xi_g^\chi}{\partial \mathbb{A} \partial \mathbb{B}}
\quad (\texttt{de\_xc})
\end{aligned}
$$

At this point, we can already split off the first term (the `fxc` contribution). The second term is the `vxc` contribution; depending on the concrete partial-derivative procedure, we split it into a diagonal and an off-diagonal part. This splitting is not trivial, and will be discussed in detail later.

:::{warning}
**Missing grid-shift derivative contributions**

The discussion above holds for ideal integration over the full spatial coordinates. In reality, in quite a few cases it does give a reasonable Hessian (referring to hybrid LDA/GGA functionals).

But as a numerical integration method, we have not taken into account the effect of nuclear-coordinate changes on the grid points $\boldsymbol{r}_g$ and the weights $w_g$. Hessian calculations for meta-GGA generally require introducing the grid-shift contributions. We will discuss this issue in other documents.
:::

:::{note}
**Column-major shape convention**

For tensors involved in DFT grid derivatives, generally speaking, the dimensions ordered from the most contiguous (the leftmost dimension under column-major) are

- the grid dimension `g`, of size `ngrids` $n_\mathrm{grids}$;
- the basis dimensions `u, v`, of size `nao` $n_\mathrm{basis}$;
- the variable dimensions `x, y`, of size `nvar` $n_\mathrm{var}$;
- the spatial dimensions `t, s`, of size 3; as a special case, the symmetry-packed dimension `ts`, of size 6;
- the basis-function grid-derivative dimension `*`, of size `ncomp`;
- the atomic dimensions `A, B`, of size `natm` $n_\mathrm{atom}$.
:::

## 2. A Common Technique for Nuclear-Coordinate Partial Derivatives of Electron Integrals

The discussion above applies to arbitrary properties. Now we specialize: we consider only the case $\mathbb{A} = A_t$; here $A_t$ denotes the three-dimensional spatial component $t$ of atom $A$. Normally this would be written as the vector $\boldsymbol{R}_A$ (or its three-dimensional components $R_{At}$), but for convenience of matching the program implementation, we sometimes simplify to the notation $A_t$.

A common technique is that, under an **atom-centered orbital basis**,

$$
\partial_{A_t} \phi_\mu = - \partial_t \phi_\mu \delta_{\mu \in A} = - \phi_\mu^t \delta_{\mu \in A}
$$

Roughly speaking, this exploits the fact that a basis function depends on the electronic coordinates and the nuclear coordinate only through the relative coordinate $\boldsymbol{r} - \boldsymbol{R}_A$; by the chain rule, the nuclear-coordinate partial derivative can be converted into an electronic-coordinate partial derivative, but the object being differentiated can only be the basis functions expanded on the specific atom, not all basis functions.

We usually do not call DFT grid integration an electron integral (electron integrals generally refer to analytic integrals of simple operators with 1-2 electrons and 2-4 centers); but this does not affect the fact that DFT is in essence a one-electron integral. The technique here applies not only to ordinary electron integrals, but equally so to DFT.

:::{note}
**This trick applies only to atom-centered orbital bases**

This trick applies only to atom-centered orbital bases. The atomic orbitals here are not limited to Gaussian basis functions; they can also be Slater basis functions, numerical basis functions, and so on, as well as their corresponding pseudopotential bases. As long as the basis is an atom-centered orbital basis, the nuclear-coordinate partial derivatives of the electron integrals can always be converted into electronic-coordinate partial derivatives.

But another major class of basis functions is plane waves. For plane waves, basis-function derivatives should generally be recast as operator manipulations in momentum space, for which there are other implementation strategies; moreover, the present technique does not apply to plane-wave basis functions at all.
:::

Note that, in terms of scaling, atoms, although few in number, still constitute a scale. In second-order gradient skeleton-derivative problems, a fairly effective strategy is to handle the electronic derivatives first, and to handle the nuclear derivatives only for the parts that cannot be processed further. This also varies from problem to problem: here, for `fxc` we handled the atomic derivatives first, whereas for `vxc` we handled the electronic derivatives first. This will show up in the implementation details below.

## 3. Implementation Details of the `fxc` Contribution

### 3.1 The final contraction of `fxc`

We first examine the final contraction of `fxc`. To begin with, its contribution consists of 4 terms in total:

$$
\frac{\partial^2 E^\text{xc}}{\partial A_t \partial B_s} \leftarrow \sum_{g \chi \chi'} w_g f_g^{\chi \chi'} (\partial_{A_t} \xi_g^\chi) (\partial_{B_s} \xi_g^{\chi'})
\quad (\texttt{de\_fxc})
$$

**Function `get_de_fxc`**

| Variable | Meaning | Index order | Shape | Notes |
|--|--|--|--|--|
| `wf` | $w_g f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` | $\mathrm{sym} (\chi, \chi')$ |
| `drho` | $\partial_{A_t} \xi_g^\chi$ | $(g, \chi, t, A)$<br>`[g, x, t, A]` | `[ngrids, nvar, 3, natm]` | |
| `de_fxc`<br>(output) | | $(t, s, A, B)$ <br>`[t, s, A, B]` | `[3, 3, natm, natm]` | $\mathrm{sym} (tA, sB)$ |

Its implementation can be given quite conveniently with einsum, in NumPy, or in RSTSR with tblis incorporated.

```python
# Please note the indices interchange of row-major / col-major
# This is numpy code, so row-major is used, different to col-major of above table
de_fxc = np.einsum("g, xyg, Atxg, Bsyg -> ABts", weights, fxc, drho, drho)
```

Doing it without einsum is also easy; one just needs to be a bit careful with broadcasting of tensor dimensions.

Note that $\partial_{A_t} \xi_g^\chi$ and $\partial_{B_s} \xi_g^{\chi'}$ are actually the same thing, merely written with different indices.

The dominant cost of this computation is $2 \times 3^2 n_\mathrm{atom}^2 n_\mathrm{var} n_\mathrm{grids}$ FLOPs; it is an $O(N^3)$-complexity step that is computationally rather cheap but a fairly large memory bottleneck.

The `wf` in the formula above is easy to obtain: $w_g$ is a basic grid parameter, and $f_g^{\chi \chi'}$ can be obtained immediately through the wrapped LibXC interface. The difficult part is the computation of `drho` $\partial_{A_t} \xi_g^\chi$.

### 3.2 Overview and implementation decisions for the density-on-grid first-order skeleton gradient `drho`

First recall the definition of the density on grids:

$$
\xi_g^\chi = \sum_{\mu \nu} D_{\mu \nu} \xi_{g \mu \nu}^\chi
$$

Taking its skeleton partial derivative with respect to the nuclear coordinates, $\partial_{A_t}$, differentiates only $\xi_{g \mu \nu}^\chi$. $\xi_{g \mu \nu}^\chi$ has different definitions for LDA/GGA/mGGA; but roughly, it always takes the form
1. The density matrix $D_{\mu \nu}$ first makes a matrix multiplication with a ket $\phi_{g \nu}^*$ and contracts the index $\nu$. The computational complexity is $O(n_\mathrm{basis}^2 n_\mathrm{grid})$, i.e. $O(N^3)$.
2. It is then element-wise scaled with a bra $\phi_{g \mu}^*$ and the index $\mu$ is contracted; but how to carry out the contraction must be handled separately for LDA (RHO) / GGA (SIGMA) / MGGA (TAU). This step makes many passes over memory, actually takes quite some time, and is the difficult part of the implementation; but its computational complexity is really only $O(n_\mathrm{basis} n_\mathrm{grid})$, i.e. $O(N^2)$.

Step 1 above is a simple matrix multiplication whose implementation is a single line, so no new function is introduced:

$$
\bar{\phi}_{g \mu}^* = \sum_{\nu} D_{\mu \nu} \phi_{g \nu}^*
\quad (\texttt{ao\_dm0})
$$

- `ao_dm0` $\bar{\phi}_{g \mu}^*$ is close to `ao` $\phi_{g \mu}^*$; it is a tensor of dimensions $(g, \mu, *)$ `[ngrids, nao, ncomp]`.
- In actual GGA/MGGA second-order gradient calculations, the superscript $*$ of `ao_dm0` $\bar{\phi}_{g \mu}^*$ goes at most to the 4 cases of first-order gradients (no derivative, $x$/$y$/$z$); whereas the superscript $*$ of `ao` $\phi_{g \mu}^*$ goes up to the 20 cases of third-order gradients. See also the function `get_hess_ncomp_ao_dm0`.
- For GGA/MGGA, the FLOPs of this step are $8 n_\mathrm{basis}^2 n_\mathrm{grids}$.

```rust
let ao_dm0 = index!(ao, ..ncomp_ao_dm0) % &dm0;
```

:::{admonition} Implementation decision: contract the orbitals with the density matrix in the AO basis
:class: note

Decision: use the variable `ao_dm0` $\bar{\phi}_{g \mu}^*$ below as an important intermediate.

Rejected alternative: use the occupied molecular orbital basis $\bar{\phi}_{g i}^* = \sum_{\nu} C_{\nu i} \phi_{g \nu}^*$ for the pre-contraction.

No detailed analysis or benchmarking was done here; but for the following considerations, this possibility was rejected:

- Primary reason: the `vxc` diagonal contribution `dao_vxc_diag` likewise needs `ao_dm0` $\bar{\phi}_{g \mu}^*$. This term cannot exploit the fact that there are fewer occupied orbitals than basis functions to gain speed: it must be expanded back to the atom-centered orbital basis set for the subsequent computation (both the atom-sliced contraction over $\delta_{\mu \in A}$ and the reuse of $\mathscr{T}_{\mu}^{(ts)}$ by the grid-shift contribution require keeping the $\mu$ index).
- Primary reason: the memory footprint is acceptable. If we allow the storage of `ao` $\phi_{g \mu}^*$ of size $20 n_\mathrm{basis} n_\mathrm{grids}$, then `ao_dm0` $\bar{\phi}_{g \mu}^*$ of size $4 n_\mathrm{basis} n_\mathrm{grids}$ should also be acceptable.
- Secondary reason: the program would become more complicated to write. For energy or response problems, contracting the basis functions on grids, `ao`, with the occupied orbitals might be more advantageous (in the nimatmul module, there are some `_bra_trans` functions dedicated to such problems). In gradient problems, however, partial basis contractions targeting a specific atom, $\delta_{\mu \in A}$, appear frequently, which makes the programming rather cumbersome.
- Possible improvement: using one contraction over occupied orbitals together with one expansion back, for fairly large basis sets (beyond 6-31G), the FLOPs are the smaller $16 n_\mathrm{occ} n_\mathrm{basis} n_\mathrm{grids}$. The possibility of this improvement comes from the index structure of `drho`: the indices of `drho` are $(g, \chi, t, A)$ and contain no $\mu$, so the side being contracted can equally be taken as occupied orbitals or as the atom-centered orbital basis; this differs from the case of `dao_vxc_diag`, which must keep the $\mu$ index. However, considering that this cost is not large compared with the first-order Fock skeleton gradient, and that it would require explicitly introducing occupied orbitals rather than the more general density matrix, this minor performance optimization, which introduces code complexity, was not adopted.
:::

### 3.3 Computation of `drho`

**Function `get_drho`**

| Variable | Meaning | Index order | Shape | Notes |
|--|--|--|--|--|
| `xc_type` | | | `LDA` / `GGA` / `MGGA` | |
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` | `ncomp`<br>4/10/10 |
| `ao_dm0` | $\bar{\phi}_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` | `ncomp`<br>1/4/4 |
| `aoslices` | | | `natm` |
| `drho`<br>(output) | $\partial_{A_t} \xi_g^\chi$ | $(g, \chi, t, A)$ <br>`[g, x, t, A]` | `[ngrids, nvar, 3, natm]` | |

In the table, the `ncomp` column of the `ao`/`ao_dm0` rows (given for the three cases LDA/GGA/MGGA) denotes the maximum number of derivative components actually read by this function; the actually allocated `ao` tensor takes the maximum requirement over all functions, as decided by `get_hess_ao_deriv` (for LDA, up to second derivatives, 10 components in total; for GGA/MGGA, up to third order, 20 components in total). The `ncomp` column in the function tables below follows the same convention.

**`drho`: LDA (RHO)**

$$
\begin{aligned}
\partial_{A_t} \xi_g^{\chi = \rho}
&= \sum_{\mu \nu} \partial_{A_t} \big( \phi_{g \mu} \phi_{g \nu} D_{\mu \nu} \big) \quad (\texttt{drho[:, 0, :, :]}) \\
&= - 2 \sum_{\mu} \delta_{\mu \in A} \phi_{g \mu}^t \bar{\phi}_{g \mu}
\end{aligned}
$$

The factor of 2 in the formula above comes from merging the $\mu, \nu$ symmetric terms produced by the chain rule of the partial derivatives.

**`drho`: GGA (SIGMA)**

$$
\begin{aligned}
\partial_{A_t} \xi_g^{\chi = \rho_r}
&= 2 \sum_{\mu \nu} \partial_{A_t} \big( \phi_{g \mu}^r \phi_{g \nu} D_{\mu \nu} \big) \quad (\texttt{drho[:, 1:4, :, :]}) \\
&= - 2 \sum_{\mu} \delta_{\mu \in A} \left( \phi_{g \mu}^{t r} \bar{\phi}_{g \mu} + \phi_{g \mu}^t \bar{\phi}_{g \mu}^r \right)
\end{aligned}
$$

Note that the derivation of the formula above skips steps. We need to use some rotations of the dummy indices $\mu \leftrightarrow \nu$ to simplify to the formula above (left as an exercise: why is the second term of the formula above not $\phi_{g \mu}^r \bar{\phi}_{g \mu}^t$?). The factor of 2 comes from the symmetry of $\phi_{g \mu}^r \phi_{g \nu} + \phi_{g \mu} \phi_{g \nu}^r$, which is slightly different from the LDA (RHO) case.

**`drho`: MGGA (TAU)**

$$
\begin{aligned}
\partial_{A_t} \xi_g^{\chi = \rho_\tau}
&= \frac{1}{2} \sum_{r \mu \nu} \partial_{A_t} \big( \phi_{g \mu}^r \phi_{g \nu}^r D_{\mu \nu} \big) \quad (\texttt{drho[:, 4, :, :]}) \\
&= - \sum_{r \mu} \delta_{\mu \in A} \phi_{g \mu}^{t r} \bar{\phi}_{g \mu}^r
\end{aligned}
$$

Here the final expression no longer carries a factor of 2; this is simply because the definition of $\tau$ contains a coefficient of $\frac{1}{2}$.

:::{admonition} Implementation decision: when to apply $\delta_{\mu \in A}$ in the `drho` implementation
:class: note

The above computation has two possible implementation strategies:

1. Implement directly following the formula, performing the contraction for the indices $\mu \in A$ in place;
2. First generate a tensor that carries the index $\mu$, of dimensions $(g, \mu, \chi, t)$, and on that basis apply $\delta_{\mu \in A}$ to perform the contraction, obtaining the final `drho` tensor of dimensions $(g, \chi, t, A)$.

We ended up choosing the first implementation. The reason is that the first implementation has a smaller write memory footprint, without several-fold large $(\mu, g)$ memory writes. Although the second is the strategy used by RI-JK or by the vxc contribution computations later on, here it is actually less suitable.
:::

## 4. Implementation Details of the Diagonal Part of `vxc`

The core of the `vxc` contribution is how to build and handle the following contraction problem:

$$
\frac{\partial^2 E^\text{xc}}{\partial A_t \partial B_s} \leftarrow \sum_{g \chi} w_g f_g^\chi \frac{\partial^2 \xi_g^\chi}{\partial A_t \partial B_s}
$$

### 4.1 Definitions of the diagonal and off-diagonal parts of `vxc`

Again, first recall the definition of the density on grids:

$$
\xi_g^\chi = \sum_{\mu \nu} D_{\mu \nu} \xi_{g \mu \nu}^\chi
$$

We are going to differentiate $\xi_{g \mu \nu}^\chi$ twice. According to the different combinations of $\mu, \nu$, we split the `vxc` contribution into two parts:

- **Diagonal part (`diag`)**: the two derivatives are taken on one and the same atomic orbital $\mu$; by the rules of nuclear derivatives, these two derivatives must be taken with respect to the same atom. By symmetry, the derivative contribution of $\nu$ is identical to that of $\mu$.
- **Off-diagonal part (`nondiag`)**: one derivative is taken on each of the two atomic orbitals $\mu, \nu$; these two derivatives may act on different atoms.

It should be pointed out that the diagonal part indeed contributes only to the $A = B$ part of the Hessian tensor of dimensions $(t, s, A, B)$, but **the off-diagonal part also contributes to the $A = B$ part**; it is just that the specific terms contributed differ. Using "diagonal" and "off-diagonal" to distinguish parts with different meanings is perhaps unfortunate, but for now there is no better naming.

Compared with J/K integrals, DFT:
- DFT is 2-center while J/K is 4-center, so there are far fewer atomic-orbital derivative combinations to consider;
- DFT must distinguish LDA/GGA/MGGA, and its expressions are nonlinear, so there are quite a few more terms to consider.

### 4.2 Introduction of the intermediate `dao_vxc_diag` $\mathscr{T}_{\mu}^{(ts)}$ for the diagonal part

Here we adopt the strategy of handling the electronic derivatives first and the nuclear derivatives afterwards. Note that, since we differentiate the nuclear coordinate twice, the two minus signs cancel each other out.

Unlike the way we obtained the first-order density just now, we do not contract all atomic-orbital indices in advance; instead, the atom-dependent $\delta_{\mu \in A}$ is applied at the very end. Recall that

$$
\xi_g^\chi = \sum_{\mu \nu} D_{\mu \nu} \xi_{g \mu \nu}^\chi
$$

Note that we are differentiating with respect to a single atom; considering the $(\mu, \nu)$ symmetry, the formula below carries a factor of 2, but the derivatives act only on $\mu$. Since this is hard to express in a formula, we can only describe the behavior of the restricted derivatives in words here.

$$
\frac{\partial^2 E^\text{xc}}{\partial A_t \partial A_s} \leftarrow 2 \sum_{g \chi \mu \nu} w_g f_g^\chi \frac{\partial^2 \xi_{g \mu \nu}^\chi}{\partial t \partial s} D_{\mu \nu} \delta_{\mu \in A} \quad \text{(restrict $\partial$ to $\mu$, \texttt{de\_vxc\_diag})}
$$

We will notice that $\nu$ in the formula above does not participate in the differentiation; therefore, marginalizing $\nu$ in advance (this is just the intermediate `ao_dm0` $\bar{\phi}_{g \mu}^{*}$ from the `fxc` computation earlier) turns the original `gu, gv -> uv` problem into a `gu, gu -> u` problem, saving one sizable $O(n_\mathrm{basis}^2 n_\mathrm{grids})$ matrix multiplication.

At the same time, this derivative is symmetric in $t, s$. Hence, the original $3 \times 3$ problem can be reduced to a 6-component problem over $(xx, xy, xz, yy, yz, zz)$, saving a little computational cost.

Therefore, we will introduce an intermediate `dao_vxc_diag` $\mathscr{T}_{\mu}^{(ts)}$ (dimensions $(\mu, (ts))$, of size $(n_\mathrm{basis}, 6)$):

$$
\mathscr{T}_{\mu}^{(ts)} = 2 \sum_{g \chi \nu} w_g f_g^\chi \frac{\partial^2 \xi_{g \mu \nu}^\chi}{\partial t \partial s} D_{\mu \nu} \quad \text{(restrict $\partial$ to $\mu$, \texttt{dao\_vxc\_diag})}
$$

Then the atom-dependent summation is introduced:

$$
\frac{\partial^2 E^\text{xc}}{\partial A_t \partial A_s} \leftarrow \sum_\mu \mathscr{T}_{\mu}^{(ts)} \delta_{\mu \in A} \quad (\texttt{de\_vxc\_diag})
$$

Obviously, this last step carries hardly any computational cost. At this point, the problem is reduced to how to evaluate the intermediate `dao_vxc_diag` $\mathscr{T}_{\mu}^{(ts)}$. In terms of computational-complexity analysis, it is actually quite similar to `fxc`: its dominant cost lies in the computation of `ao_dm0` $\bar{\phi}_{g \mu}$, and the rest consists of many intricate $O(N^2)$ contractions.

Before continuing, we point out that $w_g f_g^\chi$ always appears as a pair. Therefore, we can first multiply $w_g$ with $f_g^\chi$ and store the product in the variable `wv` or `wvxc`.

**Function `make_dao_vxc_diag`**

| Variable | Meaning | Index order | Shape | Notes |
|--|--|--|--|--|
| `xc_type` | | | `LDA` / `GGA` / `MGGA` | |
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` | `ncomp`<br>10/20/20 |
| `ao_dm0` | $\bar{\phi}_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` | `ncomp`<br>1/4/4 |
| `wv` | $w_g f_g^\chi$ | $(g, \chi)$<br>`[g, x]` | `[ngrids, nvar]` | |
| `dao_vxc_diag`<br>(output) | $\mathscr{T}_{\mu}^{(ts)}$ | $(\mu, (ts))$<br>`[u, ts]` | `[nao, 6]` | |

**Function `get_de_vxc_diag`**

| Variable | Meaning | Index order | Shape | Notes |
|--|--|--|--|--|
| `dao_vxc_diag` | $\mathscr{T}_{\mu}^{(ts)}$ | $(\mu, (ts))$<br>`[u, ts]` | `[nao, 6]` | |
| `aoslices` | | | `natm` |
| `de_vxc_diag`<br>(output) | | $(t, s, A, B)$<br>`[t, s, A, B]` | `[3, 3, natm, natm]` | $A \neq B$ zero values |

:::{admonition} Implementation decision: implement `dao_vxc_diag` and `de_vxc_diag` as separate functions
:class: note

If it were only for the second-order skeleton derivative of the energy, `dao_vxc_diag` and `de_vxc_diag` could in fact be merged into a single function. But considering the following two factors, the two functions are implemented separately:

- Implementing them separately introduces no extra computational cost.
- In the computations involving the grid shift, the function `contract_pvxc` needs the intermediate $\mathscr{T}_{\mu}^{(ts)}$ of `dao_vxc_diag`. The same applies to the intermediate $\mathscr{T}_{\mu \nu}^{ts}$ of `dao_vxc_off`. For code reuse, implementing `dao_vxc_diag` and `de_vxc_diag` separately is quite reasonable.
:::

### 4.3 Concrete formulas and implementation of `dao_vxc_diag`

Below we give the concrete formulas of `dao_vxc_diag` for LDA (RHO), GGA (SIGMA), and MGGA (TAU), respectively. Note that for GGA and MGGA tasks, the contributions below all need to be added together: GGA needs RHO + SIGMA, and MGGA needs RHO + SIGMA + TAU.

In the computations below, apart from `ao_dm0`, which as an input parameter itself has computational complexity $O(n_\mathrm{basis}^2 n_\mathrm{grid})$, i.e. $O(N^3)$ (but, being an input parameter, it is not counted toward the actual cost of the current function), all the remaining costs are fairly sizable $O(n_\mathrm{basis} n_\mathrm{grid})$, i.e. $O(N^2)$, computations.

**`dao_vxc_diag`: LDA (RHO)**

The LDA part corresponds to $\chi = \rho$, with $\xi_{g \mu \nu}^{\chi = \rho} = \phi_{g \mu} \phi_{g \nu}$. Differentiating only with respect to $\mu$,

$$
\frac{\partial^2 \xi_{g \mu \nu}^{\chi = \rho}}{\partial t \partial s} = \phi_{g \mu}^{ts} \phi_{g \nu} \quad \text{(restrict $\partial$ to $\mu$)}
$$

Substituting into the definition of $\mathscr{T}_\mu^{(ts)}$ and contracting $\nu$ onto `ao_dm0` $\bar{\phi}_{g \mu} = \sum_\nu \phi_{g \nu} D_{\mu \nu}$:

$$
\mathscr{T}_{\mu}^{(ts)} \mathrel{+}= 2 \sum_g w_g f_g^{\rho} \phi_{g \mu}^{ts} \bar{\phi}_{g \mu}
$$

In the program, $\phi_{g \mu}^{ts}$ is the 6 components of the `ao` tensor over $ts \in \{xx, xy, xz, yy, yz, zz\}$; $w_g f_g^\rho$ is just `wv[0]`.

**`dao_vxc_diag`: GGA (SIGMA)**

The GGA part corresponds to $\chi = \rho^r$, with $\xi_{g \mu \nu}^{\chi = \rho^r} = \phi_{g \mu}^r \phi_{g \nu} + \phi_{g \mu} \phi_{g \nu}^r$. Differentiating only with respect to $\mu$,

$$
\frac{\partial^2 \xi_{g \mu \nu}^{\chi = \rho^r}}{\partial t \partial s} = \phi_{g \mu}^{tsr} \phi_{g \nu} + \phi_{g \mu}^{ts} \phi_{g \nu}^r \quad \text{(restrict $\partial$ to $\mu$)}
$$

Substituting into the definition of $\mathscr{T}_\mu^{(ts)}$. In the first term, $\nu$ is contracted via $\bar{\phi}_{g \mu} = \sum_\nu \phi_{g \nu} D_{\mu \nu}$; in the second term, $\nu$ is contracted via $\bar{\phi}_{g \mu}^r = \sum_\nu \phi_{g \nu}^r D_{\mu \nu}$:

$$
\mathscr{T}_{\mu}^{(ts)} \mathrel{+}= 2 \sum_{g r} w_g f_g^{\rho^r} \left( \phi_{g \mu}^{tsr} \bar{\phi}_{g \mu} + \phi_{g \mu}^{ts} \bar{\phi}_{g \mu}^r \right)
$$

In detail:

- The first term needs $\phi_{g \mu}^{tsr}$, i.e. the third-order electronic-coordinate derivatives of the atomic orbitals. Concretely, for each $(ts)$ component, the third-order derivative components are taken from `ao` (e.g. $(ts) = (xx)$ corresponds to $\{xxx, xxy, xxz\}$); they are multiplied by `wv[1:4]` $w_g f_g^{\rho^r}$ ($r \in \{x, y, z\}$) and then contracted with $\bar{\phi}_{g \mu}$.
- In the second term, $\bar{\phi}_{g \mu}^r$ is just the $r \in \{x, y, z\}$ components of `ao_dm0`; they are multiplied by `wv[1:4]` and summed over, and then contracted with $\phi_{g \mu}^{ts}$. In fact, this contraction can share the access to $\phi_{g \mu}^{ts}$ with the LDA (RHO) part, so that the two merge into one contraction, saving computational cost.

**`dao_vxc_diag`: MGGA (TAU)**

The MGGA (TAU) part corresponds to $\chi = \tau$, with $\xi_{g \mu \nu}^{\chi = \tau} = \frac{1}{2} \sum_r \phi_{g \mu}^r \phi_{g \nu}^r$. Differentiating only with respect to $\mu$,

$$
\frac{\partial^2 \xi_{g \mu \nu}^{\chi = \tau}}{\partial t \partial s} = \frac{1}{2} \sum_r \phi_{g \mu}^{tsr} \phi_{g \nu}^r \quad \text{(restrict $\partial$ to $\mu$)}
$$

Substituting into the definition of $\mathscr{T}_\mu^{(ts)}$, contracting $\nu$ via $\bar{\phi}_{g \mu}^r = \sum_\nu \phi_{g \nu}^r D_{\mu \nu}$, and cancelling against the factor $2$ in the definition of $\mathscr{T}_\mu^{(ts)}$:

$$
\mathscr{T}_{\mu}^{(ts)} \mathrel{+}= \sum_{g r} w_g f_g^\tau \phi_{g \mu}^{tsr} \bar{\phi}_{g \mu}^r
$$

The correspondence in the program is similar to the first GGA term: for each $(ts)$ component, the 3 third-order derivative components $\phi_{g \mu}^{tsr}$ are taken from `ao`, paired one by one with the $r \in \{x, y, z\}$ components of `ao_dm0`, weighted by `wv[4]` $w_g f_g^\tau$, and contracted. This is likewise of $O(n_\mathrm{basis} n_\mathrm{grid})$ complexity.

## 5. Implementation Details of the Off-Diagonal Part of `vxc`

### 5.1 Introduction of the intermediate `dao_vxc_off` $\mathscr{T}_{\mu \nu}^{ts}$ for the off-diagonal part

$$
\frac{\partial^2 E^\text{xc}}{\partial A_t \partial B_s} \leftarrow \sum_{g \chi \mu \nu} w_g f_g^\chi \frac{\partial^2 \xi_{g \mu \nu}^\chi}{\partial t \partial s} D_{\mu \nu} \delta_{\mu \in A} \delta_{\nu \in B} + \text{swap} (A_t, B_s) \quad \text{(restrict $\partial_t$ to $\mu$, $\partial_s$ to $\nu$)}
$$

The off-diagonal part requires one derivative on each of the two atomic orbitals $\mu, \nu$. Unlike the diagonal case, it is no longer appropriate to marginalize one of the atomic orbitals in advance; instead, all the atomic orbitals need to be kept. We will define the following intermediate `dao_vxc_off` $\mathscr{T}_{\mu \nu}^{t s}$ (dimensions $(\mu, \nu, t, s)$, of size $(n_\mathrm{basis}, n_\mathrm{basis}, 3, 3)$; in the row-major Python reference implementation it is $(t, s, \mu, \nu)$, i.e. `[3, 3, nao, nao]`):

$$
\mathscr{T}_{\mu \nu}^{t s} = \sum_{g \chi} w_g f_g^\chi \frac{\partial^2 \xi_{g \mu \nu}^\chi}{\partial t \partial s} + \text{swap} (t \mu, s \nu) \quad \text{(restrict $\partial_t$ to $\mu$, $\partial_s$ to $\nu$, \texttt{dao\_vxc\_off})}
$$

Then the final second-order derivative contribution can be written as:

$$
\frac{\partial^2 E^\text{xc}}{\partial A_t \partial B_s} = \sum_{\mu \nu} \mathscr{T}_{\mu \nu}^{t s} D_{\mu \nu} \delta_{\mu \in A} \delta_{\nu \in B} \quad (\texttt{de\_vxc\_off})
$$

Note that the definition of $\mathscr{T}_{\mu \nu}^{t s}$ has already completed the symmetrization $\text{swap} (t \mu, s \nu)$. Therefore, the single contraction in the formula above (together with the convention that $\partial_t$ is restricted to $\mu$ and $\partial_s$ to $\nu$) already gives the complete off-diagonal contribution; there is no need to add an extra $\text{swap} (A_t, B_s)$ term. In the program, the $(B, A)$ block is filled directly by transposing the $(A, B)$ block; this is merely a mirroring operation in storage, not an extra summation term.

**Function `make_dao_vxc_off`**

| Variable | Meaning | Index order | Shape | Notes |
|--|--|--|--|--|
| `xc_type` | | | `LDA` / `GGA` / `MGGA` | |
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` | `ncomp`<br>4/10/10 |
| `wv` | $w_g f_g^\chi$ | $(g, \chi)$<br>`[g, x]` | `[ngrids, nvar]` | |
| `dao_vxc_off`<br>(output) | $\mathscr{T}_{\mu \nu}^{t s}$ | $(\mu, \nu, t, s)$<br>`[u, v, 3, 3]` | `[nao, nao, 3, 3]` | |

**Function `get_de_vxc_off`**

| Variable | Meaning | Index order | Shape | Notes |
|--|--|--|--|--|
| `dao_vxc_off` | $\mathscr{T}_{\mu \nu}^{ts}$ | $(\mu, \nu, t, s)$<br>`[u, v, 3, 3]` | `[nao, nao, 3, 3]` | |
| `dm0` | $D_{\mu \nu}$ | $(\mu, \nu)$<br>`[u, v]` | `[nao, nao]` | |
| `aoslices` | | | `natm` |
| `de_vxc_off`<br>(output) | | $(t, s, A, B)$<br>`[t, s, A, B]` | `[3, 3, natm, natm]` | |

### 5.2 Concrete formulas and implementation of `dao_vxc_off`

**`dao_vxc_off`: LDA (RHO)**

The LDA part corresponds to $\chi = \rho$, with $\xi_{g \mu \nu}^{\chi = \rho} = \phi_{g \mu} \phi_{g \nu}$. $\partial_t$ acts only on $\mu$ and $\partial_s$ only on $\nu$,

$$
\frac{\partial^2 \xi_{g \mu \nu}^{\chi = \rho}}{\partial t \partial s} = \phi_{g \mu}^t \phi_{g \nu}^s \quad \text{(restrict $\partial_t$ to $\mu$, $\partial_s$ to $\nu$)}
$$

Substituting into the definition of $\mathscr{T}_{\mu \nu}^{ts}$:

$$
\mathscr{T}_{\mu \nu}^{ts} \mathrel{+}= \sum_g w_g f_g^{\rho} \phi_{g \mu}^t \phi_{g \nu}^s + \text{swap} (t \mu, s \nu)
$$

In the program, this is a set of 9 matrix multiplications of $(n_\mathrm{basis}, n_\mathrm{grids})$ matrices that marginalize out the grid index. Here $\phi^t$ is the 3 components of the `ao` tensor over $t \in \{x, y, z\}$; $w_g f_g^\rho$ is just `wv[0]`.

We note that $\mathscr{T}_{\mu \nu}^{ts}$ is symmetric under the joint exchange of $(t \leftrightarrow s, \mu \leftrightarrow \nu)$ (because $\xi_{g \mu \nu}^\chi = \xi_{g \nu \mu}^\chi$, while $(t, s)$ serve only as symmetry labels). Therefore, of the original $3 \times 3 = 9$ matrix multiplications, in principle only the $6$ components $(ts) \in \{xx, xy, xz, yy, yz, zz\}$ need to be computed, and the remaining $(yx, zx, zy)$ are filled in by transposing the AO indices. In the current implementation, however, the LDA part is still computed in full with $9$ components, so as to stay consistent with the loop structure of the GGA part.

The dominant cost of this term is $3^2 \times 2 n_\mathrm{basis}^2 n_\mathrm{grid}$ FLOPs, an $O(N^3)$ complexity.

**`dao_vxc_off`: GGA (SIGMA)**

The GGA part corresponds to $\chi = \rho^r$, with $\xi_{g \mu \nu}^{\chi = \rho^r} = \phi_{g \mu}^r \phi_{g \nu} + \phi_{g \mu} \phi_{g \nu}^r$. $\partial_t$ acts only on $\mu$ and $\partial_s$ only on $\nu$,

$$
\frac{\partial^2 \xi_{g \mu \nu}^{\chi = \rho^r}}{\partial t \partial s} = \phi_{g \mu}^{t r} \phi_{g \nu}^s + \phi_{g \mu}^t \phi_{g \nu}^{s r} \quad \text{(restrict $\partial_t$ to $\mu$, $\partial_s$ to $\nu$)}
$$

Substituting into the definition of $\mathscr{T}_{\mu \nu}^{ts}$:

$$
\mathscr{T}_{\mu \nu}^{ts} \mathrel{+}= \sum_{g r} w_g f_g^{\rho^r} \left( \phi_{g \mu}^{t r} \phi_{g \nu}^s + \phi_{g \mu}^t \phi_{g \nu}^{s r} \right) + \text{swap} (t \mu, s \nu)
$$

Note that $\mathscr{T}_{\mu \nu}^{ts}$ is symmetric under the joint exchange $(t \leftrightarrow s, \mu \leftrightarrow \nu)$, and the second term above is precisely the image of the first term under this joint exchange. Therefore **only the first term needs to be computed**; a final transpose-symmetrization yields the contribution of the second term.

Furthermore, note that the first term shares the same right factor $\phi_{g \nu}^s$ as LDA (RHO); only the left factor changes from $\phi_{g \mu}^t$ to $\phi_{g \mu}^{tr}$ with an additional summation over $r$. Hence, we can first merge the LDA and GGA left factors into a single **weighted ket-side intermediate**, and then carry out the matrix multiplication with $\phi_{g \nu}^s$ in one go:

$$
\widetilde{\phi}_{g \mu}^{t} = \frac{1}{2} w_g f_g^\rho \phi_{g \mu}^t + \sum_r w_g f_g^{\rho^r} \phi_{g \mu}^{t r} \quad (\texttt{aowv})
$$

Note the origin of the $\frac{1}{2}$ in the first term: the LDA term $2 w_g f_g^\rho \phi_{g \mu}^t \phi_{g \nu}^s$ would be counted twice by the symmetrization below, so it must be pre-divided by $2$ to be merged with the GGA term into a single symmetrization. The GGA term $2 w_g f_g^{\rho^r} \phi_{g \mu}^{tr} \phi_{g \nu}^s$ itself contributes only once (the second term is supplied by the symmetrization), so no $\frac{1}{2}$ factor is needed inside $\widetilde{\phi}$.

Then carry out the matrix multiplication (note the correspondence order between the indices $t, s$ and $\mu, \nu$), and symmetrize:

$$
\mathscr{T}_{\mu \nu}^{ts} = 2 \sum_g \widetilde{\phi}_{g \mu}^t \phi_{g \nu}^s + \text{swap} (t \mu, s \nu)
$$

At this point, the LDA + GGA part of $\mathscr{T}_{\mu \nu}^{ts}$ is fully complete. Note that we computed LDA (RHO) and GGA (SIGMA) together here in order to save one set of matrix multiplications with the ket-side $\phi_{g \nu}^s$; this is the same strategy as the optimization in `dao_vxc_diag` of Section 4, where the LDA term and the second GGA term were merged into a single contraction.

The dominant cost of this term lies in the 9 matrix multiplications $(\widetilde{\phi}^t)^\dagger \phi^s$, i.e. $3^2 \times 2 n_\mathrm{basis}^2 n_\mathrm{grid}$ FLOPs, an $O(N^3)$ complexity. Here (as in the actual LDA implementation) the $(ts) \leftrightarrow (st)$ symmetry is not exploited, because $\widetilde{\phi}^t$ internally mixes the LDA and GGA terms, making it hard to exploit the $(ts)$ symmetry directly; the symmetrization is instead placed after the matrix multiplications and done by transposition.

**`dao_vxc_off`: MGGA (TAU)**

The MGGA (TAU) part corresponds to $\chi = \tau$, with $\xi_{g \mu \nu}^{\chi = \tau} = \frac{1}{2} \sum_r \phi_{g \mu}^r \phi_{g \nu}^r$. $\partial_t$ acts only on $\mu$ and $\partial_s$ only on $\nu$,

$$
\frac{\partial^2 \xi_{g \mu \nu}^{\chi = \tau}}{\partial t \partial s} = \frac{1}{2} \sum_r \phi_{g \mu}^{t r} \phi_{g \nu}^{s r} \quad \text{(restrict $\partial_t$ to $\mu$, $\partial_s$ to $\nu$)}
$$

Substituting into the definition of $\mathscr{T}_{\mu \nu}^{ts}$:

$$
\mathscr{T}_{\mu \nu}^{ts} \mathrel{+}= \frac{1}{2} \sum_{g r} w_g f_g^\tau \phi_{g \mu}^{t r} \phi_{g \nu}^{s r} + \text{swap} (t \mu, s \nu)
$$

Just like LDA, the formula above has the joint-exchange symmetry $(t \leftrightarrow s, \mu \leftrightarrow \nu)$ (note that the individual $(t, s)$ blocks are not themselves symmetric; rather, $\mathscr{T}^{ts}_{\mu \nu} = \mathscr{T}^{st}_{\nu \mu}$, with the transpose happening on the AO indices at the same time), so the 6 components $(ts) \in \{xx, xy, xz, yy, yz, zz\}$ can still be used for the computation, and the remaining $(yx, zx, zy)$ are filled in by transposing the AO indices. This computation is fairly expensive, involving $6 \times 3 \times 2 n_\mathrm{basis}^2 n_\mathrm{grid}$ FLOPs, an $O(N^3)$ complexity; here, $6$ comes from the number of $(ts)$ components, and $3$ from the number of $r$ components.

:::{admonition} Implementation decision: avoid direct computation of the second-derivative grids
:class: note

The off-diagonal part of `vxc` should be the most expensive part. What we adopt here is the $(\mu, \nu, t, s)$ intermediate strategy, which under mGGA requires 27 matrix multiplications of $2 n_\mathrm{basis}^2 n_\mathrm{grid}$ FLOPs each. This cost is considerably larger than those of `fxc` and of the diagonal part of `vxc` above.

Is there another way to solve this problem? The most straightforward way is to adopt the same strategy as `fxc`: directly compute the second-derivative density grids $\partial_{A_t} \partial_{B_s} \xi^\chi_g$, i.e. a tensor of shape $(n_\mathrm{atm}, n_\mathrm{atm}, 3, 3, n_\mathrm{var}, n_\mathrm{grids})$, and then contract it with vxc. In the concrete implementation, the components of $\partial_{A_t} \partial_{B_s} \xi^\chi_g$ can of course be given out in batches, but they involve the number of grid points, and the derivative density grids would lead to long, thin matrix multiplications whose contracted dimension is very small (the number of atomic orbitals of one atom), which may not be particularly cache friendly. Meanwhile, I believe the number of FMAs is about the same as in the intermediate strategy above.
:::
