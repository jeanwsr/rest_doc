# TDDFT MO mode: matrix-vector products

This page describes the matrix-vector products in `tddft_mode = "mo"`. Following the notation of the [RI-JK documentation](../ri.md), the core idea is to **pre-transform the three-center RI integrals into the MO basis**; afterwards every matrix-vector product is a pure DGEMM/DGEMV contraction over MO-basis tensors.

## Notation

In addition to the conventions of the [RI documentation](../ri.md), this page uses:

- occupied MO indices $i, j$, virtual MO indices $a, b$, auxiliary-basis indices $P, Q$;
- excitation amplitudes $z_{ia}^{\mathbb{A}}$, where the set index $\mathbb{A}$ runs over the trial vectors of the Davidson subspace (for the per-vector interface $n_\mathrm{set} = 1$);
- the MO-basis three-center integrals $B_{ia, P} = \sum_{\mu\nu} C_{\mu i} Y_{\mu\nu, P} C_{\nu a}$, with the metric absorbed into $Y$;
- the dimensions $n_\mathrm{occ}$, $n_\mathrm{vir}$, $n_\mathrm{aux}$, $n_\mathrm{basis}$, $n_\mathrm{grid}$, $n_\mathrm{var}$ (number of density variables: 1 for LDA, 4 for GGA).

Underlined subscripts denote indices that are batched during evaluation.

## Function `prepare_mo_data`

Function path: `ri_tddft::tddft::prepare_mo_data`

Builds the `TDDFTData`: it calls `prepare_fxc_data` (or the spin-resolved `prepare_fxc_data_unrestricted` for an unrestricted reference, stored as `fxc_u`) for the fxc kernel table, and extracts three MO-basis RI submatrices from `scf.rimatr` via `tddft_get_submatrix`, reshaping three further tensors as required by the exchange contractions. A restricted reference carries one bundle; an unrestricted reference carries one per spin sector (stored in `TDDFTData.ri_terms: Vec<RITensorTerms>`):

$$
\begin{aligned}
B_{ia, P} &\;\leftarrow\; \texttt{tddft\_get\_submatrix}(\texttt{'O'}, \texttt{'V'}) && \text{(eq.1)}
\end{aligned}
$$

| Variable (`RITensorTerms` member) | Meaning | Index order | Dimensions | Remarks |
|--|--|--|--|--|
| `coulomb` | $B_{ia, P}$ | $(P, ia)$ | $(n_\mathrm{aux}, n_\mathrm{occ} n_\mathrm{vir})$ | Coulomb |
| `oo_exch` | $B_{ij, P}$ | $(jP, i)$ | $(n_\mathrm{occ} n_\mathrm{aux}, n_\mathrm{occ})$ | A-block exchange |
| `vv_exch` | $B_{ab, P}$ | $(P a, b)$ | $(n_\mathrm{aux} n_\mathrm{vir}, n_\mathrm{vir})$ | A-block exchange |
| `ov_exch` | $B_{ia, P}$ | $(P i, a)$ | $(n_\mathrm{aux} n_\mathrm{occ}, n_\mathrm{vir})$ | B-block exchange |
| `oo_sr`/`vv_sr`/`ov_sr` | $B^{\mathrm{SR}}_{\cdots, P}$ | as above | as above | RSH only: the short-range $\mathrm{erfc}(\omega r_{12})/r_{12}$ exchange triple, built only when $|c_{SR}-c_{LR}| > 10^{-12}$ |
| `fxc` (TDDFTData) | `FXCMatvecData` | — | see below | XC kernel table (singlet); `fxc_u` for unrestricted; `None` for an HF reference (no kernel, J/K only) |

The exchange mixing coefficients are not stored at data-preparation time: they are derived per matrix-vector product from `scf.mol.xc_data` (RSH → $(c_{LR},\, c_{SR}-c_{LR})$; ordinary hybrids → $(c_x, 0)$); see "Range-separated hybrids" in [index](index.md).

## Function `a_matvec`

Function path: `ri_tddft::matvec::a_matvec` (signature `(scf, data, z, xlet)`; `xlet: char` labels the spin channel)

Computes the A-block matrix-vector product $K^\mathbb{A}_{ia} = \sum_{jb} A_{ia,jb} z^\mathbb{A}_{jb}$ in four steps:

$$
\begin{aligned}
K^\mathbb{A}_{ia} &\mathrel{{+}{=}} (\varepsilon_a - \varepsilon_i)\, z^\mathbb{A}_{ia}
&& \text{(eq.1 diagonal)} \\
\mathscr{T}_{P}^{\mathbb{A}} &= \sum_{jb} B_{jb, P}\, z_{jb}^{\mathbb{A}}
&& \text{(eq.2)} \\
K^\mathbb{A}_{ia} &\mathrel{{+}{=}} \kappa_c \sum_{P} B_{ia, P}\, \mathscr{T}_{P}^{\mathbb{A}}
&& \text{(eq.3 Coulomb)} \\
\mathscr{T}_{(Pa), j}^{\mathbb{A}} &= \sum_{b} B_{ab, P}\, z_{jb}^{\mathbb{A}}
&& \text{(eq.4)} \\
K^\mathbb{A}_{ia} &\mathrel{{+}{=}} -c_{LR} \sum_{jP} B_{ij, P}\, \mathscr{T}_{(Pa), j}^{\mathbb{A}} - (c_{SR}-c_{LR}) \sum_{jP} B^{\mathrm{SR}}_{ij, P}\, \mathscr{T}^{\mathbb{A}}_{(Pa), j}
&& \text{(eq.5 exchange, RSH)}
\end{aligned}
$$

where $\kappa_c$ is the Coulomb coupling factor determined by `xlet` (2 for `'S'` singlet, 1 for `'R'` unpolarized, 0 for `'T'` triplet; unrestricted references use unit weight). For ordinary hybrids the exchange coefficients are $(c_x, 0)$, i.e. the second term of eq.5 vanishes. eq.2–3 and eq.4–5 are each a DGEMV + DGEMM chain; the fxc contribution is added by `fxc_matvec` (see below). The B-block exchange index ordering $(ib|aj)$ is realized with the same pair of tensors plus the transposed amplitude and is not repeated here.

| Variable | Meaning | Index order | Dimensions | Remarks |
|--|--|--|--|--|
| `z` | $z_{ia}^{\mathbb{A}}$ | $(i, a)$ | $(n_\mathrm{occ} n_\mathrm{vir})$ | column-major, consistent with the amplitude printing |
| `coulomb` | $B_{ia,P}$ | $(P, ia)$ | $(n_\mathrm{aux}, n_\mathrm{occ}n_\mathrm{vir})$ | |
| `oo_exch` | $B_{ij,P}$ | $(jP, i)$ | $(n_\mathrm{occ}n_\mathrm{aux}, n_\mathrm{occ})$ | |
| `vv_exch` | $B_{ab,P}$ | $(Pa, b)$ | $(n_\mathrm{aux}n_\mathrm{vir}, n_\mathrm{vir})$ | |

| Memory type | Equation | Expression | Index order | Memory | Remarks |
|--|--|--|--|--|--|
| fixed | (eq.2) | $\mathscr{T}_{P}^{\mathbb{A}}$ | $(P)$ | $n_\mathrm{aux}$ | |
| fixed | (eq.4) | $\mathscr{T}_{(Pa),j}^{\mathbb{A}}$ | $(Pa, j)$ | $n_\mathrm{aux} n_\mathrm{vir} n_\mathrm{occ}$ | exchange temporary |
| fixed | (eq.5) | $K_{ia}^{\mathbb{A}}$ | $(i, a)$ | $n_\mathrm{occ} n_\mathrm{vir}$ | |

## Function `b_matvec`

Function path: `ri_tddft::matvec::b_matvec` (signature `(scf, data, z, xlet)`)

$$
\begin{aligned}
K^\mathbb{A}_{ia} &= \kappa_c \sum_{jb} (ia|jb)\, z_{jb}^{\mathbb{A}} - c_{LR} \sum_{jb} (ib|aj)\, z_{jb}^{\mathbb{A}} - (c_{SR}-c_{LR}) \sum_{jb} (ib|aj)_{\mathrm{SR}}\, z_{jb}^{\mathbb{A}} + f^{\mathrm{xc}}_{ia,jb} z_{jb}^{\mathbb{A}}
&& \text{(eq.1)}
\end{aligned}
$$

Differences from the A block: no diagonal term; the Coulomb part reuses the same `coulomb` contraction; the exchange index ordering is realized through `ov_exch` (plus `ov_sr` for RSH) with the transposed amplitude.

## fxc kernel table and `FXCMatvecData`

Function path: `dft::num_int::prepare_fxc_data` / `dft::num_int::fxc_matvec`

In MO mode the XC kernel is applied as "MO values on grid × kernel table". `prepare_fxc_data` evaluates the ground-state density on the grid, calls libxc for the second-order kernel, and pre-projects the occupied/virtual orbitals onto the grid:

| Member (`FXCMatvecData`) | Meaning | Dimensions | Remarks |
|--|--|--|--|
| `nvar` | number of density variables | | 1 LDA / 4 GGA |
| `ngrids` | number of grid points $n_\mathrm{grid}$ | | |
| `nocc` / `nvir` | occupied/virtual counts after frozen-core treatment | | |
| `mo_occ` | $\varphi_i(g)$ | $(n_\mathrm{occ}, n_\mathrm{grid})$ | |
| `mo_vir` | $\varphi_a(g)$ | $(n_\mathrm{vir}, n_\mathrm{grid})$ | |
| `mo_occ_grad` / `mo_vir_grad` | $\nabla\varphi(g)$ | $(n_\mathrm{occ/vir}, n_\mathrm{grid})\times 3$ | GGA only |
| `wfxc` | kernel table $w(g) f_{\alpha\beta}^{\mathrm{xc}}(g)$ | $n_\mathrm{grid} n_\mathrm{var}^2$ | f-contiguous, $g + \alpha n_\mathrm{grid} + \beta \cdot 4 n_\mathrm{grid}$ |

The kernel application (LDA shown for clarity):

$$
\begin{aligned}
\rho_z^{\mathbb{A}}(g) &= \sum_{ia} z_{ia}^{\mathbb{A}} \varphi_i(g) \varphi_a(g)
&& \text{(eq.1)} \\
v^{\mathbb{A}}(g) &= w(g) f^{\mathrm{xc}}(g)\, \rho_z^{\mathbb{A}}(g)
&& \text{(eq.2)} \\
K^\mathbb{A}_{ia} &= \sum_g \varphi_i(g)\, v^{\mathbb{A}}(g)\, \varphi_a(g)
&& \text{(eq.3)}
\end{aligned}
$$

Spin channels: the MO-mode `wfxc` table is prepared per channel by the `run_spin` closure of `tddft_main` — the singlet kernel $f_s = 2 f_u$ (`SINGLET_FXC_FACTOR`); the unpolarized ('R') response uses the bare kernel $f_u$ (factor 1). **Triplets are not supported in MO mode** — `tddft_main` raises an explicit error; use AO mode instead (see [tddft-ao](tddft-ao.md)). Unrestricted references use the spin-resolved table `fxc_u: FXCMatvecDataUnrestricted` (storing $f_{\sigma_1\sigma_2}$ per grid point, no singlet/triplet factors; see "Unrestricted references" in [index](index.md)).

## Solver interface

The MO-mode Davidson iteration (restricted and unrestricted references) uses the per-vector interface (one trial vector applied at a time):

- `solvers::davidson::davidson_solver` (alias `tda_davidson_solver`), closure type `FnMut(&Vec<f64>) -> Vec<f64>`;
- `solvers::davidson::lr_davidson_solver` for the symmetrized Casida equation of full linear response.

Both are per-column adapters around the batched cores `davidson_solver_batched` / `lr_davidson_solver_batched`. AO mode uses the batched interface to amortize the on-grid evaluation cost (see [tddft-ao](tddft-ao.md)).

Dense small-system path: for restricted references with $\dim \le 15$, full LR first builds A and B and attempts the $(\mathbf{A}-\mathbf{B})$ symmetrized reduction (`dense_lr_eigenpairs`), falling back to TDA (diagonalizing A only) when $\mathbf{A}-\mathbf{B}$ is not positive-definite; unrestricted MO references have separate relaxed dense thresholds (TDA $\dim \le 15$, LR $\dim \le 80$, the latter diagonalizing the non-Hermitian matrix directly, see [index](index.md)).
