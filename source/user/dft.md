# DFT Settings

## Method

- `xc_parser`: the parser used to interpret `xc`. The default is `parse_xc`, which supports more and more complex functional inputs (including an expression syntax, see [parse_xc](../contributor/dft/parse_xc.md)); the old `legacy` parser remains available via `xc_parser=legacy`, with a smaller set of supported functionals.
- `xc`: String. The electronic-structure method to be evaluated. With `xc_parser=legacy`, the following are supported:
    1. Wavefunction methods: HF, MP2
    2. Local density functional approximations: LDA
    3. Generalized gradient approximations: BLYP, PBE, xPBE, XLYP
    4. Meta-GGA (kinetic-energy density) functionals: SCAN, M06-L, MN15-L, TPSS
    5. Hybrid functionals: B3LYP, X3LYP, PBE0, M05, M05-2X, M06, M06-2X, SCAN0, MN15
    6. Fifth-rung functionals: XYG3, XYGJOS, XYG7, xDH-PBE0, sBGE2, ZRPS, scsRPA, R-xDH7, RPA@PBE, RPA@B3LYP
    - HF, LDA, BLYP, PBE, B3LYP and PBE0 are self-consistent-field methods; if the user does not specify a basis set, def2-TZVPP is used (`basis_path = {basis_set_pool}/def2-TZVPP`)
    - MP2, XYG3, XYGJOS, XYG7, xDH-PBE0, sBGE2, ZRPS, scsRPA, R-xDH7, RPA@PBE and RPA@B3LYP are post-self-consistent-field methods; if the user does not specify a basis set, def2-QZVPP is used (`basis_path = {basis_set_pool}/def2-QZVPP`)
    - RPA@PBE and RPA@B3LYP mean that the post-SCF RPA calculation uses the orbitals of PBE and B3LYP, respectively

`xc_parser=parse_xc` (the default) is a superset of the inputs above and supports more functionals together with an expression syntax (e.g. `xc = "0.5*gga_x_pbe + 0.5*gga_c_pbe + hf*0.2"`); see [parse_xc](../contributor/dft/parse_xc.md).
- Note: the old `xc_type = "NonStandard"` input style (the `xc_namelist`/`xc_paralist`/`xc_hybrid_para` keywords) has been removed; custom linear-combination functionals should be written as parse_xc expressions.
- `empirical_dispersion`: String. Empirical dispersion correction for lower-rank density functionals (including LDA, BLYP, PBE, B3LYP, PBE0, etc.). Currently D3, D3BJ and D4 are supported. XYG3-type double hybrids such as XYG3, XYG7, XYGJOS, scsRPA, R-xDH7 and RPA do not need an empirical dispersion correction.
- `post_ai_correction`: String. AI-assisted correction. Currently only SCC15 is supported, and it can only be combined with the renormalized double-hybrid functional R-xDH7. See: Wang, Y.; Lin, Z.; Ouyang, R.; Jiang, B.; Zhang, I. Y.; Xu, X. Toward Efficient and Unified Treatment of Static and Dynamic Correlations in Generalized Kohn-Sham Density Functional Theory. JACS Au 2024, 4 (8), 3205-3216. https://doi.org/10.1021/jacsau.4c00488
- `post_xc`: Vec<String>. Using the converged orbitals and density, evaluate alternative exchange-correlation functionals. The allowed methods include all "xc" methods supported by REST.
- `post_correlation`: Vec<String>. Using the converged orbitals and density, evaluate post-SCF high-level correlation methods. Allowed methods include PT2, sBGE2, RPA, scsRPA, etc.

## Grid

- `grid_gen_level`: usize. Grid accuracy level; larger is more accurate. Default is 3.
- `pruning`: String. Grid pruning scheme used by DFT methods or the SAP initial guess. REST supports nwchem, sg1 and none (no pruning). Default is nwchem.
- `radial_grid_method`: String. Radial grid generation method. REST supports truetler, gc2nd, delley, becke, mura_knowles and lmg. Default is truetler.
