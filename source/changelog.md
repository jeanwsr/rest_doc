# Changelog

## Generating the Changelog

The `scripts/analyze_changes.py` script in the rest-feedstock project generates commit lists across sub-repositories between two versions:

```bash
cd rest-feedstock
python3 scripts/analyze_changes.py --from <old-version> --to <new-version> --no-fetch
```

The script outputs git revisions and `git log --oneline` for each sub-repository, which can be used as the basis for curating changelog entries.

---

## v2026.1.0.7 → v2026.1.0.8

### rest (ded45015..3b2af4ff)

- !191 IYZ: fix a bug in the DFT grids memory optimization
- !190 IYZ: expose a new key `frac` to `geometric_pyo3`
- !189 Fix fchk for fractional occupation
- !188 Organize lib_rint R2 observables and regression tests; add Hirshfeld charge analysis
- !187 Bug: fix sanity check for basis projection
- !186 IYZ: free grids before the PT2 loop
- !185 Remove `statrs` dependency
- !184 IYZ: fix a bug in streaming PT2
- !182 Absorb `check_norm` into `scf_io`
- !181 IYZ: propose a new PT2 algorithm to reduce memory usage
- !180 Fix: explicitly release memory before geometry optimization
- !179 Fix occupation in r2u (RHF/ROKS → UKS initial guess)
- !177 Bug Fix: fix `scf_io` trait import issue
- !176 Unify all physical constants
- !175 Use logging in `ctrl_io` and `scf`
- !174 Add `analytic_hessian` option for geomeTRIC to use REST's analytical Hessian
- !173 IYZ: fix a bug in EDIIS
- !172 Bug Fix: ri-direct VK calculation under ROHF
- !171 Fix issue that caused `cargo test` to fail
- !170 IYZ: add molecular symmetry detection functionality
- !169 Add new SCF convergence criteria (e.g. gradient w.r.t. density matrix)
- !168 Fix a bug in dSCF (when using `guessfile`)
- !154 Hybrid functional analytical Hessian alternate module / analytical gradient module `analdrv`

### rest_regression

- !45 Add test for `basisproj` and `parsexc`
- !44 Bug Fix (regression): ri-direct VK calculation under ROHF
- !43 Clean basis file
- !42 Add Hessian/Frequencies regression tests; modify TDDFT regression tasks; add automated regression tasks for TRIC optimization

---

## v2026.1.0.6 → v2026.1.0.7

### rest (09ff1610..ded45015)

- !167 IYZ: evaluate thermo properties
- !166 Updated README.md for GW, BSE, TDDFT and Hessian

---

## v2026.1.0.5 → v2026.1.0.6

### rest (a9051c01..09ff1610)

- !165 Changed TDDFT output format; solved GW unnecessary output file issue
- !164 Fix cross-platform issue in Hessian

---

## v2026.1.0.4 → v2026.1.0.5

### rest (d5227617..a9051c01)

- !163 Fixed bug in hybrid functional Hessian
- !162 RHF and RKS analytic Hessian and frequencies are ready
- Fix libc usage in Hessian

---

## v2026.1.0.3 → v2026.1.0.4

### rest

- !160 fix some bugs in ROHF initial guess
- !159 fix a bug in RO-fchk
- !161 fix a bug in save_hamiltonian

---

## v2026.1.0.2 → v2026.1.0.3

### rest

- !158 fix for reading chkfile from fch2rest
- !157 fix a bug in saving ROHF fchk file
- !156 Feat: support for spin-free X2C method
- !155 read cwd for xc json
- !153 Support using solvent name; Support SMD and grad

---

## v2026.1.0.1 → v2026.1.0.2

### rest

- !151 Fix bugs for mGGA sparse grid, implement parallel grid iterator for KS
- !150 update and reorganize README.md
- !101 save and load basis in chkfile, basis projection
- !148 fix a bug in PT2

### rest_regression

- !40 modified H2O_dSCF_B3LYP and C4H4_mix_yamaguchi
- !39 redirect stderr, add durations
- !37 add regression tests for RHF(HI) and UKS(O2) external field gradients

---

## v2026.1.0 → v2026.1.0.1

### rest

- !147 less debug output for pcm
- !145 Implement ExtField Grad API for rhf/rks/uhf/uks cases
- !146 Gradient fix: external field + ghost atom full-scenario analytic force
- IYZ: fix some force bugs for embedded systems
- IYZ: invoke contribution of ghost potential into analytic gradient

### rest_regression

- !38 Revise files about pcm
- !36 add skip option

---

## v2025.02.7 → v2026.1.0

### rest

- !136 PCM: UFF radii added and set as default option; revised cavity radii settings in gradient; added illustration about solvent cavity radii
- !125 Solvation: RI-based single-point PCM, initial PCM gradient, solvent code refactoring and cleanup
- !112 Add lib_rint RI-r2 integration
- !132 Renormalized Doubles Scheme for FEAST calculations
- !131 IYZ: implement EDIIS and others; made EDIIS+smear_anneal depend on E(0); moved level-shift effect into DIIS
- !128 feat: add smearing (Fermi/Gaussian) to accelerate SCF convergence
- !126 windows support
- !122 make mpi optional
- !127 Sparse AO-grid-based density storage and VXC evaluation (non0tab)
- IYZ: add VXC screening with default threshold 1e-15 to accelerate DFT for large systems
- fxc optimized
- !133 Update: Support for all 118 elements and update atomic weights to IUPAC 2021 standard
- !134 IYZ: expose more geometric optimizer parameters for complex QM/MM calculations
- IYZ: turn on constraint optimization
- working hessian partial (CP-HF not integrated yet)
- !138 IYZ:: remove unexpected redundant basis set files
- !124 less output for DFA info
- !121 remove star in filename, and add alias for them
- !135 IYZ: fix a bug with VSAP when meeting with non0tab
- !129 fix bug for guess-mix starting from some cycle
- fix: ghost basis set analytic force dimension mismatch
- !109 fix: bug in ROKS case of meta-GGA
- IYZ: fix a bug in using different grid setting
- !137 fix a tiny bug in dft test

### rest_tensors

- !11 windows support

### rest_libcint

- !6 windows support

### rest_regression

- Add PCM gradient tests (CPCM and COSMO)
- !30 Add NH3_roSCAN
- !35 Update molecular structure and numerical results for NH3 RMP2 numerical dipole test
- !33 Remove bundled basis set files, use basis set names instead of absolute paths
- !34 add thread option, filter option
- fix convergence minor difference introduced in rest (commit 2ff1afc)
- !31 remove star in filename

---

## v2025.02.6 → v2025.02.7

### rest

- !120 Fix RI-JK gradient memory estimation
- !111 RSH functional first-order analytic gradient
- !110 Cholesky decomposition support for RI-JK analytic gradient
- !108 Introduce TBLIS (rstsr einsum support)
- !107 libxc refactoring: replace hand-rolled wrapper with external libxc crate
- !106 moved fxc related modules to dft/num_int.rs
- !105 complete support of RSH in parse_xc
- !104 calc nw based on occupation instead of homo in SCF
- !103 Fix partial logic issues in !102
- !102 RSH functional energy calculation
- !100 Add DFTD3/4 gradient
- !99 Fix unimported compile error
- !98 Working TDDFT and Stable Damped BSE Response; NL FEAST for dynamical BSE; GW Contour Deformation with Low-Rank Decomposition
- !86 3c-2e ERI decomposition integral (cderi/rimatr) algorithm and default replacement
- !85 for mac support
- !84 Add XYG2; Add R-XYG2
- !83 Add PCM documentation to README
- !81 PCM initial implementation (no gradient)
- !80 Add R-XYG7 and R-XYGJOS
- !78 Fix orbital indexing bug in sBGE2 for non-aufbau occupation
- !77 Update RHF/RKS/UHF/UKS gradient code; add ghost point charge force analysis; replace grid integration with analytic integration for point charge force
- !74 feature: output alpha/beta correlation orbital ranges for unrestricted DH calculations; demote some sBGE2/RPA output from level 2 to level 3
- !73 fix m06-l; panic when xc unknown
- !72 Add ghost point charge force calculation
- !70 feature: RO-sBGE2
- !69 refactor chkfile save
- !68 add def2-universal-jkfit as default; use lowercase file name for basis
- !67 search conda share path for basis
- !66 Fixed Bugs and Improved Efficiency of BSE
- !65 REST changes due to rest_libcint refactoring
- !64 fix: modify Yamaguchi spin contamination correction logic
- !63 correct DSD-PBEP95 to DSD-PBEB95; support for several DH/bDH functionals
- !62 add start_mix_cycle
- Closed shell CP-HF

### rest_regression

- !29 Fix thread count settings in some tests
- !28 RSH functional first-order analytic gradient (regression)
- !27 RSH functional energy calculation (regression)
- !26 Add DFTD3/4 gradient (regression)
- !25 Added TDDFT(TDA) Regressions for SVWN, PBE and PBE0
- !24 add test for guessfile

### rest_tensors

- !10 Modify pseudo-inverse based on truncated QR decomposition
- !9 Comment out rerun-if-changed on shared library binary

---

## v2025.02.3 → v2025.02.6

### rest

- !97 rm def2-svp-jkfit, update readme for basis input
- !96 refactor external_init_guess
- !95 Compile process optimization
- !94 address warnings, unused imports, unused deps
- !92 fix: handle element name case sensitivity for basis file lookup
- !93 fix some typos in README.md
- !91 Fix sysinfo::System::new_all hang on some machines
- !90 README compliance fix
- !88 new xc_parser
- !82 Allow using guessfile as initial chkfile
- !81 PCM initial implementation (no gradient)
- !80 Add R-XYG7 and R-XYGJOS
- !78 Fix orbital indexing bug in sBGE2 for non-aufbau occupation
- !77 Update RHF/RKS/UHF/UKS gradient code; add ghost point charge force analysis; replace grid integration with analytic integration for point charge force
- IYZ: Some modification to print energy densities of different DFT components

### rest_libcint

- !4 Fix potential rayon thread-local safety issue

### rest_tensors

- !8 fix warnings

---

## v2025.02.2 → v2025.02.3

### rest

- !74 feature: output alpha/beta correlation orbital ranges for unrestricted DH calculations; demote some sBGE2/RPA output from level 2 to level 3
- !73 fix m06-l; panic when xc unknown
- !72 Add ghost point charge force calculation
- !70 feature: RO-sBGE2

---

## v2025.02.1 → v2025.02.2

### rest

- !69 refactor chkfile save
- !68 add def2-universal-jkfit as default; use lowercase file name for basis
- IYZ:: Optimize grid integral code (powf→opowf); fix Cartesian GTO (l>=2) normalization inconsistency with program grid-valued basis definition
