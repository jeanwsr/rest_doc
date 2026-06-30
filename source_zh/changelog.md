# 更新日志

## 生成更新日志

通过 rest-feedstock 项目的 `scripts/analyze_changes.py` 脚本生成各仓库在两个版本之间的提交列表：

```bash
cd rest-feedstock
python3 scripts/analyze_changes.py --from <旧版本> --to <新版本> --no-fetch
```

脚本输出各子仓库的 git revision 和 `git log --oneline` 提交列表，可据此整理为更新日志条目。

---

## v2026.1.0.3 → v2026.1.0.4

### rest

- !160 fix some bugs in ROHF initial guess
- !159 fix a bug in RO-fchk
- !161 fix a bug in save_hamiltonian
- fix: add dimension check when reading chkfile
- fix a bug in chkfile in ROHF case

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
- !146 梯度修复：外电场 + 赝原子全场景解析力
- IYZ: fix some force bugs for embeded systems
- IYZ: to invole the contribution of ghost potential into the analytic gradient

### rest_regression

- !38 Revise files about pcm
- !36 add skip option

---

## v2025.02.7 → v2026.1.0

### rest

- !136 PCM: UFF radii added and set as default option; cavity can use UFF radii; revise cavity radii setting in grad; add illustration about solvent cavity radii
- !125 溶剂化单点RI;初步溶剂化梯度;溶剂化代码拆分与删除不必要的库
- !112 Add lib_rint RI-r2 integration
- !132 Renormalized Doubles Scheme for FEAST calculations
- !131 IYZ: implement EDIIS and others; IYZ: make ediis+smear_anneal depends on E(0); IYZ: move the level-shift effect into DIIS
- !128 feat: add smearing (Fermi/Gaussian) to accelerate SCF convergence
- !126 windows support
- !122 make mpi optional
- !127 基于AO格点稀疏化的密度存储与VXC计算 (non0tab)
- IYZ: add a screening to the generation of vxc: vxc_screen_threshold = 1.0e-15 (default)
- fxc optimized
- !133 Update: Support for all 118 elements and update atomic weights to IUPAC 2021 standard
- !134 IYZ: 为了复杂QM/MM计算，暴露更多的geometric优化参数
- IYZ: turn on the contraint optimization
- working hessian partial(CP-HF not integrated yet)
- !138 IYZ:: remove unexpected redundance basis set files
- !124 less output for DFA info
- !121 remove star in filename, and add alias for them
- !135 IYZ: fix a bug with vsap when meeting with non0tab
- !129 fix bug for guess-mix starting from some cycle
- fix: ghost basis set analytic force dimension mismatch
- !109 fix: bug in ROKS case of meta-GGA
- IYZ: fix a bug in using differnt grid setting
- !137 fix a tiny bug in dft test

### rest_tensors

- !11 windows support

### rest_libcint

- !6 windows support

### rest_regression

- 添加PCM梯度测试，包括CPCM/COSMO
- !30 Add NH3_roSCAN
- !35 对 NH3_RMP2_NUMDIPOLE 更新分子结构与数值结果
- !33 删掉基组文件，基组使用名称而非绝对路径
- !34 add thread option, filter option
- fix convergence minor difference introduced in rest (commit 2ff1afc)
- !31 remove star in filename

---

## v2025.02.6 → v2025.02.7

### rest

- !120 修正 RI-JK 梯度内存估算问题
- !111 RSH 密度泛函一阶解析梯度计算
- !110 对 RI-JK 解析梯度作 Cholesky 分解支持
- !108 引入 TBLIS (rstsr 的 einsum 支持)
- !107 libxc 使用重构：移入手写 wrapper，直接使用外部 libxc crate
- !106 moved fxc related modules to dft/num_int.rs
- !105 complete support of RSH in parse_xc
- !104 calc nw based on occupation instead of homo in scf
- !103 修正 !102 中部分代码逻辑问题
- !102 RSH 密度泛函能量计算
- !100 增加 DFTD3/4 梯度计算
- !99 修正 unimported 编译报错
- !98 Working TDDFT and Stable Damped BSE Response; NL FEAST for dynamical BSE; GW Contour Deformation with Low-Rank Decomposition
- !86 3c-2e ERI 分解积分 (cderi/rimatr) 算法实现与默认替换
- !85 for mac support
- !84 Add XYG2; Add R-XYG2
- !83 增加了README中的PCM说明
- !81 PCM初步实现，不包括梯度
- !80 Add R-XYG7 and R-XYGJOS
- !78 Fix orbital indexing bug in sbge2 for non-aufbau occupation
- !77 更新了 RHF/RKS/UHF/UKS 的梯度计算代码，增加了存在ghost点电荷情况下的原子受点电荷作用力分析；更新了点电荷受力计算代码，用解析积分替换了原来的格点积分
- !74 feature: 支持非限制性双杂化计算中输出 alpha/beta 关联轨道范围；部分原等级2的sBGE2和RPA的输出内容改为等级3输出
- !73 fix m06-l; panic when xc unknown
- !72 加入ghost点电荷的受力计算分析
- !70 feature: RO-sBGE2
- !69 refactor chkfile save
- !68 add def2-universal-jkfit as default; use lowercase file name for basis
- !67 search conda share path for basis
- !66 Fixed Bugs and Improved Efficiency of BSE
- !65 REST 受 rest_libcint 重构而作的改变
- !64 fix: 修改Yamaguchi自旋污染校正的代码逻辑
- !63 correct DSD-PBEP95 to DSD-PBEB95; support for several DH/bDH functionals
- !62 add start_mix_cycle
- Closed shell CP-HF

### rest_regression

- !29 修正部分测试线程数设置
- !28 RSH 密度泛函一阶解析梯度计算 (regression)
- !27 RSH 密度泛函能量计算 (regression)
- !26 增加 DFTD3/4 梯度计算 (regression)
- !25 Added TDDFT(TDA) Regressions for SVWN, PBE and PBE0
- !24 add test for guessfile

### rest_tensors

- !10 修改基于截断QR分解的伪逆求解
- !9 注释 build.rs 中 rerun-if-changed=librestmatr.so

---

## v2025.02.3 → v2025.02.6

### rest

- !97 rm def2-svp-jkfit, update readme for basis input
- !96 refactor external_init_guess
- !95 编译流程优化
- !94 address warnings, unused imports, unused deps
- !92 fix: handle element name case sensitivity for basis file lookup
- !93 fix some typos in README.md
- !91 修复 sysinfo::System::new_all 引起的程序卡顿
- !90 规避 README.md 合规性问题
- !88 new xc_parser
- !82 Allow using guessfile as initial chkfile
- !81 PCM初步实现，不包括梯度
- !80 Add R-XYG7 and R-XYGJOS
- !78 Fix orbital indexing bug in sbge2 for non-aufbau occupation
- !77 更新了 RHF/RKS/UHF/UKS 的梯度计算代码，增加了存在ghost点电荷情况下的原子受点电荷作用力分析；更新了点电荷受力计算代码，用解析积分替换了原来的格点积分
- IYZ: Some modification to print energy densities of different dft components

### rest_libcint

- !4 修复潜在的 rayon 并行下 thread-local 安全性问题

### rest_tensors

- !8 fix warnings

---

## v2025.02.2 → v2025.02.3

### rest

- !74 feature: 支持非限制性双杂化计算中输出 alpha/beta 关联轨道范围；部分原等级2的sBGE2和RPA的输出内容改为等级3输出
- !73 fix m06-l; panic when xc unknown
- !72 加入ghost点电荷的受力计算分析
- !70 feature: RO-sBGE2

---

## v2025.02.1 → v2025.02.2

### rest

- !69 refactor chkfile save
- !68 add def2-universal-jkfit as default; use lowercase file name for basis
- IYZ:: 优化了格点积分的生成代码 (powf→opowf)；修复 Cartesian GTO (l>=2) 归一化与程序格点化基组定义不一致的问题
