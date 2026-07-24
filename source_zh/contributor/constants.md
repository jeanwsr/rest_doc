# 常量模块 (constants)

## 概述

常量模块集中定义了 REST 程序在整个代码库中使用的全局常量，包括物理与数学常量、数值阈值、libcint 接口所需的索引常量，以及元素周期表相关的数据。模块入口位于 `src/constants/mod.rs`，另有若干子模块存放规模较大的数据表。

将常量集中管理的目的在于：

- **单一数据来源**：同一物理量在全代码库中只有一个定义，避免出现数值不一致的重复硬编码；
- **可追溯性**：物理常量附带来源注释（如 CODATA 年份），便于核对与更新；
- **易于维护**：更新一处即可影响全部使用位置。

## 模块结构

| 子模块 | 文件 | 内容 |
|---|---|---|
| `mod`（根） | `src/constants/mod.rs` | 物理与数学常量、单位换算常量、数值阈值、规模常量、libcint 索引常量 |
| `element` | `src/constants/element.rs` | 元素名、原子质量与电荷、电子组态、原子半径、壳层信息 |
| `vsap` | `src/constants/vsap.rs` | VSAP（叠加原子势）初猜所用的数值势能表 |
| `c2s` | `src/constants/c2s.rs` | 笛卡尔基到球谐基的转换矩阵 |
| `cartesian_gto` | `src/constants/cartesian_gto.rs` | 笛卡尔高斯型轨道 (GTO) 的排布信息 |
| `solvent` | `src/constants/solvent.rs` | 溶剂模型所用的溶剂参数数据 |

各子模块通过 `pub use` 在根模块重导出，因此在使用时统一从 `crate::constants` 引用即可。

## 物理与数学常量

物理常量定义于 `src/constants/mod.rs`，每一项均附有来源注释。需要注意，这些常量来自**不同的 CODATA 发布年份**，并非完全内部一致（详见下表“来源”列）。

| 常量 | 数值 | 含义 | 来源 |
|---|---|---|---|
| `EV` | 27.2113845 | Hartree → eV | CODATA 2002 |
| `FQ` | 1822.8884861920776 | 原子质量单位 u 与电子质量 mₑ 之比（= 1 / electron mass in u） | CODATA 2014 |
| `E` | `std::f64::consts::E` | 自然常数 e | 数学常量 |
| `PI` | `std::f64::consts::PI` | 圆周率 π | 数学常量 |
| `LIGHT_SPEED` | 137.03599967994 | 精细结构常数的倒数 α⁻¹（原子单位下的光速） | CODATA 2006 |
| `BOHR` | 0.52917721092 | 玻尔半径（单位 Å），用于 Å↔Bohr 转换 | CODATA 2010 |
| `BOHR_SI` | `BOHR * 1e-10` | 玻尔半径（单位 m） | 派生 |
| `G_ELECTRON` | 2.00231930436182 | 电子 g 因子 | CODATA 2014 |
| `E_MASS` | 9.10938356e-31 | 电子质量（kg） | CODATA 2014 |
| `AVOGADRO` | 6.022140857e23 | 阿伏伽德罗常数 | CODATA 2014 |
| `PLANCK` | 6.626070040e-34 | 普朗克常数（J·s） | CODATA 2014 |
| `BOLTZMANN` | 1.380649e-23 | 玻尔兹曼常数（J/K） | CODATA 2018（精确值，SI 2019） |
| `CLIGHT_CMS` | 2.99792458e10 | 光速（cm/s） | SI 定义值（精确） |
| `R_GAS` | `BOLTZMANN * AVOGADRO` | 理想气体常数（J/(mol·K)） | 派生 |
| `E_CHARGE` | 1.6021766208e-19 | 元电荷（C） | CODATA 2014 |
| `DEBYE` | 3.335641e-30 | 1 德拜对应的 C·m | 由光速定义 |

```{note}
由于历史原因，上述常量取自多个不同的 CODATA 版本，其中大多数（`G_ELECTRON`、`E_MASS`、`AVOGADRO`、`PLANCK`、`E_CHARGE`）为 CODATA 2014。较新的常量（`HARTREE2J`、`HARTREE2WAVENUMBERS`、`AMU2KG`、`SPEED_OF_LIGHT`）采用 CODATA 2022。若需要统一到单一版本，请谨慎评估对回归测试基线的影响。
```

## 单位换算常量

以下常量用于不同单位制之间的转换，同样定义于 `src/constants/mod.rs`，集中在 `unit conversion` 注释分隔块内。

| 常量 | 数值 | 含义 | 来源 |
|---|---|---|---|
| `HARTREE2J` | 4.359744722206e-18 | Hartree → joule | CODATA 2022 |
| `HARTREE2KJMOL` | `HARTREE2J * AVOGADRO / 1000.0` | Hartree → kJ/mol（≈2625.5） | 派生 |
| `HARTREE2KCALMOL` | `HARTREE2KJMOL / CALORIE2J` | Hartree → kcal/mol | 派生 |
| `HARTREE2KCAL` | 627.5094841703362 | Hartree → kcal/mol（= HARTREE2KCALMOL） | 派生 |
| `HARTREE2WAVENUMBER` | 219474.63 | Hartree → cm⁻¹ | CODATA |
| `HARTREE2WAVENUMBERS` | 219474.63136314 | Hartree → cm⁻¹（更精确） | CODATA 2022 |
| `CALORIE2J` | 4.184 | calorie → joule | Wikipedia |
| `AMU2KG` | 1.66053906892e-27 | amu → kg | CODATA 2022 |
| `SPEED_OF_LIGHT` | 299792458.0 | 光速（m/s） | CODATA 2022 |
| `BOHR2ANG` | `BOHR` | Bohr → Å（别名） | 别名 |
| `AU2DEBYE` | `E_CHARGE * BOHR*1e-10 / DEBYE` | 原子单位偶极 → 德拜（≈2.541746） | 派生 |

## 数值阈值与规模常量

| 常量 | 数值 | 含义 |
|---|---|---|
| `AUXBAS_THRESHOLD` | 1.0e-10 | 辅助库仑矩阵求逆（及平方根逆）时的阈值 |
| `INVERSE_THRESHOLD` | 1.0e-10 | 矩阵求逆阈值 |
| `SQRT_THRESHOLD` | 1.0e-10 | 矩阵平方根阈值 |
| `E5` – `E9` | 1.0e5 – 1.0e9 | 常用数量级标度因子 |
| `MPI_CHUNK` | 134217728 | MPI 传输分块大小（约 1 GB） |

## libcint 索引常量

这些整型常量用于访问 libcint 的 `bas`、`atm`、`env` 数组中的固定槽位，须与 libcint 的约定保持一致。

| 类别 | 常量 |
|---|---|
| `bas` 数组索引 | `BAS_ATM`、`BAS_ANG`、`BAS_PRM`、`BAS_CTR`、`BAS_SLOTS` |
| `atm` 数组索引 | `ATM_NUC`、`ATM_ENV`、`ATM_NUC_MOD_OF`、`ATM_FRAC_CHARGE_OF`、`ATM_SLOTS` |
| ECP 相关 | `ECP_LMAX`、`NUC_ECP` |
| `env` 指针与核模型 | `PTR_EXPCUTOFF`、`PTR_COMMON_ORG`、`PTR_RINV_ORIG`、`NUC_MOD_OF`、`NUC_STAD_CHARGE`、`NUC_GAUS_CHARGE`、`NUC_FRAC_CHARGE`、`ENV_PRT_START` |

关于 libcint 数据结构的更多说明，请参考 code agent skill `rest-libcint-knowledge`（如存在）。

## 元素数据 (`element`)

元素周期表相关的数据集中在 `src/constants/element.rs`：

| 常量 / 静态量 | 类型 | 含义 |
|---|---|---|
| `SPECIES_NAME` | `[&str; 118]` | 元素符号（H 至 Og） |
| `MASS_CHARGE` | `[(f64, f64); 118]` | 各元素的（原子质量, 核电荷）；原子质量采用 IUPAC 2021 标准原子量 |
| `ELEM1ST` – `ELEM7TH` | `[&str; N]` | 按周期划分的元素符号列表 |
| `ELEMTMS` | `[&str; 40]` | 过渡金属元素列表 |
| `ATOM_CONFIGURATION` | `[[usize; 4]; 119]` | 各元素在 s/p/d/f 壳层的电子填充组态 |
| `S_SHELL`/`P_SHELL`/`D_SHELL`/`F_SHELL` | `[f64; N]` | 各角动量壳层的占据模板（用于 SAD 与 ECP） |
| `XE_SHELL`/`KR_SHELL` | `[f64; N]` | 稀有气体芯的壳层占据模板 |
| `NELE_IN_SHELLS` | `[f64; 15]` | 各壳层容纳的电子数 |
| `SPECIES_INFO` | `HashMap<&str, &(f64,f64)>` | 由元素符号查询（质量, 电荷）的映射（`lazy_static`） |
| `ATOMIC_RADII` | `HashMap<&str, &f64>` | 由元素符号查询原子半径（Å）的映射（`lazy_static`） |

## 使用约定

为保持代码库的一致性，请遵循以下约定：

1. **在文件开头 `use` 引入常量**，不要在函数体内使用完全限定路径。

    ```rust
    // 推荐
    use crate::constants::BOHR;
    // ...
    let dist_bohr = dist_ang / BOHR;

    // 不推荐
    let dist_bohr = dist_ang / crate::constants::BOHR;
    ```

2. **不要在代码中硬编码物理常量数值**（如 `0.52917721092`、`27.2114`、`1822.8885` 等）。应统一从 `constants` 模块引用。

3. 若所需常量尚不存在，请在 `src/constants/mod.rs` 中新增，并附上来源注释（例如 CODATA 年份及对应的 NIST 数值）。

## 来源与更新

物理常量的数值取自 [NIST CODATA](https://physics.nist.gov/cuu/Constants/) 数据库。历史版本可从归档文件核对，例如：

```
https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2010.txt
```

更新某一常量时，请在其行内注释中标明所采用的 CODATA 版本年份，并在提交说明中记录数值变化，以便评估其对能量、梯度、频率等计算结果及回归测试基线的影响。
