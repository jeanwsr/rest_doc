# 隐式溶剂模型

在量子化学计算中，溶剂效应通常对分子的电子结构、能量和性质有显著影响。REST 程序支持基于极化连续介质模型 (Polarizable Continuum Model, PCM) 的隐式溶剂方法，通过在溶质分子周围构建与溶剂介电性质相关的连续介质腔体，将溶剂效应以静电势的形式嵌入自洽场计算中。

在 PCM 模型中，溶质分子被置于一个由原子范德华半径确定的分子空腔内，空腔表面上分布着感应表面电荷 $q$，用以模拟溶剂对溶质的极化响应。计算的核心方程为

$$
\textbf{K} \, q = \textbf{R} \, v
$$

其中，$v$ 是溶质在腔体表面格点上产生的静电势（包含核与电子贡献），$\textbf{K}$ 和 $\textbf{R}$ 是由所选 PCM 模型和溶剂介电常数决定的响应矩阵。求解得到的表面电荷 $q$ 被转化为有效势 $\textbf{V}_{\text{solv}}$，叠加到 Fock 矩阵中参与自洽迭代。

## 输入关键字

所有溶剂化计算的关键词均在输入卡的 `[ctrl]` 区块中声明。

| 关键词 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `solvent_model` | String | `"CPCM"` | 溶剂模型。设置后自动启用溶剂化计算 |
| `solv_epsilon` | f64 | `1.0` | 溶剂的相对介电常数（$\varepsilon_r$）。1.0 对应真空 |
| `pcm_cavity_radii` | String | `"UFF"` | 空腔原子半径方案 |
| `solvent_ri` | bool | `true` | 是否使用 RI 加速溶剂化积分计算 |
| `solv_chunk` | usize | `8` | 格点并行计算的批处理大小 |


### 常见溶剂的介电常数

| 溶剂 | $\varepsilon_r$ |
|---|---|
| 水 (Water) | 78.3553 |
| 二甲基亚砜 (DMSO) | 46.826 |
| 乙醇 (Ethanol) | 24.852 |
| 二氯甲烷 (DCM) | 8.93 |
| 四氢呋喃 (THF) | 7.4257 |
| 甲苯 (Toluene) | 2.3741 |
| 正己烷 (Hexane) | 1.8819 |

## 溶剂模型

REST 目前支持四种 PCM 溶剂模型，通过 `solvent_model` 关键词选择（不区分大小写）：

| 模型 | 输入值 | $f(\varepsilon)$ | 梯度支持 |
|---|---|---|---|
| CPCM | `"CPCM"` | $\frac{\varepsilon - 1}{\varepsilon}$ | 支持 |
| COSMO | `"COSMO"` | $\frac{\varepsilon - 1}{\varepsilon + 0.5}$ | 支持 |
| IEFPCM | `"IEFPCM"` | $\frac{\varepsilon - 1}{\varepsilon + 1}$ | 不支持 |
| SS(V)PE | `"SSVPE"` | $\frac{\varepsilon - 1}{\varepsilon + 1}$ | 不支持 |

其中，$f(\varepsilon)$ 控制溶剂极化响应的强度，不同模型在响应矩阵 $\textbf{K}$ 的构建方式上有所区别。对于大多数计算场景，CPCM 是推荐的选择，也是程序的缺省设置。COSMO 适用于高介电常数溶剂，是 CPCM 的一种变体。IEFPCM 和 SS(V)PE 提供了更为严格的边界条件处理，但目前不支持解析梯度，因此无法用于构型优化或受力计算。

## 空腔构建

分子空腔由体系中各原子的范德华球面叠加而成，表面上分布 Lebedev 球面格点用于数值求解 PCM 方程。REST 支持两种原子半径方案，通过 `pcm_cavity_radii` 关键词选择：

| 方案 | 缺省缩放因子 | 来源 |
|---|---|---|
| `"UFF"` (缺省) | 1.1 | UFF 力场半径 (J. Am. Chem. Soc. 114, 10024, 1992) |
| `"Bondi"` | 1.2 | Bondi 范德华半径 |

缩放因子用于将原子半径放大至适当的空腔大小。UFF 半径覆盖更多的元素类型，是当前程序的缺省选项。

## 用法示例

### 基本溶剂化单点能计算

以甲醛分子在水溶剂中的 CPCM 计算为例：

```toml
[ctrl]
    job_type =                  "energy"
    xc =                        "b3lyp"
    basis_path =                "def2-TZVP"
    auxbas_path =               "def2-universal-jkfit"
    solvent_model =             "CPCM"
    solv_epsilon =              78.3553

[geom]
    name = "HCHO"
    unit = "angstrom"
    position = '''
     C  0.000   0.000  -0.542
     O  0.000   0.000   0.677
     H  0.000   0.935  -1.082
     H  0.000  -0.935  -1.082
    '''
```

计算完成后，程序输出中会包含溶剂化能信息：

```text
The solvent energy    :     -0.0123456789 Ha
```

### 溶剂化梯度与构型优化

进行溶剂化环境下的受力计算或构型优化时，需注意将 `solvent_ri` 设为 `false`，因为当前溶剂化梯度的实现尚未支持 RI 加速：

```toml
[ctrl]
    job_type =                  "opt"
    xc =                        "b3lyp"
    basis_path =                "def2-TZVP"
    auxbas_path =               "def2-universal-jkfit"
    solvent_model =             "CPCM"
    solv_epsilon =              78.3553
    solvent_ri =                false

[geom]
    name = "HCHO"
    unit = "angstrom"
    position = '''
     C  0.000   0.000  -0.542
     O  0.000   0.000   0.677
     H  0.000   0.935  -1.082
     H  0.000  -0.935  -1.082
    '''
```

### 使用 Bondi 半径

对于需要使用 Bondi 半径方案的计算：

```toml
[ctrl]
    solvent_model =             "CPCM"
    solv_epsilon =              78.3553
    pcm_cavity_radii =          "Bondi"
```

## 注意事项

1. **梯度计算中的 RI 限制**：当进行受力 (`job_type = "force"`) 或构型优化 (`job_type = "opt"`) 计算时，必须将 `solvent_ri` 设为 `false`。在单点能计算中使用 `solvent_ri = true`（缺省）可以加速计算。
2. **梯度支持的模型限制**：目前仅 CPCM 和 COSMO 模型支持解析梯度。如果需要进行构型优化或受力计算，请选择这两种模型之一。
3. **介电常数**：`solv_epsilon = 1.0` 对应真空环境，不产生溶剂效应。用户应根据实际使用的溶剂设置正确的介电常数。
4. **启用方式**：设置 `solvent_model` 为有效的模型名称即可自动启用溶剂化计算，无需额外设置 `solvent_enabled = true`。
