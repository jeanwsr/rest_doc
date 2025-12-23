# 自洽场计算

<!-- 介绍自洽场计算在量子化学计算程序中的重要性，引出后续章节内容 -->

自洽场 (Self-Consistent Field, SCF) 计算是量子化学计算程序最基本的功能，是后续下游方法，包括计算后自洽相关能、波函数分析、构型优化与性质计算的基石，因为它们都建立在通过程序自洽迭代获取的电子波函数之上。自洽场方法主要包括 Hartree-Fock (HF) 方法和 Kohn-Sham 密度泛函理论 (KS-DFT) 方法，它们分别通过求解单电子表示下的 Hartree-Fock-Roothaan 方程和 Kohn-Sham 方程来获得分子体系的电子结构信息，如分子轨道系数 $\textbf{C}$（本征矢量）、轨道能量 $\epsilon$（本征值）和电子密度矩阵 $P$。而在程序实现中，这被统一为如下形式的广义本征值问题：

$$
\textbf{F}[P] \textbf{C} = \textbf{S} \textbf{C} \epsilon
$$

其中，$\textbf{S}$ 是重叠积分矩阵 (Overlap)；$\textbf{F}$ 是 Fock 矩阵，即单电子哈密顿量，从形式上，一般包括如下的分量贡献

$$
\textbf{F}[P] = \textbf{H}_{core} + \textbf{J}[P] + \textbf{K}[P] + \textbf{V}_{xc}[P]
$$

$\textbf{H}_{core}$ 代表单电子的核哈密顿量，包括动能项和核-电子吸引项（外势）；$\textbf{J}$ 和 $\textbf{K}$ 分别是库仑 Coulomb 和交换 Exchange 两部分的电子-电子相互作用积分；$\textbf{V}_{xc}$ 是密度泛函理论中特有的交换-关联 (exchange-correlation, xc) 势。这些项都是在给定**基组**下的**矩阵**表示。在 Fock 矩阵中，所有的双电子项（包括交换关联矩阵）都依赖于密度矩阵 $P$，而密度矩阵又是由分子轨道系数 $\textbf{C}$ 计算得到，这形成了一个闭环的依赖关系，导致程序需要通过**迭代**的方式来求解上述方程，即**自洽**过程。

具体而言，程序会从一个**初始猜测 (initial guess)**的密度矩阵 $P$ 出发，**计算**对应的 Fock 矩阵 $\textbf{F}[P]$，然后求解广义本征值问题，得到**新一轮**的分子轨道系数 $\textbf{C}$ 和轨道能量 $\epsilon$，进而更新电子密度矩阵 $P$。这个过程会不断重复，直到密度矩阵和总能量的变化**满足预设的收敛标准**为止，从而获得最终的自洽波函数与能量。

在实际的程序实现与具体的计算场景中，自洽场计算还涉及到许多细节和优化策略，例如构造合适的初始猜测、使用辅助技术加速收敛、处理特殊的电子态计算、合理高效构建 Fock 矩阵等。这些内容将在后续章节中详细介绍。

## 自洽场初猜

作为自洽场计算的起点，初始猜测的选择对自洽过程的收敛速度和最终结果有着重要影响。合理的初始猜测可以显著加快收敛过程，而不合适的初猜可能导致收敛缓慢，甚至失败。在 REST 中，通过 `initial_guess` 关键词以字符串声明初猜方法，目前，程序支持的初猜方法包括：

1. `"sad"` 原子叠加密度矩阵 (Superposition of Atomic Densities)：对体系中的每个原子在声明的基组下进行独立自洽场计算，得到各自的原子密度矩阵，然后将这些密度矩阵按顺序置于对角位置，形成整个分子体系的初始密度矩阵。这种方法通常能够提供一个较为合理的初猜，适用于大多数分子体系，是目前程序的缺省设置。
2. `"vsap"` 原子势叠加 (Superposition of Atomic Potentials)：采用半经验方法对体系的势能项进行估计，并与由libcint生成的动能项进行加和，得到初始的 Fock 矩阵。
3. `"hcore"` 即核哈密顿：忽略所有的电子间相互作用项，直接将 $\textbf{H}_{core}$ 对角化得到初猜密度矩阵。

例如，一个采用 `vsap` 初猜的输入文件片段如下所示：

```toml
[ctrl]
    job_type =                  "energy"
    xc =                        "b3lyp"
    basis_path =                "cc-pVDZ"
    initial_guess =             "vsap"
```

### 从 chk 文件读取波函数初猜

有时，用户希望从之前任务的计算结果中读取波函数作为新任务的初猜，方便进行后续的计算，这种方式可以节省计算时间，尤其是在进行一系列具有顺序关系的计算任务时。与常见的计算化学程序一样，REST 提供了检查点文件（checkpoint file，简称 chk 文件）机制，允许用户将计算结果保存到 chk 文件中，并在后续计算中读取这些结果作为初猜。实际计算中，用户可通过 `chkfile` 关键词指定 chk 文件的路径与名字，程序会自动从该文件中提取所需的波函数信息作为初始猜测。

- `chkfile`: 取值String。给定初始猜测所在位置/路径，同时（目前）也是检查点文件的存放位置。缺省为 none。

例如，一个读取 chk 文件作为初猜的输入文件片段如下所示：

```toml
[ctrl]
    job_type =                  "energy"
    xc =                        "b3lyp"
    basis_path =                "cc-pVDZ"
    chkfile =                   "scf.chk"
```

在上述示例中，程序会从名为 `scf.chk` 的文件中读取波函数信息，并将其用作新的自洽场计算的初始猜测。在 REST 输入卡中，`chkfile` 接收的是一个路径 (path)，因此，本例中的 `scf.chk` 表示当前工作目录下的名为 `scf.chk` 的文件，它来自于前序的计算。如果 chk 文件位于其他目录，用户需要提供相应的相对路径或绝对路径。同时，用户也需要明确声明当前计算任务所选择的基组，尽管原则上它与前序计算的设置是一致的。

关于 REST 检查点文件的机制，需要额外说明的是，若在输入卡中声明了 `chkfile`，那么无论在声明的路径位置是否存在同名的 chk 文件，程序在运行结束后，都会将当前计算任务的结果保存到该路径下的 chk 文件中，意味着原有文件若存在，则会被覆盖。

借助 chk 文件，用户可以方便地连续进行一系列相关计算任务，例如从单点能计算过渡到几何优化，或从几何优化结果继续进行频率分析等，这些任务都建立在相同的波函数结果之上，因此原则上只需在第一个任务中进行完整的自洽场计算，后续任务均可通过读取 chk 文件来获得初猜波函数，从而节省计算资源与时间。在 REST 程序中，我们进一步提供了 `noiter` 关键词

- `noiter`: 取值布尔类型。是否跳过自洽场运算，缺省为false。

允许用户在程序构造完初始波函数后（通过初猜生成或读取），跳过自洽场计算，直接使用该波函数结果进行后续计算。例如，接续前述的单点能计算（波函数存在 `scf.chk` 中），直接进行几何优化的输入文件片段如下所示：

```toml
[ctrl]
    job_type =                  "opt"
    xc =                        "b3lyp"
    basis_path =                "cc-pVDZ"
    chkfile =                   "scf.chk"
    noiter  =                    true
```

利用这一机制，REST 程序通过 [mokit](https://jeanwsr.gitlab.io/mokit-doc-mdbook/) 程序接口，可实现对常见计算化学程序波函数文件的转换与读取，例如 Gaussian 的 fchk 文件、ORCA 的 gbw 文件等，这为用户在不同程序间迁移计算任务提供了便利，允许用户在 REST 中利用其他程序的自洽场计算结果，或将 REST 的计算结果导出到其他程序中进行后续分析与计算。

<!-- 举个例子 -->

## 自洽场收敛相关的设置

### 收敛判据

在常见的量子化学计算程序中，对于自洽收敛的判据建立在密度矩阵和总能量的变化上。在 REST 中，对应的收敛判据关键词为：

- `scf_acc_rho`：密度矩阵收敛标准。计算连续两次迭代中密度矩阵元素变化的 L2 范数，当该值小于设定的阈值时，认为密度矩阵已满足收敛。缺省值为1.0e-8，原子单位制 (a.u.)，使用浮点数 (f64) 声明，支持科学记数法。
- `scf_acc_etot`：总能量收敛标准。计算连续两次迭代中总能量的变化，当该值小于设定的阈值时，认为总能量已满足收敛。缺省值为1.0e-8，原子单位制 (a.u.)，使用浮点数 (f64) 声明，支持科学记数法。

需要指出，REST 程序默认的收敛标准已比较严格，对于通常的单点能计算，特别是大体系的计算，用户可以适当放宽这些收敛标准，以兼顾计算效率与结果精度，一个适当放宽的收敛标准可以是：

```toml
[ctrl]
    scf_acc_rho =                1.0e-7
    scf_acc_etot =               1.0e-6
```

此外，REST 程序还额外引入了一个关于轨道能量 $\epsilon$ 的收敛判据：

- `scf_acc_eev`：轨道能量收敛标准。计算连续两次迭代中所有轨道能量变化的平方和，当该值小于设定的阈值时，认为轨道能量已满足收敛。缺省值为1.0e-6，单位电子伏特 (ev)，使用浮点数 (f64) 声明，支持科学记数法。

引入这一判据的目的是为了更稳定的收敛本征值空间，由于 REST 程序特色的双杂化泛函与 RPA 泛函在形式上依赖于正则轨道能量，因此在自洽过程中监控轨道能量的变化，有助于确保最终收敛的波函数在轨道能量空间上也是稳定的，从而提升后续计算的准确性。需要指出的是，原理上轨道能量的收敛并非自洽场收敛的必要条件，且考虑到轨道能对于密度矩阵的依赖性，其收敛性与密度矩阵是同步的，因此这一判据可认为是额外的补充选项，用户可根据实际需求调整该判据。

### 收敛算法

在自洽场计算中，稳健的收敛算法对于确保计算的成功与效率至关重要。REST 程序目前支持三种自洽场收敛方法，可以通过 `mixer` 关键词选择合适的算法：

- `mixer`: 自洽场收敛方法，以字符串声明，选项包括
  1. "direct": 朴素的 (naive) 迭代方法，直接对角化获得新的密度矩阵。
  2. "linear": 线性混合方法 (linear mixing)，通过对新旧两套密度矩阵作线性组合作为下一轮的密度矩阵。
  3. "diis": 即 Direct Inversion in the Iterative Subspace 方法，通过构建误差向量空间并在该子空间中进行误差最小化确定迭代方向，辅助加速收敛过程。DIIS 是目前最常用且有效的自洽场收敛方法之一，也是目前程序的缺省设置。
  4. "ddiis": 即 Damped DIIS 方法，是 DIIS 方法的变种，通过引入阻尼因子来控制迭代过程中的波动，进一步提升收敛的稳定性。

针对上述的收敛方法，REST 程序还提供了一些相关的辅助关键词，以便用户根据具体计算需求进行调整，包括：

- `max_scf_cycle`: 适用于全部选项，指定自洽场计算的最大迭代循环次数，取值为整数，缺省值为 100。

- `mix_param`: 适用于 "linear" 方法，控制新旧密度矩阵的线性混合比例，缺省值为 0.2。

- `start_diis_cycle`: 适用于 "diis" 和 "ddiis" 方法，指定从第几轮迭代开始启用 DIIS 加速收敛，缺省值为 2。
- `num_max_diis`: 适用于 "diis" 和 "ddiis" 方法，指定 DIIS 子空间的维度，即使用最近多少轮次的 Fock 矩阵误差向量进行计算，缺省值为 8。

此外，对于一些特殊的场景，REST 程序还提供了额外的收敛辅助功能与机制。例如，当自洽场计算出现振荡不收敛的情况时，往往意味着波函数在两个态之间来回切换，而真实的稳定点往往位于这两个态之间。为此，REST 程序提供了振荡检测与处理机制，用户可以通过相关关键词启用该功能：

- `start_check_oscillation`: 指定从第几轮迭代开始检查自洽场计算是否发生振荡不收敛的情况，缺省值为20。

而当检测到振荡时，且当前轮次的自洽场能量上升时，程序自动开启一次线性混合方案（"linear"），线性混合的比例亦由 `mix_param` 关键词控制。

而对于分子体系中 HOMO/LUMO 存在近简并的体系，此时由于占据/未占据轨道无法充分的分离，会导致 DIIS 算法在收敛上的困难。为此，需要人为构造 HOMO-LUMO Gap 作为辅助，这通过 `level_shift` 关键词进行声明：

- `level_shift`: 人为构造的 HOMO-LUMO Gap 的数值大小，通过抬升 LUMO 轨道实现，单位 Hartree (a.u.)，缺省为 0.0。

一个使用 `level_shift` 功能辅助收敛的例子，例如对近解离的同核双原子分子 ($N_2$)：

```toml
[ctrl]
    jobtype     = "energy"
    xc          = "pbe"
    basis_path  = "def2-TZVP"
    auxbas_path = "def2-SV(P)-JKFIT"
    charge      = 0.0
    spin        = 1.0
    scf_acc_rho = 1.0e-7
    scf_acc_eev = 0.001
    max_scf_cycle = 150
    level_shift = 0.3

[geom]
    name = "N2"
    unit = "angstrom"
    position = '''
    N 0.0 0.0 0.0
    N 0.0 0.0 5.0
    '''
```

在这个例子中，不仅需要开启 `level_shift`，调整相应的收敛判据对于这一体系也有影响，需要用户在实践中探索。

## 自洽场 RI 算法设置

### RI 基本原理
<!-- RI 基本原理 -->

对于中等或大尺寸的分子体系（原子数大于 40 或重原子数大于 20），自洽场计算中，双电子积分项 $\textbf{J}$ 和 $\textbf{K}$ (对于 HF 方法和杂化泛函) 的计算消耗鉴于其 $O(N^4)$ 的计算标度将显著增长成为主要的计算瓶颈。因此，选择合适的 RI 内存算法对于这一类计算需求是重要的。在 REST 中，提供了 `algorithm_jk` 关键词，供用户指定算法选项：

- `algorithm_jk`: 设置自洽场 Fock 矩阵计算中 $\textbf{J}$ (Coulomb) 和 $\textbf{K}$ (Exchange) 两部分的算法，以字符串声明，选项包括：
  - "ri-direct": 强制使用 direct RI 算法。对于 RI-K 部分，取决于内存大小，可能会使用 semi-direct 算法 (储存相对较小的 $O(N^3)$ 的 $g_{\mu i, P}$)。
  - "ri-incore": 强制使用 incore RI 算法 (储存完整的 Cholesky decomposed 3c-2e ERI $Y_{\mu \nu, P}$)。该算法对内存需求较大，但计算速度更快。
  - "ri": 自动选择 incore 或 direct RI 算法，取决于自洽场计算前内存大小；在内存空间较大时选择更快的 incore 方法，内存空间较小时选择 ri-direct 方法。
  - "default": 目前同 "ri", REST 程序的缺省设置。

此外， REST 也提供了更细分的 $\textbf{J}$ 和 $\textbf{K}$ 部分算法选择，通过 `algorithm_j` 和 `algorithm_k` 关键词进行分别设置，这些关键词的选项与 `algorithm_jk` 相同。一般用户建议使用 `algorithm_jk` 关键词进行整体设置。需要指出的是， direct-RI 算法需要配合关键词 `max_memory` 使用。目前该关键词是以 MB 为单位指定可用内存大小；未来可能可以允许用户使用更方便的带单位字符串。

<!-- 我们以后很可能需要增加更多 J 和 K 的算法，例如 J 的 multi-grids、multipole、McMurchie-Davidson 算法，K 的 COSX、aCOSX 算法等等，或引入李之韵等人已经实现的 ISDF。
该关键词在未来可以用于扩展这些算法。未来可能会有更多的性质 (如梯度、TDDFT 等)，该关键词也可以用于指定这些性质的 J 和 K 的计算算法。 -->

## 自旋极化与自旋校正

对于自旋轨道 (spin-orbital) 的不同处理方式引向不同的自洽场方法，影响最终获得的波函数性质。常见的自旋处理方式包括**限制性闭壳层** (Restricted Closed-Shell)、**非限制性开壳层** (Unrestricted Open-Shell)、以及**限制性开壳层** (Restricted Open-Shell) 方法。总体上，限制性方法假设所有不同自旋的成对电子（$\alpha$/$\beta$）共享相同的空间轨道 (spatial orbital)，而非限制性方法允许不同自旋通道 (spin channel) 的电子占据不同的空间轨道，从而更灵活地描述自旋极化 (spin polarization) 现象。闭壳层方法适用于所有电子成对的体系，而开壳层方法则适用于存在未成对电子的体系。因此，限制性开壳层方法指的是在开壳层体系中，限制不同自旋通道的电子共享相同的空间轨道，即若体系中分别有 $N_\alpha$ 和 $N_\beta$ 个自旋电子 ($N_\alpha > N_\beta$)，则总共占据前 $N_\beta$ 个轨道的电子对共享相同的空间轨道，而剩余的 $N_\alpha - N_\beta$ 个 $\alpha$ 自旋电子则单独占据后续的空间轨道。限制性方法得到的波函数是自旋角动量算符 $\hat{S}^2$ 的本征函数，即自旋纯态 (spin pure state)，而非限制性方法得到的波函数通常不是 $\hat{S}^2$ 的本征函数，存在自旋污染 (spin contamination) 现象，但另一方面，对于相同的开壳层体系，非限制性方法通常能获得更低的基态能量，因为体系本身的真实基态波函数往往具有自旋极化特征。

在 REST 中，调用限制性与非限制性自洽场方法通过 `spin_polarization` 关键词进行设置：

- `spin_polarization`: bool 类型，设置自洽场计算中是否采用自旋非限制性方法，true 表示采用非限制性方法 (UHF/UKS)，false 表示采用限制性方法 (RHF/RKS/ROHF/ROKS)。默认缺省视用户声明的 `spin` 关键词的值（体系的自旋多重度）而定：当 `spin` 为 1 时，缺省为 false，即限制性闭壳层；当 `spin` 大于 1 时，缺省为 true，即非限制性开壳层。

因此，进行限制性开壳层计算（Spin > 1）时的输入设置为

```toml
[ctrl]
    spin =                      3.0
    spin_polarization =         false
```

为了在自旋非限制性计算中，更好地收敛到自旋极化的波函数来获得更低基态能量，使用轨道混合 (orbital mixing) 方法制备对称性破缺初猜，会有助于这一自洽过程，是量子化学计算中的常见操作。在 REST 程序中，通过 `guess_mix` 关键词启用该功能，并结合相关设置参数辅助制备对称性破缺初猜：

- `guess_mix`: 布尔类型，决定是否采用旋转混合 HOMO 和 LUMO 的方式获得对称性破缺初猜，以破坏体系的空间和自旋对称性，缺省为false。
- `guess_mix_theta_deg`: 取值为浮点数或长度为 2 的浮点数数组，设置两个自旋通道的混合角度（单位：度）。
  - 设为 0.0，则表示完全不混合
  - 取值应在 0.0-90.0 范围内，角度越大，表示破坏原始初猜效果越显著。一般建议取值 0.0-45.0。缺省为[15.0, 15.0]
- `start_mix_cycle`: 读取整数，指定从自洽场的第 N 次迭代时再启动 `guess_mix` 构造对称破缺初猜。默认值为0，即在初始初猜阶段就进行 mix。

双自由基分子的自旋单重态 (singlet) 是一类具有显著自旋极化特点的体系，例如氧气与卡宾，对应的使用 `guess_mix` 功能的 REST 输入卡如下

```toml
[ctrl]
    print_level =               2
    num_threads =               2
    xc =                        "hf"
    basis_path =                "def2-TZVP"
    auxbas_path =               "def2-SVP-JKFIT"
    charge =                    0.0
    spin =                      1.0
    spin_polarization =         true    
    guess_mix =                 true    
    guess_mix_theta_deg =       45.0 

[geom]
    name = "O2"
    unit = "angstrom"
    position = '''
    O      1.200   0.000   0.000
    O      0.000   0.000   0.000
    '''
```

这里，我们选择在初猜阶段就启用 `guess_mix`，对于当前体系的 HOMO/LUMO (类比于 $\pi_{2p_x}$ 与 $\pi_{2p_y}$，它们是简并的) 并将混合角度设置为 45 度。

```toml
[ctrl]
  job_type    = "energy"
  print_level = 2
  num_threads = 2
  xc          = "pbe0"
  basis_path  = "def2-SVP"
  auxbas_path = "def2-SVP-RIFIT"
  charge      = 0.0
  spin        = 1.0
  spin_polarization = true
  guess_mix   = true
  start_mix_cycle = 4
  guess_mix_theta_deg = [45, 30]
  
[geom]
  name = "Carbene"
  unit = "angstrom"
  position = '''
    C   0.000000000000  0.000000000000  0.000000000000
    H   1.108900000000  0.000000000000  0.000000000000
    H  -0.228470437505  0.000000000000  1.085108505720
  '''
```

在上例中，我们选择在自洽场运行的第 4 圈再作 `guess_mix`，并对不同的自旋通道施加不等价的旋转角，以更好地构造对称破缺波函数。

## Delta SCF 方法

- `force_state_occupation`: Delta SCF 方法中的电子占据数约束条件，形式为二维嵌套列表，
  
  ```text
  [
    [condition 1], [condition 2], ...
  ]
  ```

- 每一个子列表表示一个具体的约束条件，强制占据或强制空穴。具体设置格式为：

  `[reference, prev_state, prev_spin, target_spin, force_occ, force_check_min, force_check_max]`

  其中：
  
  - `reference`: （可选）取值String。Delta SCF 计算需要有一个常规的 SCF 计算结果，并以 chk 文件格式存在 `reference` 中。若省略，则默认与 `chkfile` 相同。
  - `prev_state` 和 `prev_spin`：整数，定位需要约束的轨道在 reference 中的轨道序号和自旋通道。
  - `target_spin`：（可选）整数，在 MOM 计算中，约束轨道的目标自旋通道
    - 若省略该值，则默认与 `prev_spin` 相同
    - 若给定，则会在指定自旋通道中寻找与 `prev_state`/`prev_spin` 最相似的轨道。
  - `force_occ`：浮点数，设置上述定位的轨道在计算中的占据数，1.0表示完全占据，0.0表示完全空穴。
  - `force_check_min` 和 `force_check_max`：整数范围 (range)，在 MOM 计算中设置搜索窗口，仅从这个窗口中寻找和 `prev_state`/`prev_spin` 最相似的轨道。

## 性质计算与波函数分析

在完成自洽场计算后，获取到最终的分子轨道系数 $\textbf{C}$ 和轨道能量 $\epsilon$，REST 程序也提供了部分的性质计算与波函数分析功能，如 Mulliken 电荷布居分析，它会与自洽场计算的结果一起打印输出。REST 也支持计算分子的偶极矩，除了以 `jobtype = "dipole"` 进行单独的偶极矩计算任务声明外，用户也可以通过 `outputs` 关键词在自洽场计算任务中同时请求偶极矩计算，例如：

```toml
[ctrl]
    job_type =                  "energy"
    xc =                        "b3lyp"
    basis_path =                "cc-pVDZ"
    outputs =                   ["dipole"]
```

用户也可通过  `outputs = ["force"]` 请求输出分子的受力信息，这与 `jobtype = "force"` 等价。需要指出的是，对于自洽场方法，当前 REST 程序的偶极矩与受力的计算是解析的，而对于后自洽场方法、双杂化泛函与 RPA 泛函，相关的计算是通过数值差分实现的，未来会逐步完善解析计算功能。

关键词 `outputs` 以字符串数组的形式声明所需的性质计算选项，它的主要功能是用于生成标准的波函数文件，以方便被后续的带有丰富波函数分析功能的程序调用，例如 Multiwfn、qubic-molstar 等，以补充目前 REST 程序缺乏的波函数后处理功能。关于 `outputs` 关键词的具体说明如下：

- `outputs`: 以字符串数组形式声明所需的性质计算选项，当前支持的选项包括：
  - `dipole`    偶极矩
  - `force`     分子受力
  - `fchk`　    Gaussian 程序的 fchk 文件 (formatted checkpoint file)
  - `cube_orb`  格点化的轨道文件信息 (cube file)
  - `molden`　　 molden 格式的波函数文件
  - `geometry`  分子结构文件 (xyz 格式)

另外，当输出选项为 `cube_orb` 时，用户还可以通过 `cube_orb_indices` 指定需要输出的轨道编号列表，为二维嵌套列表，

  ```text
  [
    [start, end, spin], [start, end, spin], ...
  ]
  ```

子列表格式为 `[start, end, spin]`，表示一组轨道的起始编号、截止编号和自旋通道。

例如如下苯分子的计算

```toml
[ctrl]
    jobtype     = "energy"
    print_level = 1
    num_threads = 2
    xc          = "x3lyp"
    basis_path  = "def2-SVP"
    auxbas_path = "def2-SV(P)-JKFIT"
    charge      = 0.0
    spin        = 1.0
    spin_polarization = false
    outputs     = ["cube_orb", "fchk"]
    cube_orb_indices = [[20, 21, 0]] 

[geom]
    name = "benzene"
    unit = "angstrom"
    position = '''
      C   1.217739890298750  -0.703062453466927  0.000000000000000
      H   2.172991468538160  -1.254577209307266  0.000000000000000
      C   1.217739890298750   0.703062453466927  0.000000000000000
      H   2.172991468538160   1.254577209307266  0.000000000000000
      C   0.000000000000000   1.406124906933854  0.000000000000000
      H   0.000000000000000   2.509154418614532  0.000000000000000
      C  -1.217739890298750   0.703062453466927  0.000000000000000
      H  -2.172991468538160   1.254577209307266  0.000000000000000
      C  -1.217739890298750  -0.703062453466927  0.000000000000000
    '''
```

计算完成后同时输出 fchk 和 cube 文件，其中，指定输出了苯分子的 HOMO (20) 和 LUMO (21) 轨道 (从 0 开始编号)，自旋通道为 $\alpha$ 自旋。

此外，cube 文件还接收以下两个关键词设置：

- `cube_orb_setting`: 二维浮点数数组 [margin, num_grids]。`cube_orb` 的格点参数设置。第一个值表示边界范围 (margin)，第二个值 (num_grids) 为生成格点的数目，缺省值为[3.0, 80.0]。
- `cube_orb_type`: 以字符串方式指定生成的 cube 文件类型，选项为
  1. "wavefunction": 生成轨道波函数的 cube 文件（默认）
  2. "density": 生成轨道概率密度的 cube 文件（|ψ|²）
