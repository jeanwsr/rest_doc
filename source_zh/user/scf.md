# SCF

## 自洽场计算相关关键词（Keyword）
- `initial_guess`: 取值String。分子体系进行自洽场运算所用的初始猜测方法。目前REST支持:
	1. `sad` : 对体系各原子进行自洽场计算得到自洽的密度矩阵后，将多个密度矩阵按顺序置于对角位置后得到初始的密度矩阵进行自洽场运算。缺省为sad
    1. `vsap`: Superposition of Atomic Potentials的初始猜测方法。采用半经验方法对体系势能项进行估计，与libcint生成的动能项进行加和后得到初始的Fock矩阵
	1. `hcore`: Hcore则对应单电子近似初猜，直接将由libcint生成的hcore矩阵作为初始猜测的fock矩阵进行计算
- `guess_mix`: 取值布尔类型，是否采用混合HOMO和LUMO的方法获得对称性破缺初猜。由此可以破坏体系的空间对称性和自旋对称性，有助于得到单重态UHF波函数。缺省为false。
- `guess_mix_theta_deg`: 取值`[f64;2]`或f64，分别设置两个自旋通道的混合角度（单位：度）。
    - 设为0.0，则表示完全不混合
    - 在0.0-90.0范围内，角度越大，表示破坏原始初猜效果越显著。一般建议取值0.0-45.0。缺省为[15.0, 15.0]
- `chkfile`: 取值String。给定初始猜测所在位置/路径。缺省为none
- `mixer`：取值String。辅助自洽场收敛的方法。目前REST支持direct，diis，linear及ddiis。Direct对应不使用辅助收敛方法，linear对应于线性辅助收敛方法，diis对应于direct inversion in the iterative subspace。Diis是有效的加速收敛方法。缺省为diis
- `mix_param`: 取值f64。Diis方法或linear方法的混合系数。缺省为0.2
- `start_diis_cycle`: 取值i32。开始使用diis加速收敛方法的循环数。缺省为2
- `num_max_diis`: 取值i32。diis空间大小。缺省为8
- `max_scf_cycle`: 取值i32。自洽场运算的最大迭代循环数。缺省为100
- `noiter`: 取值布尔类型。是否跳过自洽场运算。缺省为false
- `scf_acc_rho`: 取值f64。自洽场运算密度矩阵的收敛标准。缺省为1.0e-8
- `scf_acc_eev`: 取值f64。自洽场运算能量差平方和的收敛标准。缺省为1.0e-6
- `scf_acc_etot`: 取值f64。自洽场运算总能量的收敛标准。缺省为1.0e-8
- `level_shift`: 取值f64。对于发生近简并振荡不收敛的情况，可以采用level_shift的方式人为破坏简并，加速收敛。单位为hartree，缺省值为0.0
- `start_check_oscillation`: 取值i32。开始检查并自洽场计算不收敛发生振荡的循环数。当监控到自洽场发生振荡，SCF能量上升的情况，开启一次线性混合方案（linear）。缺省为20

## RI 设置

- `algorithm_jk`: 设置 Fock 矩阵计算中 J (Coulomb) 和 K (Exchange) 两部分的算法：
    - `ri-direct`: 强制使用 direct RI 算法。对于 RI-K 部分，取决于内存大小，可能会使用 semi-direct 算法 (储存相对较小的 $O(N^3)$ 的 $g_{\mu i, P}$)。
    - `ri-incore`: 强制使用 incore RI 算法 (储存完整的 Cholesky decomposed 3c-2e ERI $Y_{\mu \nu, P}$)。该算法对内存需求较大，但计算速度更快。
    - `ri`: 自动选择 incore 或 direct RI 算法，取决于自洽场计算前内存大小；在内存空间较大时选择更快的 incore 方法，内存空间较小时选择 ri-direct 方法。
    - `default`: 目前同 `ri`。
- `algorithm_j`: 设置 Fock 矩阵计算中 J (Coulomb) 部分的算法；该关键词是高级选项，一般用户建议使用`algorithm_jk`关键词进行整体设置。
    - 下述选项同 `algorithm_jk`：`ri-direct`, `ri-incore`, `ri`, `default`。
- `algorithm_k`: 设置 Fock 矩阵计算中 K (Exchange) 部分的算法；该关键词是高级选项，一般用户建议使用`algorithm_jk`关键词进行整体设置。
    - 下述选项同 `algorithm_jk`：`ri-direct`, `ri-incore`, `ri`, `default`。

<!-- 我们以后很可能需要增加更多 J 和 K 的算法，例如 J 的 multi-grids、multipole、McMurchie-Davidson 算法，K 的 COSX、aCOSX 算法等等，或引入李之韵等人已经实现的 ISDF。
该关键词在未来可以用于扩展这些算法。未来可能会有更多的性质 (如梯度、TDDFT 等)，该关键词也可以用于指定这些性质的 J 和 K 的计算算法。 -->

direct-RI 算法需要配合关键词 `max_memory` 使用。目前该关键词是以 MB 为单位指定可用内存大小；未来可能可以允许用户使用更方便的带单位字符串。
