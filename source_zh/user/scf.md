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
