# Input Card Structure

## Molecular structure

## Basis sets

- `eri_type`: 取值String类型。自洽场运算中的四中心积分计算方法，目前REST支持：
    1. `analytic`: 四中心积分的解析计算方法，使用libcint库实现
	1. `ri-v`: 全称为resolution of identity，又名density fitting，是对四中心积分进行张量分解后的近似算法。(缺省)
	- **注意：REST中的analytic算法并未被充分优化，仅供程序开发测评使用，不建议在实际计算中使用**
- `basis_type`: 取值String类型。使用高斯基组的类型，有Spheric及Cartesian两种选择。Spheric对应球谐型基函数，Cartesian对应笛卡尔型基函数。缺省为spheric
- `basis_path`: 取值String类型，无缺省值。计算所使用的基组所在位置。若所用基组为cc-pVTZ, 则应为`{basis_set_pool}/cc-pVTZ`；若所用基组为STO-3G, 则应为`{basis_set_pool}/STO-3G`。其中`{basis_set_pool}`是具体基组文件夹所在的根目录。**注意：基组信息高度依赖于具体的计算体系，因此没有缺省值，必须在输入卡中声明**。
- `auxbas_path`: 取值String类型，无缺省值。计算所使用的辅助基组所在位置。辅助基组通常与常规基组放置在相同的文件夹下(`{basis_set_pool}`)。使用最广泛的辅助基组为`def2-SV(P)-JKFIT`，则申明方式应为`auxbas_path={basis_set_pool}/def2-SV(P)-JKFIT`。**注意：辅助基组信息高度依赖于具体的计算体系，因此没有缺省值。如果使用RI-V的近似方法，则必须在输入卡中声明。若`eri_type=anlaytic`，则无需使用辅助基组，也就不用申明auxbas_path**。

`basis_path`和`auxbas_path`中，`{basis_set_pool}`可以省略，例如只写 `cc-pVTZ`。REST 将自动从某些默认路径搜索该文件，优先级为：
  1. 环境变量 `REST_BASIS_DIR`；
  2. 内置路径 `/opt/rest_workspace/rest/basis-set-pool/` （REST docker 的默认基组路径）；
  3. `$REST_HOME/rest/basis-set-pool/`

若以上搜索都不成功则尝试使用 bse 下载。

`basis_path` 申明的路径本身已经存在时，其优先级高于以上所有的搜索路径。可以是完整路径，也可以是相对路径，例如`./my_basis`、`my_basis`。
REST程序对于基组的使用是高度自由和自定义的，可以根据具体的计算任务，从基组网站上下载、修改或者混合使用不同的基组。你所需要做的是：
  1. 在`{basis_set_pool}`基组文件夹下创建一个新的基组文件夹。比如你想使用混合基组，并取名这个混合基组名称为mix_bs_01。则需要创建一个基组文件夹为：`mkdir {basis_set_pool}/mix_bs_01`
  2. 然后将这些基组以”元素名称.json”放置在`{basis_set_pool}/mix_bs_01`的文件夹内
  3. 在输入卡内申明`basis_path = {basis_set_pool}/mix_bs_01`（自定义辅助基组则申明`auxbas_path`）
