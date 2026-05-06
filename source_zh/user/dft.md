# DFT Settings

## Method

- `xc_parser`：用于解析 `xc` 的 parser。默认选项为 `legacy`，支持的泛函较少；新的 parser 支持更多、更复杂的泛函输入，通过 `xc_parser=parse_xc` 调用。
- `xc`：取值String类型。调用的电子结构计算方法。目前`xc_parser=legacy` （默认情况）时，支持
    1. 波函数方法：HF、MP2
    2. 局域密度泛函近似：LDA
    3. 广义梯度泛函近似：BLYP、PBE、xPBE、XLYP
    4. 动能密度泛函近似：SCAN、M06-L、MN15-L、TPSS
    5. 杂化泛函近似：B3LYP、X3LYP、PBE0、M05、M05-2X、M06、M06-2X、SCAN0、MN15
    6. 第五阶泛函近似：XYG3、XYGJOS、XYG7、xDH-PBE0、sBGE2、ZRPS、scsRPA、R-xDH7、RPA@PBE、RPA@B3LYP
    - HF、LDA、BLYP、PBE、B3LYP、PBE0是自洽场计算方法，若用户未申明具体基组，则使用def2-TZVPP基组 (`basis_path = {basis_set_pool}/def2-TZVPP`)
    - MP2、XYG3、XYGJOS、XYG7、xDH-PBE0、sBGE2、ZRPS、scsRPA、R-xDH7、RPA@PBE、RPA@B3LYP为后自洽场计算方法。若用户未申明具体基组，则使用def2-QZVPP基组 (`basis_path = {basis_set_pool}/def2-QZVPP`)
    - RPA@PBE、RPA@B3LYP表示后自洽场RPA计算使用PBE、B3LYP方法的轨道
`xc_parser=parse_xc` 时，兼容以上输入，并支持更多泛函，详见 [parse_xc](../contributor/parse_xc.md)。
- `empirical_dispersion`:　取值为String。针对低级别密度泛函方法（包括LDA、BLYP、PBE、B3LYP、PBE0等）的经验色散校正方法。目前支持D3, D3BJ和D4。对于XYG3型双杂化泛函比如XYG3、XYG7、XYGJOS、scsRPA、R-xDH7、RPA等不需要经验色散校正
- `post_ai_correction`：取值String。AI辅助的校正方法。目前仅支持SCC15，并只能和R-xDH7重整化双杂化泛函方法相匹配。相关文章见：Wang, Y.; Lin, Z.; Ouyang, R.; Jiang, B.; Zhang, I. Y.; Xu, X. Toward Efficient and Unified Treatment of Static and Dynamic Correlations in Generalized Kohn–Sham Density Functional Theory. JACS Au 2024, 4 (8), 3205–3216. https://doi.org/10.1021/jacsau.4c00488
- `post_xc`：取值Vec\<String\>。采用自洽收敛的轨道和密度，进行不同的交换－关联泛函(xc)的计算。允许的方法包括REST支持的"xc"方法
- `post_correlation`：取值Vec\<String\>。采用自洽收敛的轨道和密度，进行后自洽场高等级相关能方法计算。允许的方法包括PT2、sBGE2、RPA、scsRPA等

## Grid

- `grid_gen_level`: 取值usize。格点精度等级，数值越大越精确。缺省为3
- `pruning`: 取值String。DFT方法或sap初猜所选用格点筛选。目前，REST支持nwchem，sg1以及none。其中none为不筛选。缺省为nwchem
- `radial_grid_method`: 取值String。径向格点的生成方法。目前REST支持truetler，gc2nd， delley, becke, mura_knowles及lmg。缺省为truetler