# parse_xc

为了解决任意的泛函输入的问题，仿照 pyscf 的 parse_xc 和 ajz34/dh 的扩展写了一个新的 xc_parser。在 REST 程序中可以如下使用
```
xc = "B3LYP"
xc_parser = "parse_xc"
```

以下是规则介绍。整个输入字符串**大小写不敏感**。

## 一级结构

xc String 可以由`,`分割为两个区域，或只有一个区域，不可以有大于两个区域。（`()`内的`,`不计入）

例如 `0.5*PBE + 0.5*B88, PBE` 或 `B3LYP`。

当有两个区域时，左边为 X(exchange) 区域，右边为 C(correlation) 区域。只有一个区域时，称为 XC 区域。

## 二级结构
对于每个区域，可以写成泛函 component 的线性组合。每一项中系数在前，component 名称在后，中间用`*`连接。可以没有系数，视为`1.0`。

例如 `PBE - 0.1*B88 + 23.45*PBEsol`。

## 三级结构
对于每个泛函 component，做如下规定。
### 名字
按以下优先级解析
1. WHITELIST_NONLIBXC，如 HF,MP2,SCSRPA,RPA,SBGE2,DFTD3,DFTD4 等。
2. 预定义的 CODES。如 `LDA -> 1, B3LYPG -> 402`。此时会直接得到 libxc name 和 id。
3. ALIAS，递归展开 alias，例如 `PBE -> PBE,PBE`。有两个区域时，X/C 区域的 alias 不能展开成其他类型的泛函，例如 "BLYP,"

而 `BLYP` 是合法的
```
Parsing xc: blyp
SCF components:
  type: Libxc, factor: 1, func: B88, id: 106, func_full_name: GGA_X_B88, 
    reference: ["A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)"]
  type: Libxc, factor: 1, func: LYP, id: 131, func_full_name: GGA_C_LYP, 
    reference: ["C. Lee, W. Yang, and R. G. Parr.,  Phys. Rev. B 37, 785 (1988)", "B. Miehlich, A. Savin, H. Stoll, and H. Preuss.,  Chem. Phys. Lett. 157, 200 (1989)"]
  Total hybrid: 0
```
4. 其他常规 libxc 泛函
* 可以是 libxc full name。例如`GGA_X_PBE`。但在 X 区域不能出现 C 类型的名字，反之亦然。XC 区域可以出现 X 或 C 类型。
* 可以带前缀 `X_`, `C_`, `XC_`。同样，不能出现和区域不符的前缀。
* 可以不带前缀，例如`PBE`。根据其所处的区域，视为带有 `X_`, `C_` 或 `XC_`前缀。
对于此类和上一类情况，将从所有可能的 full name 中查询匹配的名字。若有多个匹配的结果则报错。

以下是若干示例
```
Parsing xc: .2*hf + 0.08*lda + 0.72*b88, 0.81*lyp + 0.19*vwn3
Detected REST_DATA_DIR environment variable: /home/wsr/rest_test/xc
Loaded JSON file: "/home/wsr/rest_test/xc/xc_user.json"
SCF components:
  type: HF, factor: 0.2, func: HF, 
    reference: Not available yet
  type: Libxc, factor: 0.08, func: LDA, id: 1, func_full_name: LDA_X, 
    reference: ["P. A. M. Dirac.,  Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)", "F. Bloch.,  Z. Phys. 57, 545 (1929)"]
  type: Libxc, factor: 0.72, func: B88, id: 106, func_full_name: GGA_X_B88, 
    reference: ["A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)"]
  type: Libxc, factor: 0.81, func: LYP, id: 131, func_full_name: GGA_C_LYP, 
    reference: ["C. Lee, W. Yang, and R. G. Parr.,  Phys. Rev. B 37, 785 (1988)", "B. Miehlich, A. Savin, H. Stoll, and H. Preuss.,  Chem. Phys. Lett. 157, 200 (1989)"]
  type: Libxc, factor: 0.19, func: VWN3, id: 8, func_full_name: LDA_C_VWN_RPA, 
    reference: ["S. H. Vosko, L. Wilk, and M. Nusair.,  Can. J. Phys. 58, 1200 (1980)"]
  Total hybrid: 0.2
```
不合法的情况
```
Parsing xc: 0.5*pbe + 0.5*b88, x_pbe

thread 'main' (1822079) panicked at rest/src/dft/parse_xc/parse.rs:130:17:
Error: functional X_PBE has illegal prefix for type C
```

### NONLIBXC component

| | ComponenType| alias | param  |
| --- | --- | --- | --- |
| HF | HF | | |
| SR_HF | RSHF | | 1 positional (`omega`) |
| LR_HF | RSHF | | 1 positional (`omega`) |
| RSH | RSHF | | 3 positional (`omega`, `alpha`, `beta`) |
| MP2 | PT2 | | `os=, ss=` or 2 positional |
| MP2_OS | PT2 | | 1 positional |
| MP2_SS | PT2 | | 1 positional |
| RPA | RPA | | |
| SCSRPA | SCSRPA| | like MP2 |
| SBGE2 | SBGE2 | | like MP2 |
| DFTD3 | Disp | D3 | toml-param of dftd3-rs |
| DFTD4 | Disp | D4 | toml-param of dftd4-rs |
| VV10 | Disp | | |

注：VV10 目前只支持解析，对应的非局域相关计算暂未实现。

### 参数
每个 component 可以有参数，置于`()`内。形式只能是 positional 或 keyword 其中的一种（前者为`(0.3,0.4)`，后者为`(a=0.3,b=0.4)`）。

目前 libxc component 的参数只支持解析，尚不能发挥作用。

MP2, SCSRPA 等可以指定一些参数，具体规则见上一节的表格。

D3,D4 也支持自定义参数，具体规则见 [dftd3-rs](https://github.com/RESTGroup/dftd3-rs#example-custom-parameters-by-toml)。

### Dispersion

支持色散后缀的简洁输入。在泛函名称后附加后缀，会自动展开为对应的 Disp component：

| 后缀 | 展开结果 |
| --- | --- |
| `-D3ZERO` | `D3(version=zero, xc=...)` |
| `-D3BJ` | `D3(version=bj, xc=...)` |
| `-D4` | `D4(version=bj, xc=...)` |
| `-V` / `-VV10` | `VV10(xc=...)` |
| `-RVV10` | `rVV10(xc=...)` |

例如 `B3LYP-D3BJ` 等价于 `B3LYP + D3(version=bj, xc=b3lyp)`（然而这种等价并不总是成立，见以下的特殊情况）。
当 `xc` 含有上述后缀时，不能与 `empirical_dispersion` 共存。

对于某些特殊的色散泛函，还支持以下直接别名展开：

| 名称 | 展开结果 |
| --- | --- |
| `CF22D` | `CF22D + D3(version=zero, xc=cf22d)` |
| `WB97X-D3BJ` | `wb97x_v + D3(version=bj, xc=wb97x)` |
| `WB97X-D3` | `wb97x_d3 + D3(version=zero, xc=wb97x)` |
| `WB97X-V` | `wb97x_v + VV10(xc=wb97x_v)` |
| `WB97M-D3BJ` | `wb97m_v + D3(version=bj, xc=wb97m)` |
| `WB97M-V` | `wb97m_v + VV10(xc=wb97m_v)` |
| `SCAN-VV10` | `SCAN,SCAN_VV10 + VV10(xc=scan_vv10)` |
| `SCAN-RVV10` | `SCAN,SCAN_RVV10 + VV10(xc=scan_rvv10)` |
| `REVSCAN-VV10` | `REVSCAN,REVSCAN_VV10 + VV10(xc=revscan_vv10)` |

<!-- ## 杂化泛函
目前已实现了对于显式定义`HF`的杂化泛函的解析，例如上面`".2*HF + 0.08*LDA + 0.72*B88, 0.81*LYP + 0.19*VWN3"`的例子。
此时 `HF` 也是一个独立的 component。 -->


<!-- 结构体可以有 hybrid 参数 -->

## 多步泛函
目前采用内置 json （与ajz34/dh类似的格式）的方式支持了一些多步泛函。其解析优先级为最高。即首先会判断是否是多步泛函，然后分别执行多步或一步的解析。

例如
```
Parsing xc: xyg3
Detected REST_DATA_DIR environment variable: /home/wsr/rest_test/xc
Loaded JSON file: "/home/wsr/rest_test/xc/xc_user.json"
Detected multi-step functional xyg3, which is parsed to:
Step for SCF         : B3LYPg
Step for final energy: 0.8033*HF - 0.0140*LDA + 0.2107*B88, 0.6789*LYP + 0.3211*MP2
SCF components:
  type: Libxc, factor: 1, func: B3LYPG, id: 402, func_full_name: HYB_GGA_XC_B3LYP, 
    reference: ["P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and M. J. Frisch.,  J. Phys. Chem. 98, 11623 (1994)"]
  Total hybrid: 0.2
Final energy components:
  type: HF, factor: 0.8033, func: HF, 
    reference: Not available yet
  type: Libxc, factor: -0.014, func: LDA, id: 1, func_full_name: LDA_X, 
    reference: ["P. A. M. Dirac.,  Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)", "F. Bloch.,  Z. Phys. 57, 545 (1929)"]
  type: Libxc, factor: 0.2107, func: B88, id: 106, func_full_name: GGA_X_B88, 
    reference: ["A. D. Becke.,  Phys. Rev. A 38, 3098 (1988)"]
  type: Libxc, factor: 0.6789, func: LYP, id: 131, func_full_name: GGA_C_LYP, 
    reference: ["C. Lee, W. Yang, and R. G. Parr.,  Phys. Rev. B 37, 785 (1988)", "B. Miehlich, A. Savin, H. Stoll, and H. Preuss.,  Chem. Phys. Lett. 157, 200 (1989)"]
  type: PT2, factor: 1, func: MP2, 
    param_positional: [0.3211, 0.3211], reference: Not available yet
  Total hybrid: 0.8033
References:
  10.1073/pnas.0901093106
```

支持的双杂化泛函列表见 `src/dft/parse_xc/*json`。

若用户需要自定义多步泛函，可以指定环境变量 REST_DATA_DIR，在该目录下创建文件xc_*.json，其内容如
```
{
    "B2GP-PLYP-D4": {
        "code": "0.65*HF + 0.35*B88, 0.64*LYP + 0.36*MP2 + D4(xc=b2gpplyp)"
    },
    "R-XYG3": {
        "code_scf": "B3LYPg",
        "code": "0.8033*HF - 0.0140*LDA + 0.2107*B88, 0.6789*LYP + 0.3211*SBGE2"
    }
}
```

关于如何支持其他形式的多步泛函输入，例如 `String`, `Vec<String>`，待议。

## 结果
解析结果为结构体 `DFAdef`。
```
pub struct DFAdef {
    pub xc_scf: Option<Vec<DFAComponent>>,
    pub xc_nscf: Option<Vec<DFAComponent>>,
    pub reference: Vec<String>,
    pub spin_channel: usize,
    pub dfa_hybrid_scf: f64,                    // SCF 步的 HFX 系数
    pub dfa_rsh_scf: (Option<f64>, f64, f64),   // SCF 步的 RSH 参数 (omega, alpha, beta)
    pub dfa_hybrid_nscf: Option<f64>,           // 最终能量步的 HFX 系数（仅两步泛函）
}
```
同时有 `xc_scf` 和 `xc_nscf` 时为两步泛函；一步泛函只有前者非None。
这个结构体已经包含了 DFA4REST 中的所有信息，并且可以包含更多。

Component 结构体的定义如下
```rust
pub struct DFAComponent {
    pub factor: f64,
    pub func: String,
    pub func_full_name: String,
    pub id: usize,
    pub param_positional: Vec<f64>,
    pub param_keyword: HashMap<String, Value>,   // Value 可以是 f64 或 String
    pub component_type: ComponentType,
    pub reference: Vec<String>,
}

pub enum ComponentType {
    HF,
    RSHF,      // SR_HF, LR_HF, RSH
    PT2,
    RPA,
    SCSRPA,
    SBGE2,
    Disp,
    Libxc,
    Unknown,
}
```