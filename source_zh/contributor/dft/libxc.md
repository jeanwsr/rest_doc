# libxc 接口

REST 使用外部 crate [`libxc-rs`](https://github.com/RESTGroup/libxc-rs)（[docs.rs](https://docs.rs/libxc/latest/libxc/)）作为 libxc C 库的安全 Rust 封装。通过 Cargo feature `dynamic_loading` 在运行时动态链接，无需编译时静态链接。

## 架构

```
┌─────────────────────┐
│  DFA4REST (mod.rs)   │  核心 DFT 数据结构，通过 DFA 调用 libxc
│  libxc_itrf.rs       │  数值积分接口，为 num_int/response 提供高阶导数
└────────┬────────────┘
         │ xc_func_init / lda/gga/mgga_exc_vxc / eval_libxc_func_new
┌────────▼────────────┐
│  libxc_helper.rs     │  薄封装层，屏蔽 libxc-rs 细节
└────────┬────────────┘
         │ LibXCFunctional / compute_lda/gga/mgga
┌────────▼────────────┐
│  libxc (外部 crate)  │  安全 Rust 封装，FFI 调用 C libxc
└─────────────────────┘
```

## 关键数据结构

### `LibXCFunctional`（外部 crate `libxc`）

构造：`LibXCFunctional::from_number(id, spin)`，或通过 `xc_func_init(code, spin)`。

实现了 `Drop`，离开作用域时自动释放 C 端资源，不再需要手动调用 `xc_func_end()`。

| 方法 | 返回 | 说明 |
|---|---|---|
| `family()` | `LibXCFamily` | 泛函家族（`LDA`/`GGA`/`MGGA`/`HybGGA`/`HybMGGA`/`HybLDA`） |
| `hyb_exx_coef()` | `Option<f64>` | Exact-exchange 系数，非杂化返回 `None` |
| `cam_coef()` | `Option<(f64,f64,f64)>` | CAM 参数 `(omega, alpha, beta)` |
| `is_hyb_cam()` | `bool` | 是否为 range-separated 杂化泛函 |
| `needs_tau()` | `bool` | 是否需要动能密度（MGGA） |
| `needs_laplacian()` | `bool` | 是否需要拉普拉斯量 |
| `flags()` | `LibXCFlags` | 特性标志位（如 `VV10`） |
| `references()` | `&[LibXCReference]` | 引用文献列表 |
| `describe()` | `String` | 人类可读的描述信息 |

### `LibXCFamily` 枚举

| 旧名称（手写） | 新名称（外部 crate） |
|---|---|
| `HybridGGA` | `HybGGA` |
| `HybridMGGA` | `HybMGGA` |
| （无） | `HybLDA` |

### `DFA4REST`（`src/dft/mod.rs`）

| 方法 | 说明 |
|---|---|
| `init_libxc(code)` | 从 func ID 初始化 `LibXCFunctional` |
| `init_libxc_and_set_param(code)` | 带参数（如 omega）初始化 |
| `get_hybrid_libxc(...)` | 遍历 SCF 泛函计算总 HFX 系数 |
| `get_rsh_libxc(...)` | 遍历 SCF 泛函获取 RSH 参数 |
| `xc_exc_vxc(...)` | 按 family 分发调用 `lda/gga/mgga_exc_vxc` |

## `libxc_helper.rs` 公共 API

| 函数 | 说明 |
|---|---|
| `xc_func_init(func_id, spin)` | 从 libxc 编号初始化泛函，spin=1 非极化 / 2 极化 |
| `xc_code_fdqc(name)` | 将泛函名称（如 `"pbe"`）映射为 `[corr, exch, mixed]` 编号数组 |
| `lda_exc_vxc(func, rho)` | LDA 能量 + 势 |
| `gga_exc_vxc(func, rho, sigma)` | GGA 能量 + 势 |
| `mgga_exc_vxc(func, rho, sigma, lapl, tau)` | MGGA 能量 + 势 |
| `eval_libxc_func_new(func, spin, deriv, ...)` | 通用计算入口，输出到预分配缓冲区 |

## `libxc_itrf.rs` 数值积分接口

`eval_xc_eff(func_ids, factors, xc_type, spin, rho_array, np, deriv)` 是 `num_int.rs` 和 `response.rs` 的主要调用入口。接受多个 libxc 功能 ID 及其系数，计算 XC 能量、一阶势、二阶/三阶导数，并将 libxc 原始输出转换为 REST 内部格式。

## 接口迁移对照

| 旧接口（`XcFuncType`） | 新接口（`libxc-rs`） |
|---|---|
| `XcFuncType::xc_func_init(code, spin)` | `xc_func_init(code, spin)`（保持兼容） |
| `xc_func.xc_func_end()` | （删除，RAII 自动释放） |
| `xc_func.xc_func_family`（字段） | `xc_func.family()`（方法） |
| `xc_func.xc_hyb_exx_coeff()` | `xc_func.hyb_exx_coef()` 返回 `Option` |
| `xc_func.xc_hyb_cam_coef()` | `xc_func.cam_coef()` 返回 `Option` |
| `xc_func.is_rsh()` | `xc_func.is_hyb_cam()` |
| `xc_func.use_laplacian()` | `xc_func.needs_laplacian()` |
| `xc_func.use_kinetic_density()` | `xc_func.needs_tau()` |
| `xc_func.get_libxc_family()` | `xc_func.family()` |
| `xc_func.get_libxc_references()` | `xc_func.references()` |
| `xc_func.xc_func_info_printout()` | `println!("{}", xc_func.describe())` |
| `unsafe { xc_version(...) }` | `libxc::util::libxc_version()` |
| `unsafe { xc_functional_get_name(id) }` | `libxc::util::libxc_functional_get_name(id)` |

## 依赖配置

```toml
libxc = { version = "0.1", features = ["api-v7_0", "dynamic_loading"] }
```
