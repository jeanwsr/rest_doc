# 开发者安装指南

程序开发者可以使用 Visual Studio Code 进行开发。

下述流程在 Linux/Mac 可以跑通。在 Windows 下，建议使用 WSL (Windows Subsystem of Linux) 开发。

在编译结束时，应该看到文件目录大致如下。其中 `.vscode` 与 `debug` 文件夹会在 VSCode 设置时修改与生成，在初次编译时不会出现。

```
$REST_HOME
├─ deps              # external dynamic libraries
├─ rest              # main code of rest
├─ rest_libcint      # libcint-based integral support
├─ rest_tensors      # 2d/3d-tensor support
├─ rest_workspace    # dummy, just for storing Cargo.toml
├─ target            # compiled binary for CLI
├─ Cargo.toml        # cargo workspace configuration
├─ rest_regression   # regression tests
│  └─ bench_pool
│      └─ NH3_X3LYP
│          └─ ctrl.in
├─ .vscode           # vscode configuration files
│  └─ settings.json
└─ debug             # compiled binary for VSCode
```

## 一般编译流程

1. **开发环境准备，必要的环境变量声明**

    ```bash
    export REST_HOME="<specify by yourself>"
    export HDF5_DIR="<fill in value of $CONDA_PREFIX>"
    export REST_EXT_DIR=$REST_HOME/deps
    export LD_LIBRARY_PATH=$REST_EXT_DIR:$LD_LIBRARY_PATH
    mkdir -p $REST_EXT_DIR
    ```

    对于 MacOS 系统，使用 `DYLD_LIBRARY_PATH`。

    `HDF5_DIR` 环境变量比较特殊。该变量是由 `hdf5-metno` 库引入的，并非 REST 本身引入。在一些包管理系统 (如 Ubuntu 的 apt 等) 下，开发者可能不需要指定该环境变量也能编译；但如果包管理系统不可用，就需要开发者手动设置。

    REST 编译需要使用 Fortran 编译器。默认情况下会使用 gfortran；但如果有特殊情况，请设置环境变量 `REST_FORTRAN_COMPILER` 或 `FC` 来指定 Fortran 编译器。

    除此之外，建议开发者安装比要的系统库：
    - 对于 Linux 开发者，一般需要安装 gcc, gfortran, g++, libclang. 例如，在 ubuntu 系统中 `sudo apt-get install build-essential libclang-dev`；
    - 对于 MacOS 开发者，一般需要安装 Xcode Command Line Tools (包含 clang 等)。可以通过 `xcode-select --install` 来安装。
    - 请准备好 conda 环境。一般建议 Miniconda。

2. **下载必要的代码**

    ```
    mkdir -p $REST_HOME; cd $REST_HOME
    git clone https://gitee.com/restgroup/rest_tensors.git
    git clone https://gitee.com/restgroup/rest_libcint.git
    git clone https://gitee.com/restgroup/rest.git
    git clone https://gitee.com/restgroup/rest_regression.git
    git clone https://gitee.com/restgroup/rest_workspace.git
    cp rest_workspace/Cargo.toml $REST_HOME
    ```

    如果开发者对程序性能有需求，可以在 `Cargo.toml` 设置 debug 模式的优化级别、或创建新的 target。

3. **安装库依赖**

    REST 所依赖的库一般都可以从主流的平台上下载 (conda, apt, yum, brew 等)。

    conda 目前是最方便的库平台；它可以安装几乎所有 REST 的依赖。

    一般建议开发者在一个新的虚拟环境 (名称可以是 rest-dev) 下安装依赖，以免与其他项目的依赖发生冲突。

    ```bash
    conda create -n rest-dev
    conda activate rest-dev
    ```

    随后安装依赖：

    ```bash
    conda install                         \
       python=3.11                        \
       numpy scipy h5py                   \
       "libopenblas=*=*openmp*"           \
       libxc libcint                      \
       openmpi gfortran libclang          \
       simple-dftd3 dftd4 geometric mokit \
       -c conda-forge                     \
       -c mokit                            
    ```

4. **链接库文件到 `$REST_EXT_DIR`**

    REST 在编译期的库依赖需要放置到特定文件夹 `$REST_EXT_DIR`。我们不一定必须要将库文件本身复制到当前文件夹，而只是需要在当前文件夹下建立指向库文件的链接即可。

    ```bash
    cd $REST_EXT_DIR
    for i in \
        $CONDA_PREFIX/lib/libcint*.so \
        $CONDA_PREFIX/lib/libxc*.so \
        $CONDA_PREFIX/lib/libmpi*.so \
        $CONDA_PREFIX/lib/libhdf5*.so \
        $CONDA_PREFIX/lib/libopenblas*.so \
        $CONDA_PREFIX/lib/libgomp*.so \
        $CONDA_PREFIX/lib/libpython*.so \
        $CONDA_PREFIX/lib/python3.11/site-packages/mokit/lib/librest2fch.so
    do ln -s $i; done
    ln -s libopenblas.0.so libopenblas.so
    ```

    上述命令会复制非常多冗余的链接 (不过由于是链接而不是复制，因此不占硬盘空间；当然既然是 rust debug 开发，这点外部动态库的硬盘空间可以看做微扰量)。开发者可以自己选择要具体链接哪些库。

    这里比较特殊的地方有
    - `librest2fch.so` 并非直接安装到 `$CONDA_PREFIX/lib`，需要从特定 python 版本的 site-packages 下链接。
    - `libopenblas.so` 有可能在一些 conda 发行版下是没有的，用户需要手动从 `libopenblas.0.so` 链接过来。
    - 在 MacOS 系统下，将后缀 `.so` 更换为 `.dylib`；将 `libgomp` 更换为 `libomp`。

5. **编译 REST**

    ```bash
    cargo build -p rest
    ```

    其中，`-p rest` 是为了只编译 rest 包；这也等同于

    ```bash
    cd $REST_HOME/rest
    cargo build
    ```

    产生的可执行程序是 `$REST_HOME/target/debug/rest`。

    请注意该编译产物是 debug 模式的，其性能是无法满足实际应用的。如果需要 release 模式的编译产物，请使用 `cargo build -p rest --release`，这会显著增加编译时间；release 编译产生的可执行程序是 `$REST_HOME/target/release/rest`。

    如果想执行单元测试，以 `rest/src/tests/test_parse_xc.rs` 为例，可以通过下述命令实现：

    ```bash
    cargo test --test test_parse_xc
    ```

6. **确认 REST 可以运行**

    到 rest_regression 运行任意一个测试样例，确认 rest 是否能正常工作。

    运行下述代码前，请留意 `REST_EXT_DIR` 文件夹是否已经包含在 `LD_LIBRARY_PATH` 中。

    ```bash
    cd $REST_HOME/rest_regression/bench_pool/NH3_X3LYP
    $REST_HOME/target/debug/rest  # use binary built from debug mode
    ```

    如果想确认全部测试样例，可以运行下述代码：

    ```bash
    cd $REST_HOME/rest_regression
    cargo run -- -c debug
    ```

## Visual Studio Code 设置

1. 首先确保能在命令行下编译与运行成功。

2. 建议在 workspace 级别打开 vscode (即 `$REST_HOME`)。不建议在 `rest` 主代码文件夹下打开 vscode。

3. 在 vscode 中安装 Rust 相关插件 (rust-analyzer)。

4. 在 vscode 的设置文件夹 `.vscode` 下，创建 `settings.json`，并添加如下内容：

    ```json
    {
        "rust-analyzer.cargo.extraEnv": {
            "REST_HOME": "<your $REST_HOME path>",
            "HDF5_DIR": "<fill in value of $CONDA_PREFIX>",
            "REST_EXT_DIR": "<fill in value of $REST_EXT_DIR>",
            "CARGO_TARGET_DIR": "debug",
            "LD_LIBRARY_PATH": "<fill in value of $LD_LIBRARY_PATH>",
            "PyO3_PYTHON": "<fill in value of your current python executable path>",
        },
        "rust-analyzer.runnables.extraEnv": {
            "REST_HOME": "<your $REST_HOME path>",
            "HDF5_DIR": "<fill in value of $CONDA_PREFIX>",
            "REST_EXT_DIR": "<fill in value of $REST_EXT_DIR>",
            "CARGO_TARGET_DIR": "debug",
            "LD_LIBRARY_PATH": "<fill in value of $LD_LIBRARY_PATH>",
            "PyO3_PYTHON": "<fill in value of your current python executable path>",
        },
    }
    ```

    - `rust-analyzer.cargo.extraEnv` 与 `rust-analyzer.runnables.extraEnv` 的设置一般是一样的，前者是针对 rust-analyzer 插件的编译环境设置，后者是针对 vscode 运行环境设置。
    - `HDF5_DIR` 与 `PyO3_PYTHON` 不是必须要设置的环境变量；但如果 rust-analyzer 或遇到编译困难，则需要设置。
    - 建议 `CARGO_TARGET_DIR` 设置为 `debug` (至少要与 `target` 文件夹有所区分)。这会产生两倍的硬盘占用，但可以避免 vscode 与命令行之间的编译产物冲突[^1]。

    [^1]: 由于各种可能的原因（特别是 Linux 情景下），VSCode 的编译环境会与 CLI 编译环境有些许差异。即使是非常小的差异，也可能会导致 cargo 认为需要重新编译一部分 dependency，从而 CLI 与 VSCode 共同使用同一个 target 目录很容易破坏增量编译（接近于重新编译），拖累一些开发时间。

5. 到 `rest/tests` 下的任意一个测试文件，点击 `▸ Run Test` 按钮，确认 vscode 能否编译与运行测试。

## MacOS 特殊设置

MacOS 系统默认情况下有较为严格的安全设置，因此需要作如下设置：

1. 如果遇到运行非常慢 (对比配置类似的 Linux 服务器或笔记本)
    - 确认所使用的 CLI 应用名称 (譬如 Terminal 或 iTerm)；
    - 随后前往 `Privacy & Security` → `Developer Tools`，将 CLI 应用添加到允许列表中。

2. 使用 VSCode 的情况下，如果编译通过但运行时总是遇到 `dyld: Library not loaded: ...` 的错误，建议将 System Integrity Protection (SIP) 功能关闭：
    - 重启 Mac 并进入恢复模式。
    - 打开恢复模式终端，输入 `csrutil disable`；
    - SIP 之所以会导致上述问题，是因为 SIP 会阻止子进程使用 `DYLD_LIBRARY_PATH` 环境变量 (参考文档 [Allow DYLD environment variables entitlement](https://developer.apple.com/documentation/BundleResources/Entitlements/com.apple.security.cs.allow-dyld-environment-variables))。
