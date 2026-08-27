# ADR-0001：OpenBLAS 的并行 BLAS 线程调度

- **状态**：使用中
- **日期**：2026-08-13
- **适用范围**：REST 本体 (crate `rest`) 的全部模块 (SCF、post-SCF、GW、微扰、TD 等)，不限于特定模块。线程控制函数由 crate `rest_tensors` 提供实现；使用 RSTSR 张量运算时的处理见第 5 节。存量代码存在尚未遵循之处，将在重构中逐步对齐。
- **作者**：祝震予 (ajz34@outlook.com)

## 决定

- REST 目前使用 **OpenMP 编译的 OpenBLAS** 作为 BLAS 后端；本 ADR 仅对该情形有约束。其他数学库目前是实验性支持 (见第 7 节备注)。
- 在并发 (rayon) 线程内部调用 BLAS 前，调用 `omp_set_num_threads_wrapper(1)`；该函数是 thread-safe 与 thread-local 的，因此无需保存/恢复线程数。
- 使用 RSTSR 的矩阵乘法与 einsum 时，无需手动设置线程数 (见第 5 节)。
- 主线程初始化时，与 rayon 全局线程池的构建一起，调用一次 `omp_set_num_threads_global_wrapper` 设置全局的 BLAS 线程数；输入卡 `num_threads` 选项覆盖 `OPENBLAS_NUM_THREADS` 等环境变量的设置 (见 3.3 节)。
- pthreads 并行 OpenBLAS 的支持已废止 (见 4.2 节)。

决定的核心模式：

```rust
(0..niter).into_par_iter().for_each(|_| {
    omp_set_num_threads_wrapper(1); // 并发线程内：BLAS 单线程
    let c = dgemm(a, b);
});
```

## 1. 设计动机

电子结构计算大量地需要使用 BLAS 函数。REST 目前主要使用 OpenBLAS 作为 BLAS 库。

BLAS 函数有两种调用情形。一种情形较为简单，即普通的矩阵乘法或矩阵分解；普通地调用即可，BLAS 以全局设定的满线程数并行。

但另一种情形可以归结为分批矩阵乘法 $C_{i j}^p = \sum_k A_{i k}^p B_{k j}^p$ 或基于此的算子融合。当指标 $p$ 的维度足够大时，我们习惯对 $p$ 进行并行迭代，而不是使用普通的矩阵乘法：
```python
# 伪代码
parallel for p in 0..N {
    C[p] = A[p] @ B[p]
}
```
在更复杂的问题中，对 $p$ 的迭代将会涉及多步矩阵乘法、甚至是复杂的算子融合。如果不对 $p$ 的迭代作并行，或者性能会下降、或者会产生较大内存占用的中间张量。

我们希望这两种情形都能有效地利用多核 CPU 进行并行。基于 OpenMP 并行的多线程 OpenBLAS 恰能同时支撑这两种用法：普通 BLAS 计算使用满线程数，并发批量计算在每个并发线程内将 BLAS 限制为单线程。如何作此调度，即本决定的内容。

## 2. 用法

下文以 `dgemm` 泛指任一多线程 BLAS/LAPACK 调用。在 rayon 并行区域之外，无需任何设置，BLAS 按初始化时设定的全局线程数运行 (见 3.3 节)。对于多线程并发、单线程 BLAS 计算任务，为了防止线程超支 (thread oversubscription)，需要在每个并发线程中，将 BLAS 的线程数设置为 1。设置方法示例如下 (基于 rayon 的并发)：

```rust
// 半伪代码

// 主线程中无需重新设置；全局线程数已在初始化时设置
let c = dgemm(a, b);      // 使用满线程数的 BLAS
(0..niter).into_par_iter().for_each(|_| {
    // 调用 BLAS 之前设置
    omp_set_num_threads_wrapper(1);
    let c = dgemm(a, b);  // 限制并行区域内 BLAS 线程数为 1
});
let c = dgemm(a, b);      // 使用满线程数的 BLAS
```

示例中，并行区域之外前后两处 `dgemm` 均使用满线程数的 BLAS；并发线程内经 `omp_set_num_threads_wrapper(1)` 设置后，`dgemm` 仅使用单线程。`omp_set_num_threads_wrapper` 应当是 **thread-safe 与 thread-local** 的；对于 OpenMP 编译的 OpenBLAS，它即是对 `omp_set_num_threads` 的简单封装。全局线程数则由 `omp_set_num_threads_global_wrapper` 在主线程初始化时设置一次，不应在并发线程中调用。

两个线程控制函数定义于 `rest_tensors::matrix::matrix_blas_lapack`，在 crate rest 中经 `extern crate rest_tensors as tensors` 以 `tensors::matrix_blas_lapack::*` 的路径引用。两者的性质与实际调用如下：

- **线程控制函数**：`omp_set_num_threads_wrapper(n)`
  - 实际调用 (OpenBLAS 情形)： `omp_set_num_threads(n)`
  - 性质：thread-local、thread-safe；在并发线程内调用
- **线程控制函数**：`omp_set_num_threads_global_wrapper(n)`
  - 实际调用 (OpenBLAS 情形)： `openblas_set_num_threads(n)` 与 `omp_set_num_threads(n)`
  - 性质：thread-unsafe；仅主线程初始化时调用一次

## 3. 技术细节

### 3.1 其他软件的做法

首先，对于 C/C++ 程序，通常来说如果是 OpenBLAS 则会使用 OpenMP 编译的版本。OpenBLAS 自动会识别函数是否在 OpenMP 并行环境中被调用，并且会自动将线程数设置为 1[^1]；在这种情况下，程序线程是不会超支的：

```C
#pragma omp parallel for ...
{
    // OpenBLAS 会自动在并行区域将线程数设置为 1
    // 用户只需要正常地调用 dgemm(...) 即可
    dgemm(...);
}
```

[^1]: 参考程序：[common_thread.h](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/common_thread.h#L165) 与 [openblas_set_num_threads.c](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/driver/others/openblas_set_num_threads.c#L66-L68)。

对于 ORCA 等以 MPI 为主要或唯一并行方式的程序，它们通常使用串行的 BLAS，因此不需要考虑线程超支的问题。但我猜测这些程序也难以直接使用多线程的 BLAS；这会在一部分计算问题中带来代码编写上、或性能上的困难。

### 3.2 Rust 语言的困难与解决方案

Rust 语言并不使用 OpenMP。这导致的问题是，如果我们在并行子线程里调用 OpenBLAS 的多线程计算函数，OpenBLAS 并不会自动将线程数设置为 1：

```rust
(0..niter).into_par_iter().for_each(|p| { // nth 并发 rayon 线程
    let c = dgemm(a, b);                  // 每个并发再引入 nth OpenBLAS 线程
}); // 整个程序消耗 nth * nth 线程，导致线程超支
```

- 如果在并行环境中不正确地设置线程数，可能会导致线程超支 (thread oversubscription)，进而导致性能下降、或影响队列系统其他用户的作业。
- 但是如果退回到单线程的 BLAS 计算，那么并行区域外的 BLAS 计算就会被限制为单线程，性能会严重下降。

这不止是 Rust 语言的问题；只要并行框架不使用或不兼容 OpenMP，其他语言也是一样的。只是 Rust 语言的并行框架一般使用 rayon，而 rayon 与 OpenMP 不兼容。

理论上，如果我们假设被迭代的 $p$ 的维度足够大，那么一种解决方案是将并行迭代内部的线程数设置为 1；也就是说，使用 OpenBLAS 的多线程计算时，外层迭代的线程数为 $N$，而内层迭代的线程数为 1，总线程使用是 $N$。

现实里，对于 OpenBLAS 而言，确实有 `omp_set_num_threads` 可以控制 thread-local 线程。我们不确定 `omp_set_num_threads()` 是否严格线程安全，但目前的实际使用没有遇到过问题，一般也可以假设这是安全的。

这里同时指出，由于 `omp_set_num_threads` 是 thread-local 的，因此只需要在 rayon 并发线程中调用一次即可；一般不需要在并发前重新设置全局的 BLAS 线程数，也不需要在并发循环内最后、或并发结束后恢复线程数。

```rust
// 并发前不需要重新设置全局 BLAS 线程数
(0..niter).into_par_iter().for_each(|_| {
    omp_set_num_threads(1); // 并发开始时设置
    let c = dgemm(a, b);
    // 并发循环内最后不需要恢复线程数
}); // 并发结束后不需要恢复线程数
```

需要指出，这一设置按操作系统线程生效。rayon 线程池的工作线程是长驻的，在线程执行的任务闭包开头设置即可覆盖该线程后续的 BLAS 调用 (习惯上直接在每个并发任务闭包开头调用，幂等且廉价)。但对于原生线程 (如 `std::thread::scope` 启动的工作线程)，线程设置不会从启动线程继承，需在线程闭包内部、首次调用 BLAS 前设置。

### 3.3 主线程初始化线程池时使用 `openblas_set_num_threads`

一些细节参考 4.1 节 (`openblas_set_num_threads` 存在 thread-safety 问题)。我们可以在主线程中调用 `openblas_set_num_threads` 来设置全局的 BLAS 线程数，但不应该在并发线程中调用该函数。

之所以推荐使用 `openblas_set_num_threads` 函数，是因为用户有可能在环境变量中设置了 `OPENBLAS_NUM_THREADS`；但 REST 程序的输入卡 `num_threads` 的选项应该要覆盖掉环境变量的设置。如果仅使用 `omp_set_num_threads` 设置全局线程数、输入卡的 `num_threads` 比环境变量设置的 `OPENBLAS_NUM_THREADS` 更大，则 OpenBLAS 有可能会无法以预期的线程数运行或生成 thread-buffer。

实践中，这一初始化由 `omp_set_num_threads_global_wrapper` 完成，与 rayon 全局线程池的构建一起在主线程中调用一次；它会同时调用 `openblas_set_num_threads` 与 `omp_set_num_threads` (见第 2 节的对照表)。

## 4. 注意事项与已废止的设计决定

### 4.1 `openblas_set_num_threads` 存在 thread-safety 问题

函数 [`openblas_set_num_threads`](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/driver/others/blas_server_omp.c#L122-L125) 会调用[`adjust_thread_buffers`](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/driver/others/blas_server_omp.c#L77-L98)，而函数 `adjust_thread_buffers` 会对 OpenBLAS 的全局共享的 thread-buffer 进行 heap alloc/free 操作。因此，该函数是 thread-unsafe 的。函数 `goto_set_num_threads` 情况相同。

我们发现若将该函数用作并发线程控制函数，会导致对 thread-buffer 的竞争 alloc/free，从而以非常小的概率影响其他线程正在进行的 BLAS 计算 (如 dgemm)，最终导致数值不稳定。但该函数需要在主线程中调用一次，以设置全局的 BLAS 线程数与 thread-buffer。

### 4.2 已废止：pthreads 并行的 OpenBLAS 支持

REST 目前要求使用 OpenMP 编译的 OpenBLAS (见决定)。但 REST 在早年曾同时允许 pthreads 与 OpenMP 并行版本的 OpenBLAS。实际上，在本 ADR 确立当前方案之前，同时使用 `openblas_set_num_threads` 与 `omp_set_num_threads` 是为了同时保证 pthreads 与 OpenMP 并行版本的 OpenBLAS 都能有效地并行运行：pthreads 并不支持函数 `omp_set_num_threads` 的设置，因此必须也只能使用 `openblas_set_num_threads` 来设置线程数。

**下述伪代码的做法对于 pthreads 的 OpenBLAS 是有效的**。这是因为 pthreads 并行的 OpenBLAS 不会在 `openblas_set_num_threads` 中调用 `adjust_thread_buffers`。但需要留意，`openblas_set_num_threads` 是 thread-safe 但并非 thread-local 的，意味着并行线程内部的改动是全局改动，会影响主线程的 BLAS 线程数设置。

```python
# 伪代码
num_threads_orig = openblas_get_num_threads() # 获取全局线程数
openblas_set_num_threads(1)        # 设置全局线程数为 1
parallel for tid in 0..NTHREADS:
    # openblas_set_num_threads(1)  # 线程数设置放在内部也是有效的
    C = dgemm(A, B)                # 执行 BLAS
openblas_set_num_threads(num_threads_orig) # 恢复全局线程数
```

建议停止 pthreads 的 OpenBLAS 支持，出于两方面考虑。
- 效率，对程序使用者而言这通常是更重要的问题：pthreads 构建的 OpenBLAS 性能一般要比 OpenMP 差一些，特别是在较大的矩阵问题上。
- 线程控制：RSTSR 的并行线程控制依赖于 thread-local 的线程控制；OpenBLAS 只有 OpenMP 版本的 `omp_set_num_threads` 是 thread-local 的，而 pthreads 版本仅有全局的 `openblas_set_num_threads`。

需要说明的是，这并非认为 pthreads 技术本身劣于 OpenMP，而是 OpenBLAS 的 OpenMP 在上述两方面更符合 REST 的需要。

如果线程控制函数只能保证 thread-safe 而不能保证 thread-local (即只能作全局设置)，则需要在并发前保存线程数、并发后恢复线程数；上文 pthreads 情形的伪代码正是这一模式。该模式在 thread-local 的情形也适用，但冗余。

### 4.3 部署验证：确认 OpenBLAS 为 OpenMP 构建

本 ADR 的线程控制方案依赖 OpenMP 构建的 OpenBLAS。若实际链接到 pthreads 构建的 OpenBLAS，`omp_set_num_threads_wrapper` 的设置将静默失效 (仅 `omp_set_num_threads_global_wrapper` 中的全局设置仍然有效)：程序一般仍能运行，但并发批量 BLAS 计算会发生线程超支，即退化为 4.2 节的情形。

部署时应验证所链接的 OpenBLAS 的并行方式。规范的方法是用 `ldd` (Linux) 或 `otool -L` (macOS) 检查库的动态依赖：OpenMP 构建的 OpenBLAS 会依赖 OpenMP 运行时 (如 `libgomp`)，pthreads 构建则没有此类依赖。也可以在程序中调用 `openblas_get_parallel()` 判断 (返回 2 为 OpenMP，1 为 pthreads)。注意不要依赖库文件名判断：文件名中的字母 (如 `libopenblasp` 中的 `p`) 在不同发行方式下含义并不一致，未必表示 pthreads。

## 5. RSTSR 的处理模式

RSTSR 数学库从设计上，希望库使用者不需要手动调整线程数设置。RSTSR 的库函数内部会判断是否在 rayon 并发线程中被调用 (基于 `rayon::current_thread_index`)。如果是，则自动将 BLAS 线程数设置为 1；如果不是，则使用全局的 BLAS 线程数。

目前的自动限制覆盖矩阵乘法 (matmul) 与基于 TBLIS 的 einsum。RSTSR 的 LAPACK 类驱动 (如 `eigh`) 尚未作此处理；这是 RSTSR 已知的待修复问题，在修复之前，这些函数应在主线程调用，或在并发线程内手动调用 `omp_set_num_threads_wrapper(1)`。

为达到这种设计目的，有必要使用 thread-local 的线程控制。RSTSR 在运行时通过 `openblas_get_parallel()` 识别 OpenBLAS 的并行方式并选择线程控制函数：OpenMP 构建使用 `omp_set_num_threads` (thread-local)；pthreads 构建则退化为 `openblas_set_num_threads` (全局的设置与恢复，即 4.2 节的模式)，无法兑现 thread-local 的语义。因此，对于下述用法，REST 规范要求 OpenMP 并行的 OpenBLAS。

```rust
// a, b, c are RSTSR tensor/tensorview
let c = &a % &b;      // 使用满线程数的 BLAS
(0..niter).into_par_iter().for_each(|_| {
    let c = &a % &b;  // 自动限制并行区域内 BLAS 线程数为 1
});
let c = &a % &b;      // 使用满线程数的 BLAS
```

## 6. 潜在影响

由于在线程内部更改了 OpenMP 的线程数设置，其他依赖 OpenMP 的库同时会受到影响。一般来说影响不大：
- 当前的线程控制策略对其他 OpenMP 库仍然是有效且有意义的；
- 目前 REST 所依赖的外部库中，只有 OpenBLAS 使用了 OpenMP，且其他库目前没有在并发线程中调用 OpenMP 的函数；
但若未来需要引入其他依赖 OpenMP 的库，当前的策略可能需要调整。

## 7. 备注 (非规范性)

OpenMP 构建的 OpenBLAS 是 REST 唯一的规范后端；本节内容仅为未来可能支持其他数学库时提供参考，不构成支持承诺。

其他主流的数学库 (MKL, BLIS/AOCL, KML 等) 有自己特有的 thread-local 线程控制函数。RSTSR 通过 crate feature 已经支持部分下述数学库；如果 REST 本体需要支持，可以选择使用下述函数作为替代。另外，rest_tensors 的 `intel-mkl` feature 目前将线程控制函数映射到 `mkl_set_num_threads` (全局设置)；MKL 支持为实验性质。

| 数学库 | thread-local setter | thread-local getter |
| --- | --- | --- |
| OpenBLAS-pthread | No | No |
| OpenBLAS-OpenMP | `omp_set_num_threads` | `omp_get_max_threads` |
| MKL (Intel) | `MKL_Set_Num_Threads_Local` | `MKL_Get_Max_Threads` |
| BLIS/AOCL (AMD) | `bli_thread_set_num_threads` | `bli_thread_get_num_threads` |
| KML (Huawei) | `BlasSetNumThreadsLocal`, `KmlSetNumThreads` | `BlasGetNumThreadsLocal` |

本决定所述机制的最终形式、含全局线程初始化，随 rest_tensors PR #13 同时给出。
