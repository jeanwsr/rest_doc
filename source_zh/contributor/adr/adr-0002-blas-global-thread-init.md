# ADR-0002：全局 BLAS 线程数初始化仅使用 omp_set_num_threads

- **状态**：使用中
- **日期**：2026-09-16
- **适用范围**：主线程的全局 BLAS 线程数初始化：crate `rest` 的程序初始化流程，以及 crate `rest_tensors` 中 `omp_set_num_threads_global_wrapper` 的实现。本篇取代 [ADR-0001](adr-0001-blas-threads) 第 3.3 节；ADR-0001 的其余决定 (并发线程内的 thread-local 线程控制、RSTSR 的处理模式等) 继续有效。
- **作者**：GLM-5.3 (依王石嵘与祝震予的测试与分析结论起草)
- **审阅**：祝震予 (ajz34@outlook.com)
- **取代**：ADR-0001 第 3.3 节

## 决定

- 主线程初始化时，全局 BLAS 线程数仅经 `omp_set_num_threads` 设置，不再调用 `openblas_set_num_threads`。
- `omp_set_num_threads_global_wrapper` 的实现相应改为仅封装 `omp_set_num_threads`。它与 `omp_set_num_threads_wrapper` 的目的不同 (主线程的一次性全局初始化，而非并发线程内的 thread-local 设置)，但实现恰巧相同；两者仍作为两个函数分别保留。
- 并发线程内的线程控制不变：在调用 BLAS 前设置 `omp_set_num_threads_wrapper(1)`，仍遵循 [ADR-0001](adr-0001-blas-threads)。

## 1. 背景与动机

[ADR-0001](adr-0001-blas-threads) 第 3.3 节的原方案，是主线程初始化时同时调用 `openblas_set_num_threads` 与 `omp_set_num_threads`。当时推荐 `openblas_set_num_threads` 的原因是：用户有可能在环境变量中设置了 `OPENBLAS_NUM_THREADS`，而 REST 程序的输入卡 `num_threads` 的选项应该要覆盖掉环境变量的设置。如果仅使用 `omp_set_num_threads` 设置全局线程数、输入卡的 `num_threads` 比环境变量设置的 `OPENBLAS_NUM_THREADS` 更大，则 OpenBLAS 有可能会无法以预期的线程数运行或生成 thread-buffer。

OpenBLAS v0.3.34 更新了线程控制机制。在程序不使用 C/Fortran 的 OpenMP 时 (REST 以 rayon 作并行，即属此情形)，`openblas_set_num_threads` 会在所有线程以全局线程数执行程序，造成死锁或线程超支。这出于正确性的强制要求，原先以 `openblas_set_num_threads` 兑现环境变量覆盖的做法随之不再可用；主线程初始化改为仅使用 `omp_set_num_threads` 设置全局的 BLAS 线程数。同时测试的情况认为，用户环境变量里设置的 `OPENBLAS_NUM_THREADS` 并不会影响程序的正确性 (缓存仍然会依程序可见的实际的线程数生成)，因此不再需要 `openblas_set_num_threads`。

详细的分析与复现，见由王石嵘与祝震予提供的[测试与分析结论](https://gitee.com/ajz34/openblas-deadlock/tree/icv-only-mwe/)。

## 2. 技术细节

`omp_set_num_threads` 是 thread-local 的，也就意味着主线程的 BLAS 线程数设置不会被并发线程内的设置覆盖。因此，该函数不仅可以在并发线程内调用，也可以在主线程中调用 (如在主线程初始化时设置全局的 BLAS 线程数)。

## 3. 备注 (非规范性)

- 本篇的废止决定出于正确性 (死锁、线程超支) 的强制要求，与 `openblas_set_num_threads` 的 thread-safety 问题 (ADR-0001 第 4.1 节) 无关，后者在旧版本中亦是已知问题。
- 被取代的旧机制在 OpenBLAS v0.3.33 或更早版本适用；ADR-0001 第 3.3 节的内容留作历史记录。
