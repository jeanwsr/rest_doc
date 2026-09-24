# ADR-0002: Global BLAS thread-count initialization uses only omp_set_num_threads

- **Status**: In use
- **Date**: 2026-09-16
- **Scope**: Main-thread global BLAS thread-count initialization: the program initialization flow of crate `rest`, and the implementation of `omp_set_num_threads_global_wrapper` in crate `rest_tensors`. This record supersedes Section 3.3 of [ADR-0001](adr-0001-blas-threads); all other decisions of ADR-0001 (thread-local thread control inside concurrent threads, the RSTSR approach, etc.) remain in force.
- **Author**: GLM-5.3 (drafted from the tests and analysis by Wang Shirong and Zhu Zhenyu)
- **Reviewed by**: Zhu Zhenyu (ajz34@outlook.com)
- **Supersedes**: Section 3.3 of ADR-0001

## Decision

- During main-thread initialization, the global BLAS thread count is set solely via `omp_set_num_threads`; `openblas_set_num_threads` is no longer called.
- The implementation of `omp_set_num_threads_global_wrapper` accordingly becomes a wrapper over `omp_set_num_threads` only. Its purpose differs from that of `omp_set_num_threads_wrapper` (one-time global initialization on the main thread, rather than a thread-local setting inside concurrent threads), but the implementations happen to be identical; the two remain separate functions.
- Thread control inside concurrent threads is unchanged: set `omp_set_num_threads_wrapper(1)` before calling BLAS, still following [ADR-0001](adr-0001-blas-threads).

## 1. Background and motivation

The original scheme of Section 3.3 of [ADR-0001](adr-0001-blas-threads) was to call both `openblas_set_num_threads` and `omp_set_num_threads` during main-thread initialization. `openblas_set_num_threads` was recommended at the time because the user may have set `OPENBLAS_NUM_THREADS` in the environment, while the `num_threads` option of the REST input file should override the setting of environment variables. If the global thread count were set only via `omp_set_num_threads`, and the input-file `num_threads` is larger than the `OPENBLAS_NUM_THREADS` set by the environment, OpenBLAS might fail to run with the expected thread count or to generate thread-buffers.

OpenBLAS v0.3.34 changed its thread-control mechanism. When the program does not use C/Fortran OpenMP (REST parallelizes in Rust with rayon, which is exactly this case), `openblas_set_num_threads` makes all threads execute the program with the global thread count, causing deadlock or thread oversubscription. This is mandated by correctness, and the old practice of honoring the environment-variable override via `openblas_set_num_threads` is no longer available; main-thread initialization now sets the global BLAS thread count using `omp_set_num_threads` alone. Testing further indicates that an `OPENBLAS_NUM_THREADS` set by the user in the environment does not compromise correctness (the buffers are still generated according to the number of threads actually visible to the program), so `openblas_set_num_threads` is no longer needed.

For detailed analysis and reproduction, see the [tests and analysis](https://gitee.com/ajz34/openblas-deadlock/tree/icv-only-mwe/) provided by Wang Shirong and Zhu Zhenyu.

## 2. Technical details

`omp_set_num_threads` is thread-local, which means that the main thread's BLAS thread-count setting is not overridden by settings inside concurrent threads. The function can therefore be called not only inside concurrent threads, but also on the main thread (e.g. to set the global BLAS thread count during main-thread initialization).

## 3. Remarks (non-normative)

- The deprecation decision of this record is mandated by correctness (deadlock, thread oversubscription), and is unrelated to the thread-safety problem of `openblas_set_num_threads` (Section 4.1 of ADR-0001), which was likewise a known issue in older versions.
- The superseded old mechanism applies to OpenBLAS v0.3.33 or earlier versions; the content of Section 3.3 of ADR-0001 is kept as a historical record.
