# ADR-0001: Parallel BLAS thread scheduling for OpenBLAS

- **Status**: In use
- **Date**: 2026-08-13
- **Partially superseded**: Section 3.3 (using `openblas_set_num_threads` for main-thread initialization) has been superseded by ADR-0002 (2026-09-16); all other decisions remain in force
- **Scope**: All modules of REST proper (crate `rest`) — SCF, post-SCF, GW, perturbation, TD, etc. — not limited to any particular module. The thread-control functions are implemented in crate `rest_tensors`; the handling for RSTSR tensor operations is covered in Section 5. Legacy code partially does not comply yet and will be aligned gradually during refactoring.
- **Author**: Zhu Zhenyu (ajz34@outlook.com)

> **Update notice (2026-09-16)**: The global thread-count initialization mechanism described in Section 3.3 of this record (calling both `openblas_set_num_threads` and `omp_set_num_threads` via `omp_set_num_threads_global_wrapper` during main-thread initialization) has been superseded by [ADR-0002](adr-0002-blas-global-thread-init); all other decisions remain in force.

## Decision

- REST currently uses an **OpenMP-compiled OpenBLAS** as the BLAS backend; this ADR constrains only that case. Other math libraries are currently supported experimentally (see the remarks in Section 7).
- Before calling BLAS inside concurrent (rayon) threads, call `omp_set_num_threads_wrapper(1)`; this function is thread-safe and thread-local, so there is no need to save/restore the thread count.
- When using RSTSR matrix multiplications and einsums, no manual thread-count setting is needed (see Section 5).
- During main-thread initialization, together with the construction of the rayon global thread pool, call `omp_set_num_threads_global_wrapper` once to set the global BLAS thread count; the `num_threads` option of the input file overrides settings of environment variables such as `OPENBLAS_NUM_THREADS` (see Section 3.3; that section has been superseded by ADR-0002).
- Support for pthreads-parallel OpenBLAS is deprecated (see Section 4.2).

The core pattern of the decision:

```rust
(0..niter).into_par_iter().for_each(|_| {
    omp_set_num_threads_wrapper(1); // inside concurrent threads: single-threaded BLAS
    let c = dgemm(a, b);
});
```

## 1. Design motivation

Electronic-structure calculations use BLAS functions extensively. REST currently uses OpenBLAS as its principal BLAS library.

There are two situations for BLAS calls. The simpler one is an ordinary matrix multiplication or matrix factorization; simply calling it, BLAS runs in parallel with the globally configured full thread count.

The other situation can be reduced to batched matrix multiplication $C_{i j}^p = \sum_k A_{i k}^p B_{k j}^p$ or operator fusion built upon it. When the dimension of index $p$ is large enough, we prefer to iterate over $p$ in parallel rather than perform an ordinary matrix multiplication:
```python
# pseudo-code
parallel for p in 0..N {
    C[p] = A[p] @ B[p]
}
```
In more involved problems, the iteration over $p$ involves multiple steps of matrix multiplication or even complex operator fusion. If the iteration over $p$ is not parallelized, either performance degrades or intermediate tensors of large memory footprint are produced.

We want both situations to use multi-core CPUs effectively. A multithreaded OpenBLAS based on OpenMP parallelism supports both usages exactly: ordinary BLAS computations use the full thread count, while batched concurrent computations restrict BLAS to a single thread within each concurrent thread. How to schedule this is the content of this decision.

## 2. Usage

In the following, `dgemm` stands for any multithreaded BLAS/LAPACK call. Outside rayon parallel regions, no setting is needed: BLAS runs with the global thread count configured at initialization (see Section 3.3). For workloads with multi-threaded concurrency and single-threaded BLAS computation, the BLAS thread count must be set to 1 in each concurrent thread to prevent thread oversubscription. An example of the setting (based on rayon concurrency):

```rust
// semi-pseudo-code

// No resetting on the main thread; the global thread count was set at initialization
let c = dgemm(a, b);      // BLAS with the full thread count
(0..niter).into_par_iter().for_each(|_| {
    // Set before calling BLAS
    omp_set_num_threads_wrapper(1);
    let c = dgemm(a, b);  // Restrict the BLAS thread count to 1 inside the parallel region
});
let c = dgemm(a, b);      // BLAS with the full thread count
```

In the example, the two `dgemm` calls before and after the parallel region both use BLAS with the full thread count; inside concurrent threads, after the `omp_set_num_threads_wrapper(1)` setting, `dgemm` uses only a single thread. `omp_set_num_threads_wrapper` is expected to be **thread-safe and thread-local**; for OpenMP-compiled OpenBLAS it is a thin wrapper over `omp_set_num_threads`. The global thread count, in turn, is set once during main-thread initialization by `omp_set_num_threads_global_wrapper`, which must not be called from concurrent threads.

The two thread-control functions are defined in `rest_tensors::matrix::matrix_blas_lapack` and referenced from crate rest via `extern crate rest_tensors as tensors` under the path `tensors::matrix_blas_lapack::*`. Their properties and actual calls are as follows:

- **Thread-control function**: `omp_set_num_threads_wrapper(n)`
  - Actual call (OpenBLAS case): `omp_set_num_threads(n)`
  - Properties: thread-local, thread-safe; called inside concurrent threads
- **Thread-control function**: `omp_set_num_threads_global_wrapper(n)`
  - Actual call (OpenBLAS case): `openblas_set_num_threads(n)` and `omp_set_num_threads(n)`
  - Properties: thread-unsafe; called once during main-thread initialization

## 3. Technical details

### 3.1 How other software does it

First, for C/C++ programs, the OpenMP-compiled version of OpenBLAS is normally used. OpenBLAS automatically detects whether a function is called inside an OpenMP parallel region and automatically sets the thread count to 1 there[^1]; in that case the program threads are not oversubscribed:

```C
#pragma omp parallel for ...
{
    // OpenBLAS automatically sets the thread count to 1 inside parallel regions
    // The user simply calls dgemm(...) as usual
    dgemm(...);
}
```

[^1]: Referenced source: [common_thread.h](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/common_thread.h#L165) and [openblas_set_num_threads.c](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/driver/others/openblas_set_num_threads.c#L66-L68).

For programs whose primary or only parallel mode is MPI, such as ORCA, serial BLAS is normally used, so thread oversubscription is not a concern. But I conjecture that these programs also have difficulty using multithreaded BLAS directly; this creates difficulties in code writing or in performance for part of the computational problems.

### 3.2 Difficulties of the Rust language and their solution

The Rust language does not use OpenMP. The resulting problem is that if we call multithreaded computation functions of OpenBLAS inside parallel worker threads, OpenBLAS does not automatically set the thread count to 1:

```rust
(0..niter).into_par_iter().for_each(|p| { // nth concurrent rayon threads
    let c = dgemm(a, b);                  // each concurrency introduces nth OpenBLAS threads again
}); // the program consumes nth * nth threads in total, causing thread oversubscription
```

- If the thread count is not set correctly inside parallel regions, thread oversubscription may occur, which degrades performance or disturbs jobs of other users on queueing systems.
- But if we fell back to single-threaded BLAS computation entirely, BLAS computation outside parallel regions would be restricted to a single thread, and performance would degrade severely.

This is not a problem of the Rust language alone; it applies to any language whose parallel framework does not use, or is incompatible with, OpenMP. It is just that parallel frameworks in Rust generally use rayon, and rayon is incompatible with OpenMP.

In theory, if we assume the iterated dimension of $p$ is large enough, one solution is to set the thread count inside the parallel iteration to 1; that is, when using multithreaded computation of OpenBLAS, the outer iteration runs with $N$ threads, the inner with 1, and the total thread usage is $N$.

In practice, OpenBLAS does have `omp_set_num_threads` to control a thread-local count. We are not certain that `omp_set_num_threads()` is strictly thread-safe, but no problems have been encountered in actual use so far, and it is generally safe to assume so.

It should also be pointed out that since `omp_set_num_threads` is thread-local, calling it once in a rayon concurrent thread is enough; in general there is no need to reconfigure the global BLAS thread count before entering concurrency, nor to restore the thread count at the end of the concurrent loop or after concurrency ends.

```rust
// No need to reconfigure the global BLAS thread count before concurrency
(0..niter).into_par_iter().for_each(|_| {
    omp_set_num_threads(1); // Set at the start of the concurrent closure
    let c = dgemm(a, b);
    // No need to restore the thread count at the end of the concurrent loop
}); // No need to restore the thread count after concurrency ends
```

Note that this setting takes effect per operating-system thread. The worker threads of the rayon pool are long-lived: setting it at the start of a task closure executed by the thread covers the thread's subsequent BLAS calls (by convention it is called at the start of every concurrent task closure — idempotent and cheap). For native threads, however (e.g. workers started by `std::thread::scope`), the setting is not inherited from the spawning thread; it must be set inside the thread closure, before the first BLAS call.

### 3.3 Using `openblas_set_num_threads` when initializing the thread pool on the main thread

The mechanism described in this section was superseded by [ADR-0002](adr-0002-blas-global-thread-init) on 2026-09-16; the content is kept as a historical record.

For some details, refer to Section 4.1 (`openblas_set_num_threads` has a thread-safety problem). We may call `openblas_set_num_threads` on the main thread to set the global BLAS thread count, but must not call this function from concurrent threads.

The reason for recommending `openblas_set_num_threads` is that the user may have set `OPENBLAS_NUM_THREADS` in the environment, while the `num_threads` option of the REST input file should override the setting of environment variables. If the global thread count is set only via `omp_set_num_threads`, and the input-file `num_threads` is larger than the `OPENBLAS_NUM_THREADS` set by the environment, OpenBLAS might fail to run with the expected thread count or to generate thread-buffers.

In practice, this initialization is performed by `omp_set_num_threads_global_wrapper`, called once on the main thread together with the construction of the rayon global thread pool; it calls both `openblas_set_num_threads` and `omp_set_num_threads` (see the comparison table in Section 2).

## 4. Caveats and deprecated design decisions

### 4.1 `openblas_set_num_threads` has a thread-safety problem

The function [`openblas_set_num_threads`](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/driver/others/blas_server_omp.c#L122-L125) calls [`adjust_thread_buffers`](https://github.com/OpenMathLib/OpenBLAS/blob/8d73a856fb6727f3593b3783c51353b5a9c1eacd/driver/others/blas_server_omp.c#L77-L98), and `adjust_thread_buffers` performs heap alloc/free operations on OpenBLAS's globally shared thread buffers. The function is therefore thread-unsafe. The same applies to the function `goto_set_num_threads`.

We found that using this function as the concurrent-thread control function leads to competing alloc/free of thread buffers, which with very small probability disturbs BLAS computations (such as dgemm) in progress on other threads, eventually causing numerical instability. However, this function needed to be called once on the main thread, in order to set the global BLAS thread count and the thread buffers.

### 4.2 Deprecated: support for pthreads-parallel OpenBLAS

REST currently requires an OpenMP-compiled OpenBLAS (see the Decision). But in earlier years, REST allowed both the pthreads and the OpenMP parallel versions of OpenBLAS. In fact, before this ADR established the current scheme, using `openblas_set_num_threads` and `omp_set_num_threads` together was meant to keep both the pthreads and the OpenMP parallel versions of OpenBLAS running effectively in parallel: pthreads does not honor the setting of `omp_set_num_threads`, so `openblas_set_num_threads` was the only way to set the thread count.

**The practice of the pseudo-code below is valid for a pthreads OpenBLAS**. This is because a pthreads-parallel OpenBLAS does not call `adjust_thread_buffers` inside `openblas_set_num_threads`. Note, however, that `openblas_set_num_threads` is thread-safe but not thread-local, meaning that changes inside parallel threads are global changes that affect the main thread's BLAS thread-count setting.

```python
# pseudo-code
num_threads_orig = openblas_get_num_threads() # Get the global thread count
openblas_set_num_threads(1)        # Set the global thread count to 1
parallel for tid in 0..NTHREADS:
    # openblas_set_num_threads(1)  # Setting the thread count inside is also valid
    C = dgemm(A, B)                # Perform BLAS
openblas_set_num_threads(num_threads_orig) # Restore the global thread count
```

We suggest discontinuing the support of pthreads OpenBLAS, for two considerations.
- Efficiency, which is usually the more important issue for program users: a pthreads-built OpenBLAS generally performs somewhat worse than OpenMP, especially for larger matrix problems.
- Thread control: RSTSR's parallel thread control relies on thread-local thread control; only the OpenMP version of `omp_set_num_threads` in OpenBLAS is thread-local, while the pthreads version only offers the global `openblas_set_num_threads`.

It should be noted that this is not a judgment that pthreads as a technique is inferior to OpenMP; rather, OpenBLAS's OpenMP fits REST's needs better in the two respects above.

If a thread-control function can only guarantee thread-safety but not thread-locality (i.e. only global setting), the thread count must be saved before concurrency and restored after it; the pseudo-code for the pthreads case above is exactly this pattern. The pattern also applies in the thread-local case, but is redundant.

### 4.3 Deployment verification: confirming that OpenBLAS is built with OpenMP

The thread-control scheme of this ADR depends on an OpenMP-built OpenBLAS. If a pthreads-built OpenBLAS is actually linked, the setting of `omp_set_num_threads_wrapper` silently has no effect (only the global setting inside `omp_set_num_threads_global_wrapper` still works): the program generally still runs, but concurrent batched BLAS computation suffers thread oversubscription, i.e. it degenerates to the situation of Section 4.2.

At deployment time, the parallel mode of the linked OpenBLAS should be verified. The canonical method is to inspect the dynamic dependencies of the library with `ldd` (Linux) or `otool -L` (macOS): an OpenMP-built OpenBLAS depends on an OpenMP runtime (such as `libgomp`), whereas a pthreads build has no such dependency. Alternatively, `openblas_get_parallel()` can be called from the program (returning 2 for OpenMP, 1 for pthreads). Do not rely on the library file name: the meaning of letters in file names (such as the `p` in `libopenblasp`) is inconsistent across distributions and does not necessarily indicate pthreads.

## 5. The RSTSR approach

The RSTSR math library is designed so that its users do not need to adjust thread-count settings manually. RSTSR library functions detect internally whether they are called from a rayon concurrent thread (based on `rayon::current_thread_index`). If so, the BLAS thread count is automatically set to 1; if not, the global BLAS thread count is used.

Currently this automatic restriction covers matrix multiplication (matmul) and TBLIS-based einsum. RSTSR's LAPACK-type drivers (such as `eigh`) do not yet handle this; this is a known issue of RSTSR to be fixed — until then, these functions should be called on the main thread, or `omp_set_num_threads_wrapper(1)` should be called manually inside concurrent threads.

To achieve this design goal, thread-local thread control is necessary. At runtime, RSTSR identifies the parallel mode of OpenBLAS via `openblas_get_parallel()` and selects the thread-control function: for the OpenMP build it uses `omp_set_num_threads` (thread-local); for the pthreads build it degenerates to `openblas_set_num_threads` (global setting and restoration, i.e. the pattern of Section 4.2) and cannot deliver the thread-local semantics. Therefore, for the usage below, REST normatively requires OpenMP-parallel OpenBLAS.

```rust
// a, b, c are RSTSR tensor/tensorview
let c = &a % &b;      // BLAS with the full thread count
(0..niter).into_par_iter().for_each(|_| {
    let c = &a % &b;  // Automatically restrict the BLAS thread count to 1 inside the parallel region
});
let c = &a % &b;      // BLAS with the full thread count
```

## 6. Potential impact

Since the OpenMP thread-count setting is changed inside threads, other libraries that depend on OpenMP are affected as well. Generally the impact is small:
- The current thread-control strategy remains valid and meaningful for other OpenMP libraries;
- Among the external libraries REST currently depends on, only OpenBLAS uses OpenMP, and no other library currently calls OpenMP functions from concurrent threads;
but if other OpenMP-dependent libraries need to be introduced in the future, the current strategy may have to be adjusted.

## 7. Remarks (non-normative)

OpenMP-built OpenBLAS is the only normative backend of REST; this section only provides reference for possible future support of other math libraries and constitutes no commitment of support.

Other mainstream math libraries (MKL, BLIS/AOCL, KML, etc.) have their own thread-local thread-control functions. RSTSR already supports some of the libraries below through crate features; if REST proper needs to support them, the functions below can be used as replacements. Note that the `intel-mkl` feature of rest_tensors currently maps the thread-control functions to `mkl_set_num_threads` (a global setting); MKL support is experimental.

| Math library | thread-local setter | thread-local getter |
| --- | --- | --- |
| OpenBLAS-pthread | No | No |
| OpenBLAS-OpenMP | `omp_set_num_threads` | `omp_get_max_threads` |
| MKL (Intel) | `MKL_Set_Num_Threads_Local` | `MKL_Get_Max_Threads` |
| BLIS/AOCL (AMD) | `bli_thread_set_num_threads` | `bli_thread_get_num_threads` |
| KML (Huawei) | `BlasSetNumThreadsLocal`, `KmlSetNumThreads` | `BlasGetNumThreadsLocal` |

The final form of the mechanism described in this decision, including the global thread initialization, is given together with rest_tensors PR #13.
