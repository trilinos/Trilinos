API: Batched Dense (DLA)
########################

.. toctree::
   :maxdepth: 2
   :hidden:

   dense/batched_iamax
   dense/batched_trsv
   dense/batched_tbsv
   dense/batched_ger
   dense/batched_syr
   dense/batched_getrf
   dense/batched_getrs
   dense/batched_gbtrf
   dense/batched_gbtrs
   dense/batched_pbtrf
   dense/batched_pbtrs
   dense/batched_pttrf
   dense/batched_pttrs
   dense/batched_lacgv
   dense/batched_laswp
   dense/batched_apply_householder
   dense/batched_householder

Our implementation of batched dense linear algebra (DLA) allows user to compose various batched DLA operations.
For example, consider a case where small element matrices are created via `gemm` and those matrices are triangular solved by `lu` and `trsv`.
An approach using conventional batched BLAS interfaces would be like the following pseudo code:

.. code-block:: c++

  int N = 1000; /// # of batch items
  int m = 8;    /// a square matrix size

  Kokkos::View<double***> AA("AA", N, m, m); /// element matrices
  Kokkos::View<double**>  LL("LL", N, m);    /// left basis vectors
  Kokkos::View<double**>  RR("RR", N, m);    /// right basis vectors
  Kokkos::View<double**>  BB("BB", N, m);    /// load vector and would be overwritten by a solution

  batched_dgemm(LL, RR, AA);     /// construct element matrices via batched dgemm
  batched_dgetrf(AA);            /// perform lu decomposition of each instance of A matrix array
  batched_dtrsv("Lower", AA, BB) /// perform forward substitution
  batched_dtrsv("Upper", AA, BB) /// perform backward substitution

Clearly, a performance problem of the above code comes from the fact that the sequence of batched functions does not exploit temporal locality between DLA functions;
each batched function sweeps the entire set of batched matrices or vectors in parallel for a single DLA operation.
On the other hand, Kokkos batched APIs provide functor-level interfaces so that a user can compose a new batched function.
The following example has the same functionality as the above:

.. code-block:: c++

  int N = 1000; /// # of batch items
  int m = 8;    /// a square matrix size

  Kokkos::View<double***> AA("AA", N, m, m); /// element matrices
  Kokkos::View<double**>  LL("LL", N, m);    /// left basis vectors
  Kokkos::View<double**>  RR("RR", N, m);    /// right basis vectors
  Kokkos::View<double**>  BB("BB", N, m);    /// load vector and would be overwritten by a solution

  using namespace KokkosBatched;
  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i) {
    auto A = Kokkos::subview(AA, i, Kokkos::ALL(), Kokkos::ALL()); /// ith matrix
    auto L = Kokkos::subview(LL, i, Kokkos::ALL());                /// ith left vector
    auto R = Kokkos::subview(RR, i, Kokkos::ALL());                /// ith right vector
    auto B = Kokkos::subview(BB, i, Kokkos::ALL());                /// ith load/solution vector

    SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Unblocked>
      ::invoke(one, L, R, zero, A);
    SerialLU<Algo::LU::Unblocked>
      ::invoke(A);
    SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::UnitDiag,Algo::Trsv::Unblocked>
      ::invoke(one, A, B);
  });

In this example, a single `parallel_for` is launched to compute a sequence of DLA operations i.e., `gemm`, `lu` and `trsv`.
Then, one may ask, "How is this different from putting BLAS and LAPACK functions inside a parallel loop ?".
The main difference is that Kokkos batched APIs are very lightweight generic implementations focusing on small matrix sizes (kernels are developed and tuned from the application context).
Most vendor provided DLA libraries such as Intel MKL and NVIDIA CUBLAS perform well for large problem sizes. For tiny problem sizes like 3x3 or 5x5 matrix problems,
it is not feasible to use vendor optimized DLA libraries as even error checking already puts a quite amount of overhead for such tiny problem sizes.
Furthermore, CUBLAS (or GPU numeric libraries) cannot be nested inside of a parallel region, which is why CUBLAS provides separate batched APIs.
Kokkos batched APIs provide generic header-only functor-level implementations that can be embedded in a parallel region.
The APIs can be mapped to Kokkos hierarchical parallelism and also provide various implementations of algorithms with a template argument,
which allows users to choose (or customize) the batched routines for their application context.

BLAS support
============

Below are tables summarizing the currently supported function calls in Kokkos Kernels.

BLAS 1
------

.. list-table::
   :widths: 12 10 10 10
   :header-rows: 1

   * - BLAS Call
     - Serial
     - Team
     - TeamVector
   * - ROTG
     - --
     - --
     - --
   * - ROTMG
     - --
     - --
     - --
   * - ROT
     - --
     - --
     - --
   * - ROTM
     - --
     - --
     - --
   * - SWAP
     - `kokkos_swap(a, b) <https://kokkos.org/kokkos-core-wiki/API/core/utilities/swap.html>`_
     - --
     - --
   * - SCAL
     - `Blas::SerialScale`
     - `Blas::TeamScale`
     - `Blas::TeamVectorScale`
   * - COPY
     - `SerialCopy`
     - `TeamCopy`
     - `TeamVectorCopy`
   * - AXPY
     - `SerialAxpy`
     - `TeamAxpy`
     - `TeamVectorAxpy`
   * - DOT*
     - `SerialDot`
     - `TeamDot`
     - `TeamVectorDot`
   * - DOTU
     - --
     - --
     - --
   * - DOTC*
     - --
     - --
     - --
   * - NRM2
     - --
     - --
     - --
   * - ASUM
     - --
     - --
     - --
   * - IAMAX
     - :doc:`SerialIamax <dense/batched_iamax>`
     - --
     - --

BLAS 2
------

.. list-table::
   :widths: 12 10 10 10
   :header-rows: 1

   * - BLAS Call
     - Serial
     - Team
     - TeamVector
   * - GEMV
     - `Blas::SerialGemv`
     - `TeamGemv`
     - `TeamVectorGemv`
   * - GBMV
     - --
     - --
     - --
   * - SYMV
     - --
     - --
     - --
   * - SBMV
     - --
     - --
     - --
   * - SPMV
     - --
     - --
     - --
   * - TRMV
     - --
     - --
     - --
   * - TBMV
     - --
     - --
     - --
   * - TPMV
     - --
     - --
     - --
   * - TRSV
     - :doc:`SerialTrsv <dense/batched_trsv>`
     - :doc:`TeamTrsv <dense/batched_trsv>`
     - :doc:`TeamVectorTrsv <dense/batched_trsv>`
   * - TBSV
     - :doc:`SerialTbsv <dense/batched_tbsv>`
     - --
     - --
   * - TPSV
     - --
     - --
     - --
   * - GER
     - :doc:`SerialGer <dense/batched_ger>`
     - --
     - --
   * - GERU
     - :doc:`SerialGer <dense/batched_ger>`
     - --
     - --
   * - GERC
     - :doc:`SerialGer <dense/batched_ger>`
     - --
     - --
   * - SYR
     - :doc:`SerialSyr <dense/batched_syr>`
     - --
     - --
   * - HER
     - :doc:`SerialSyr <dense/batched_syr>`
     - --
     - --
   * - SPR
     - --
     - --
     - --
   * - SYR2
     - --
     - --
     - --

BLAS 3
------

.. list-table::
   :widths: 12 10 10 10
   :header-rows: 1

   * - BLAS Call
     - Serial
     - Team
     - TeamVector
   * - GEMM
     - `SerialGemm`
     - `TeamGemm`
     - `TeamVectorGemm`
   * - SYMM
     - --
     - --
     - --
   * - SYRK
     - --
     - --
     - --
   * - SYRK2
     - --
     - --
     - --
   * - HEMM
     - --
     - --
     - --
   * - HERK
     - --
     - --
     - --
   * - HERK2
     - --
     - --
     - --
   * - TRMM
     - `SerialTrmm`
     - `TeamTrmm`
     - `TeamVectorTrmm`
   * - TRSM
     - `SerialTrsm`
     - `TeamTrsm`
     - `TeamVectorTrsm`

LAPACK support
==============

Below are tables summarizing the currently supported function calls in Kokkos Kernels.

.. list-table::
   :widths: 12 10 10 10
   :header-rows: 1

   * - LAPACK Call
     - Serial
     - Team
     - TeamVector
   * - getrf
     - :doc:`SerialGetrf <dense/batched_getrf>`
     - --
     - --
   * - getrs
     - :doc:`SerialGetrs <dense/batched_getrs>`
     - --
     - --
   * - gesv
     - `SerialGesv`
     - `TeamGesv`
     - `TeamVectorGesv`
   * - gbtrf
     - :doc:`SerialGbtrf <dense/batched_gbtrf>`
     - --
     - --
   * - gbtrs
     - :doc:`SerialGbtrs <dense/batched_gbtrs>`
     - --
     - --
   * - pbtrf
     - :doc:`SerialPbtrf <dense/batched_pbtrf>`
     - --
     - --
   * - pbtrs
     - :doc:`SerialPbtrs <dense/batched_pbtrs>`
     - --
     - --
   * - pttrf
     - :doc:`SerialPttrf <dense/batched_pttrf>`
     - --
     - --
   * - pttrs
     - :doc:`SerialPttrs <dense/batched_pttrs>`
     - --
     - --
   * - getri
     - --
     - --
     - --
   * - trtri
     - `SerialTrtri`
     - --
     - --
   * - geequ
     - --
     - --
     - --
   * - geqrf
     - `SerialQR`
     - `TeamQR`
     - `TeamVectorQR`
   * - ungqr
     - --
     - --
     - --
   * - unmqr
     - --
     - --
     - --
   * - gesvd
     - `SerialSVD`
     - --
     - --
   * - geev
     - --
     - --
     - --
   * - lacgv
     - :doc:`SerialLacgv <dense/batched_lacgv>`
     - --
     - --
   * - laswp
     - :doc:`SerialLaswp <dense/batched_laswp>`
     - --
     - --
