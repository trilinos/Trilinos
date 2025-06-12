KokkosSparse::gauss_seidel_symbolic
###################################

Defined in header ``KokkosSparse_gauss_seidel.hpp``

.. code:: cppkokkos

  template <typename ExecutionSpace, typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void gauss_seidel_symbolic(const ExecutionSpace &space, KernelHandle *handle,
                             typename KernelHandle::const_nnz_lno_t num_rows,
                             typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                             lno_nnz_view_t_ entries, bool is_graph_symmetric = true);

  template <typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_>
  void gauss_seidel_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                             typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                             lno_nnz_view_t_ entries, bool is_graph_symmetric = true);

Performs the symbolic phase of the Gauss-Seidel preconditioner setup. The details of the setup depend on the algorithm selected by the handle.

1. Perform symbolic operations to setup the Gauss-Seidel/SOR preconditioner.
2. Same as 1. but uses the resources of the default execution space instance associated with ``KernelHandle::HandleExecSpace``.

Parameters
==========

:space: execution space instance describing the resources available to run the kernels.

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which a gs_handle will be used to extract control parameters.

:num_rows, num_cols: the number of rows and columns of the input matrix.

:rowmap: row map of the input matrix.

:entries: column indices of the input matrix.

:is_graph_symmetric: control parameters indicating if the graph of the input matrix is symmetric. This information may be used by the implementation to chose faster implementations.

Type Requirements
=================

- The types of the input parameters must be consistent with the types defined by the ``KernelHandle``:

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>``

Example
=======

.. code:: cppkokkos

  #include "Kokkos_Core.hpp"
  #include "KokkosKernels_Handle.hpp"
  #include "KokkosSparse_IOUtils.hpp"
  #include "KokkosSparse_spmv.hpp"
  #include "KokkosSparse_CrsMatrix.hpp"
  #include "KokkosSparse_gauss_seidel.hpp"
  #include "KokkosBlas1_nrm2.hpp"
  
  int main() {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using MemSpace  = typename ExecSpace::memory_space;
    using Handle    = KokkosKernels::KokkosKernelsHandle<int, int, double, ExecSpace, MemSpace, MemSpace>;
    using Matrix    = KokkosSparse::CrsMatrix<double, int, ExecSpace, void, int>;
    using Vector    = typename Matrix::values_type;
    constexpr int numRows = 10000;
    Kokkos::initialize();
    {
      // Generate a square, diagonally dominant, but nonsymmetric matrix
      // on which Gauss-Seidel should converge in a few iterations.
      // Insert about 20 entries per row.
      int nnz = numRows * 20;
      Matrix A = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<Matrix>(
          numRows, numRows, nnz, 2, 100, 1.05);
      Handle handle;
      handle.create_gs_handle(KokkosSparse::GS_DEFAULT);
      // Symbolic setup (for A's sparsity pattern).
      KokkosSparse::gauss_seidel_symbolic(
          &handle, numRows, numRows, A.graph.row_map, A.graph.entries,
          /* whether matrix is structurally symmetric */ false);
      // Numeric setup (for A's values).
      // If A's values change but sparsity pattern remains the same,
      // gauss_seidel_numeric can be called again to reuse the handle.
      KokkosSparse::gauss_seidel_numeric(
          &handle, numRows, numRows, A.graph.row_map, A.graph.entries, A.values,
          /* whether matrix is structurally symmetric */ false);
      // Create random right-hand side vector b
      Vector b(Kokkos::view_alloc(Kokkos::WithoutInitializing, "b"), numRows);
      auto bHost = Kokkos::create_mirror_view(b);
      for (int i = 0; i < numRows; i++) bHost(i) = 3 * ((1.0 * rand()) / RAND_MAX);
      Kokkos::deep_copy(b, bHost);
      // Create uninitialized left-hand side (solution) vector
      Vector x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), numRows);
      // Compute initial residual norm, (for initial guess x = 0)
      double initialRes    = KokkosBlas::nrm2(b);
      double scaledResNorm = 1.0;
      bool firstIter    = true;
      Vector res(Kokkos::view_alloc(Kokkos::WithoutInitializing, "res"), numRows);
      // Iterate until reaching the tolerance
      int numIters = 0;
      while (scaledResNorm > 1e-6) {
        KokkosSparse::forward_sweep_gauss_seidel_apply(
            &handle, numRows, numRows, A.graph.row_map, A.graph.entries, A.values, x, b,
            /* whether to zero out x */ firstIter,
            /* that we are running with a new right-hand side b */ firstIter,
            /* damping factor (omega) */ 1.0,
            /* number of iterations */ 1);
        firstIter = false;
        // Compute residual: res := Ax - b
        Kokkos::deep_copy(res, b);
        KokkosSparse::spmv("N", 1.0, A, x, -1.0, res);
        // Recompute the scaled norm
        scaledResNorm = KokkosBlas::nrm2(res) / initialRes;
        numIters++;
        std::cout << "Iteration " << numIters << " scaled residual norm: " << scaledResNorm << '\n';
      }
      std::cout << "SUCCESS: converged in " << numIters << " iterations.\n";
    }
    Kokkos::finalize();
    return 0;
  }

