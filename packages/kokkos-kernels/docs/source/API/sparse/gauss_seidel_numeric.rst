KokkosSparse::gauss_seidel_numeric
##################################

Defined in header ``KokkosSparse_gauss_seidel.hpp``

.. code:: cppkokkos

  template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
            typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(const ExecutionSpace &space, KernelHandle *handle,
                            typename KernelHandle::const_nnz_lno_t num_rows,
                            typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                            lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric = true);

  template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, typename KernelHandle,
            typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                            typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                            lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, bool is_graph_symmetric = true);

  template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
            typename KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(const ExecutionSpace &space, KernelHandle *handle,
                            typename KernelHandle::const_nnz_lno_t num_rows,
                            typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                            lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, scalar_nnz_view_t_ given_inverse_diagonal,
                            bool is_graph_symmetric = true);

  template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, typename KernelHandle,
            typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_>
  void gauss_seidel_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                            typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                            lno_nnz_view_t_ entries, scalar_nnz_view_t_ values, scalar_nnz_view_t_ given_inverse_diagonal,
                            bool is_graph_symmetric = true);

Performs the numeric phase of the Gauss-Seidel preconditioner setup. The details of the setup depend on the algorithm selected by the handle.

1. Perform symbolic operations to setup the Gauss-Seidel/SOR preconditioner.
2. Same as 1. but uses the resources of the default execution space instance associated with ``KernelHandle::HandleExecSpace``.
3. Same as 1. but takes the inverse of the diagonal values as input to accelerate the kernel.
4. Same as 3. but uses the resources of the default execution space instance associated with ``KernelHandle::HandleExecSpace``.

Parameters
==========

:space: execution space instance describing the resources available to run the kernels.

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which a gs_handle will be used to extract control parameters.

:num_rows, num_cols: the number of rows and columns of the input matrix. 

:rowmap: row map of the input matrix.

:entries: column indices of the input matrix.

:values: values of the input matrix.

:given_inverse_diagonal: (optional) user-provided vector of length ``num_rows``, containing the inverses of the diagonal elements of the input matrix. If not provided, the inverse diagonal will be computed internally.

:is_graph_symmetric: control parameters indicating if the graph of the input matrix is symmetric. This information may be used by the implementation to chose faster implementations.

Type Requirements
=================

- The types of the input parameters must be consistent with the types defined by the ``KernelHandle``:

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>``

Example
=======

This example shows how to use ``gauss_seidel_numeric``. It is called twice with the same handle:
once for a completely new matrix, and again after the matrix's values have changed.

See the examples for :doc:`symmetric_gauss_seidel_apply <gauss_seidel_apply_symmetric>`, :doc:`forward_sweep_gauss_seidel_apply <gauss_seidel_apply_forward>`,
and :doc:`backward_sweep_gauss_seidel_apply <gauss_seidel_apply_backward>` for how to apply the Gauss-Seidel preconditioner to a specific linear system.

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_gauss_seidel.cpp
  :language: c++
  :lines: 16-

