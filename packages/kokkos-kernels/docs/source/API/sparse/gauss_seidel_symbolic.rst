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

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_gauss_seidel.cpp
  :language: c++
  :lines: 16-

