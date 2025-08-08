KokkosSparse::forward_sweep_gauss_seidel_apply
##############################################

Defined in header ``KokkosSparse_gauss_seidel.hpp``

.. code:: cppkokkos

  template <class ExecutionSpace, KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS,
            class KernelHandle, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_,
            typename x_scalar_view_t, typename y_scalar_view_t>
  void forward_sweep_gauss_seidel_apply(const ExecutionSpace &space, KernelHandle *handle,
                                        typename KernelHandle::const_nnz_lno_t num_rows,
                                        typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                        lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                        x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                        bool init_zero_x_vector, bool update_y_vector,
                                        typename KernelHandle::nnz_scalar_t omega, int numIter);

  template <KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS, class KernelHandle,
            typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_, typename x_scalar_view_t,
            typename y_scalar_view_t>
  void forward_sweep_gauss_seidel_apply(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t num_rows,
                                        typename KernelHandle::const_nnz_lno_t num_cols, lno_row_view_t_ row_map,
                                        lno_nnz_view_t_ entries, scalar_nnz_view_t_ values,
                                        x_scalar_view_t x_lhs_output_vec, y_scalar_view_t y_rhs_input_vec,
                                        bool init_zero_x_vector, bool update_y_vector,
                                        typename KernelHandle::nnz_scalar_t omega, int numIter);

Apply the Gauss-Seidel preconditioner to a linear system of equations :math:`Ax=b`.

1. Perform symbolic operations to setup the Gauss-Seidel/SOR preconditioner.
2. Same as 1. but uses the resources of the default execution space instance associated with ``KernelHandle::HandleExecSpace``.

Parameters
==========

:space: execution space instance describing the resources available to run the kernels.

:handle: an instance of ``KokkosKernels::KokkosKernelsHandle`` from which a gs_handle will be used to extract control parameters.

:num_rows, num_cols: the number of rows and columns of the input matrix. 

:rowmap: row map of the input matrix.

:entries: column indices of the input matrix.

:values: values of the input matrix.

:x_lhs_output_vec, y_rhs_input_vec: left hand side and right hand side vectors of the linear system.

:init_zero_x_vector: whether ``x_lhs_output_vec`` should be zero initialized.

:update_y_vector: whether ``y_rhs_input_vec`` has changed since the last call to one of the Gauss-Seidel apply functions using ``handle``.

:omega: damping parameter for successive over-relaxations.

:numIter: the number of preconditioner iteration to perform on the linear system.

Type Requirements
=================

- The types of the input parameters must be consistent with the types defined by the ``KernelHandle``:

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename lno_row_view_t_::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename lno_nnz_view_t_::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::const_nnz_scalar_t, typename scalar_nnz_view_t_::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::const_nnz_scalar_t, typename y_scalar_view_t::const_value_type>``
  - ``std::is_same_v<typename KernelHandle::nnz_scalar_t, typename x_scalar_view_t::value_type>``

- The views describing the matrix, ``rowmap``, ``entries`` and ``values``, should not have a ``Kokkos::LayoutStride``

  - ``!std::is_same_v<typename lno_row_view_t_::array_layout, Kokkos::LayoutStride>``
  - ``!std::is_same_v<typename lno_nnz_view_t_::array_layout, Kokkos::LayoutStride>``
  - ``!std::is_same_v<typename scalar_nnz_view_t_::array_layout, Kokkos::LayoutStride>``

.. note::

   The non-strided requirement in this function should probably be checked on earlier during the symbolic and numeric phases? Also are we really happy with something `not layout stride`? Should we check on `is_contiguous` instead?

Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_gauss_seidel.cpp
  :language: c++
  :lines: 16-

