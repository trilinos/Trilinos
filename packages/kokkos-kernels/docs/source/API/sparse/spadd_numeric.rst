KokkosSparse::spadd_numeric
###########################

Defined in header ``KokkosSparse_spadd.hpp``

.. code:: cppkokkos

  template <typename ExecSpace, typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
            typename ascalar_t_, typename ascalar_nnz_view_t_, typename blno_row_view_t_, typename blno_nnz_view_t_,
            typename bscalar_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_, typename clno_nnz_view_t_,
            typename cscalar_nnz_view_t_>
  void spadd_numeric(const ExecSpace &exec, KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                     typename KernelHandle::const_nnz_lno_t n, const alno_row_view_t_ a_rowmap,
                     const alno_nnz_view_t_ a_entries, const ascalar_nnz_view_t_ a_values, const ascalar_t_ alpha,
                     const blno_row_view_t_ b_rowmap, const blno_nnz_view_t_ b_entries,
                     const bscalar_nnz_view_t_ b_values, const bscalar_t_ beta, const clno_row_view_t_ c_rowmap,
                     clno_nnz_view_t_ c_entries, cscalar_nnz_view_t_ c_values);

  template <typename KernelHandle, typename... Args>
  void spadd_numeric(KernelHandle *handle, Args... args);

  template <typename ExecSpace, typename KernelHandle, typename AScalar, typename AMatrix, typename BScalar,
            typename BMatrix, typename CMatrix>
  void spadd_numeric(const ExecSpace &exec, KernelHandle *handle, const AScalar alpha, const AMatrix &A,
                     const BScalar beta, const BMatrix &B, CMatrix &C);

  template <typename KernelHandle, typename AScalar, typename AMatrix, typename BScalar, typename BMatrix,
            typename CMatrix>
  void spadd_numeric(KernelHandle *handle, const AScalar alpha, const AMatrix &A, const BScalar beta, const BMatrix &B,
                     CMatrix &C);

Performs the numeric phase of the addition of two sparse matrices.

.. math::

   C = \beta*C + \alpha*(A+B)

1. Assuming the symbolic phase has been called, i.e. the row map of ``C`` is correct and the entries and values views are allocated, compute the entries and values of ``C`` and return unique and sorted column indices and associated values using the resources associated with ``exec``.
2. Same as 1. but uses the resources of ``KernelHandle::HandleExecSpace{}``.
3. Same as 1. but uses matrix objects as inputs.
4. Same as 3. but uses the resources of ``KernelHandle::HandleExecSpace{}``.

Parameters
==========

:exec: execution space instance.

:handle: spadd kernels handle obtained from an instance of ``KokkosKernels::KokkosKernelsHandle``.

:m, n: number of rows and column of the matrices ``A``, ``B`` and ``C``.

:alpha, beta: scalar coefficients applied to scale the sum of ``A`` and ``B`` and the ``C`` matrix.

:a_rowmap, b_rowmap: row maps of the input matrices ``A`` and ``B``.

:a_entries, b_entries: column indices of the entries in each row of ``A`` and ``B``.

:a_values, b_values: values of the entries stored in ``A`` and ``B``.

:c_rowmap: row map of the output matrix ``C``.

:c_entries: column indices of the ``C`` matrix.

:c_values:  values of the entries in ``C``.

:A, B, C: three crs matrices.

Type Requirements
-----------------

Checks are performed in the :doc:`KokkosKernelsHandle <kokkoskernelshandle>` that was used to obtain the spadd handle.


Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_spadd.cpp
  :language: c++
  :lines: 16-
