KokkosSparse::spadd_symbolic
############################

Defined in header ``KokkosSparse_spadd.hpp``

.. code:: cppkokkos

  template <typename ExecSpace, typename KernelHandle, typename alno_row_view_t_,
            typename alno_nnz_view_t_, typename blno_row_view_t_,
            typename blno_nnz_view_t_, typename clno_row_view_t_>
  void spadd_symbolic(const ExecSpace &exec, KernelHandle *handle,
                      typename KernelHandle::const_nnz_lno_t m,  //same type as column indices
                      typename KernelHandle::const_nnz_lno_t n,
		      const alno_row_view_t_ a_rowmap,
                      const alno_nnz_view_t_ a_entries, const blno_row_view_t_ b_rowmap,
		      const blno_nnz_view_t_ b_entries, clno_row_view_t_ c_rowmap)

  template <typename KernelHandle, typename... Args>
  void spadd_symbolic(KernelHandle *handle, Args... args)

  template <typename ExecSpace, typename KernelHandle, typename AMatrix, typename BMatrix,
            typename CMatrix>
  void spadd_symbolic(const ExecSpace &exec, KernelHandle *handle, const AMatrix &A,
                      const BMatrix &B, CMatrix &C);

  template <typename KernelHandle, typename AMatrix, typename BMatrix, typename CMatrix>
  void spadd_symbolic(KernelHandle *handle, const AMatrix &A, const BMatrix &B, CMatrix &C);

Performs the symbolic phase of the addition of two sparse matrices.

.. math::

   C = \beta*C + \alpha*(A+B)


1. Computes the row map of ``C`` using the resources of ``exec`` and store the total number of non-zeros of ``C`` in the handle.
2. Same as 1. but use the resources of ``KernelHandle::HandleExecSpace{}``.
3. Compute the row map of ``C``, create and allocate views to hold column indices and values of ``C`` without initializing them using the resources of ``exec``. Construct matrix ``C`` from the row map and the newly created column indices and values views.
4. Same as 3. but use the resources of ``KernelHandle::HandleExecSpace{}``.

Parameters
==========

:exec: execution space instance.

:handle: spadd kernels handle obtained from an instance of ``KokkosKernels::KokkosKernelsHandle``.

:m, n: number of rows and column of the matrices ``A``, ``B`` and ``C``.

:a_rowmap, b_rowmap: row maps of the input matrices ``A`` and ``B``.

:a_entries, b_entries: column indices of the entries in each row of ``A`` and ``B``.

:c_rowmap: row map of the output matrix ``C``.

:A, B, C: three crs matrices.

Type Requirements
-----------------

Checks are performed in the :doc:`KokkosKernelsHandle <kokkoskernelshandle>` that was used to obtain the spadd handle.


Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_spadd.cpp
  :language: c++
  :lines: 16-
