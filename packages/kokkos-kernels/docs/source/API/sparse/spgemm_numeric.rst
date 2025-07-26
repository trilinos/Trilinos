KokkosSparse::spgemm_numeric
############################

Defined in header ``KokkosSparse_spgemm.hpp``

.. note::
  The minimal header ``KokkosSparse_spgemm_numeric.hpp`` is also sufficient for version 1 below.

.. code:: cppkokkos

  template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename ascalar_nnz_view_t_,
            typename blno_row_view_t_, typename blno_nnz_view_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_,
            typename clno_nnz_view_t_, typename cscalar_nnz_view_t_>
  void spgemm_numeric(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                      typename KernelHandle::const_nnz_lno_t n, typename KernelHandle::const_nnz_lno_t k,
                      alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA, ascalar_nnz_view_t_ valuesA,
                      bool transposeA, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB, bscalar_nnz_view_t_ valuesB,
                      bool transposeB, clno_row_view_t_ row_mapC, clno_nnz_view_t_ &entriesC,
                      cscalar_nnz_view_t_ &valuesC,
                      typename KernelHandle::const_nnz_lno_t block_dim = 1);

  template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
  void spgemm_numeric(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode,
                      CMatrix& C);

Performs the numeric phase of the multiplication of two sparse matrices.

The symbolic phase (:doc:`spgemm_symbolic <spgemm_symbolic>`)
must have previously been called with the same handle and input matrices (except for scalar values).

.. math::

   C = A*B

1. Computes the column indices and values of matrix ``C``. Also computes the row map of ``C`` if it was not previously computed by :doc:`spgemm_symbolic <spgemm_symbolic>`.
2. Same as 1. but uses :doc:`CrsMatrix <crs_matrix>` objects for a more compact interface.

.. warning::
   While transpose flags are part of the interface, the algorithm is only implemented for non-transposed A and B.

``spgemm_numeric`` can reuse the analysis performed by :doc:`spgemm_symbolic <spgemm_symbolic>`,
as long as the sparsity patterns (row map and entries) of ``A`` and ``B`` do not change.
If only the values of ``A`` and ``B`` have changed, then ``spgemm_numeric`` can always be called
with the same handle to recompute the values of ``C``.

Parameters
==========

.. note::
   For all parameters that are shared with ``spgemm_symbolic`` (dimensions, row maps, entries and matrices), the parameters to ``spgemm_numeric`` must match what was previously passed into ``spgemm_symbolic``
   with the same handle. For example, it is not valid to call ``spgemm_symbolic`` with one version of ``entriesA``, and then call ``spgemm_numeric`` with a different ``entriesA`` but the same handle.

:handle: spadd kernels handle obtained from an instance of ``KokkosKernels::KokkosKernelsHandle``. Must be the same handle that was previously passed into ``spgemm_symbolic`` with the same input matrices.

:m, n, k: dimensions of the matrices. **m** is the number of rows in ``A`` and ``C``. **n** is the number of columns in ``A`` and the number of rows in ``B``. **k** is the number of columns in ``B`` and ``C``.

:row_mapA, row_mapB: row maps of the input matrices ``A`` and ``B``.

:entriesA, entriesB: column indices of the entries in each row of ``A`` and ``B``. In most cases, these should be in sorted and merged order (see note in :doc:`spgemm_symbolic <spgemm_symbolic>`).

:valuesA, valuesB: values of the entries in each row of ``A`` and ``B``.

:transposeA, transposeB: transposition operator to be applied to ``row_mapA``, ``entriesA`` and ``row_mapB``, ``entriesB`` respectively. **Must both be false**; the algorithm is not yet implemented for other cases.

:row_mapC: row map of the output matrix ``C``.

:entriesC: column indices of the entries in each row of ``C``. Must be allocated (though not initialized) by the user to size ``handle->get_spgemm_handle()->get_c_nnz()``.

:valuesC: values of the entries in each row of ``C``. Must be allocated (though not initialized) by the user to size ``handle->get_spgemm_handle()->get_c_nnz()``.

:Amode, Bmode: transposition operation to apply to ``A`` and ``B`` respectively. **Must both be false**; the algorithm is not yet implemented for other cases.

:A, B, C: three instances of :doc:`KokkosSparse::CrsMatrix <crs_matrix>`. ``A`` and ``B`` are inputs parameters only. The output is written into the values of ``C``. In most cases, ``A`` and ``B`` should be in sorted and merged order (see note in :doc:`spgemm_symbolic <spgemm_symbolic>`).

Type Requirements
-----------------

- `alno_row_view_t_`, `alno_nnz_view_t_` and `ascalar_nnz_view_t_` must be compatible with the types expected by `KernelHandle`

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename alno_row_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename alno_nnz_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_scalar_t, typename ascalar_nnz_view_t_::const_value_type> == true``

- `blno_row_view_t_`, `blno_nnz_view_t_` and `bscalar_nnz_view_t_` must be compatible with the types expected by `KernelHandle`

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename blno_row_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename blno_nnz_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_scalar_t, typename bscalar_nnz_view_t_::const_value_type> == true``

- `clno_row_view_t_`, `clno_nnz_view_t_` and `cscalar_nnz_view_t_` must be compatible with the type expected by `KernelHandle` and `clno_nnz_view_t_` and `cscalar_nnz_view_t_` must be non-const.

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename clno_row_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename clno_nnz_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_scalar_t, typename cscalar_nnz_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename clno_nnz_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename cscalar_nnz_view_t_::value_type, typename cscalar_nnz_view_t_::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_spgemm.cpp
  :language: c++
  :lines: 16-
