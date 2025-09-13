KokkosSparse::spgemm_symbolic
#############################

Defined in header ``KokkosSparse_spgemm.hpp``

.. note::
  The minimal header ``KokkosSparse_spgemm_symbolic.hpp`` is also sufficient for version 1 below.

.. code:: cppkokkos

  template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_,
            typename blno_row_view_t_, typename blno_nnz_view_t_, typename clno_row_view_t_>
  void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::const_nnz_lno_t m,
                       typename KernelHandle::const_nnz_lno_t n,
		       typename KernelHandle::const_nnz_lno_t k,
                       alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA, bool transposeA,
		       blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB, bool transposeB,
                       clno_row_view_t_ row_mapC, bool computeRowptrs = false);

  template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
  void spgemm_symbolic(KernelHandle& kh, const AMatrix& A, const bool Amode,
                       const BMatrix& B, const bool Bmode, CMatrix& C);

Performs the symbolic phase of the multiplication of two sparse matrices.

.. math::
   C = A * B

1. Computes the number of non-zeros in ``C`` and optionally the row map of ``C`` if ``computeRowptrs == true``.

.. note::
  Setting ``computeRowptrs == true`` may increase running time. Even if ``computeRowptrs == false``, ``row_mapC`` will always be populated after subsequently calling :doc:`spgemm_numeric <spgemm_numeric>`.

2. Compact interface using :doc:`CrsMatrix <crs_matrix>` objects. Computes the number of non-zeros in ``C``, and internally allocates the row map, entries and values of ``C`` without initializing them. Finally, constructs matrix ``C`` with these newly allocated views.

After ``spgemm_symbolic`` returns, the number of non-zeros in C is available through the spgemm subhandle (for example, ``handle->get_spgemm_handle()->get_c_nnz()``).

.. warning::
   While transpose flags are part of the interface, the algorithm is only implemented for non-transposed A and B.

.. note::
  Some underlying implementations require that both A and B are sorted and merged, meaning that the column indices within each row are sorted and unique.

  Expert users may want to use SpGEMM but avoid the cost of explicitly sorting and merging the input matrices.
  The implementations where this requirement doesn't apply are Kokkos Kernels, rocSPARSE (as long as there is no row with more than 8192 intermediate products), and MKL.

Parameters
==========

:handle: spadd kernels handle obtained from an instance of ``KokkosKernels::KokkosKernelsHandle``.

:m, n, k: dimensions of the matrices. **m** is the number of rows in ``A`` and ``C``. **n** is the number of columns in ``A`` and the number of rows in ``B``. **k** is the number of columns in ``B`` and ``C``.

:row_mapA, row_mapB: row maps of the input matrices ``A`` and ``B``.

:entriesA, entriesB: column indices of the entries in each row of ``A`` and ``B``. In most cases, these should be in sorted and merged order (see note above).

:transposeA, transposeB: transposition operator to be applied to ``row_mapA``, ``entriesA`` and ``row_mapB``, ``entriesB`` respectively. **Must both be false**; the algorithm is not yet implemented for other cases.

:row_mapC: row map of the output matrix ``C``. For version 1, it must already be allocated to size :math:`m + 1`. Initialization is not required.

:computeRowptrs: force the computation of ``row_mapC``. By default, only the number of non-zeros in ``C`` is guaranteed to be computed.

:Amode, Bmode: transposition operation to apply to ``A`` and ``B`` respectively. **Must both be false**; the algorithm is not yet implemented for other cases.

:A, B, C: three instances of :doc:`KokkosSparse::CrsMatrix <crs_matrix>`. ``A`` and ``B`` are input parameters. ``C`` is an output parameter only. In most cases, ``A`` and ``B`` should be sorted and merged (see note above).

Type Requirements
-----------------

- `alno_row_view_t_` and `alno_nnz_view_t_` must be compatible with the types expected by `KernelHandle`

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename alno_row_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename alno_nnz_view_t_::const_value_type> == true``

- `blno_row_view_t_` and `blno_nnz_view_t_` must be compatible with the types expected by `KernelHandle`

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename blno_row_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename KernelHandle::const_nnz_lno_t, typename blno_nnz_view_t_::const_value_type> == true``

- `clno_row_view_t_` must be non-const and compatible with the type expected by `KernelHandle`

  - ``std::is_same_v<typename KernelHandle::const_size_type, typename clno_row_view_t_::const_value_type> == true``
  - ``std::is_same_v<typename clno_row_view_t_::value_type, typename clno_row_view_t_::non_const_value_type> == true``

Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_spgemm.cpp
  :language: c++
  :lines: 16-

