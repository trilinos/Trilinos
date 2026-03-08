KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_with_rcb_sequential
############################################################################

Defined in header: :code:`KokkosSparse_Utils.hpp`

.. code:: c++

  template <typename crsMat_t, typename perm_view_type>
  void kk_extract_diagonal_blocks_crsmatrix_with_rcb_sequential(
      const crsMat_t &A, const perm_view_type &perm_rcb, const perm_view_type &reverse_perm_rcb,
      const std::vector<typename crsMat_t::non_const_ordinal_type> partition_sizes_rcb,
      std::vector<crsMat_t> &DiagBlk_v);

Extract diagonal blocks from a square CRS matrix using the pre-run RCB partition information. This function must be called after applying RCB to the coordinates associated with the rows/columns of the CRS matrix.

The number of diagonal blocks is the number of partitions for RCB.

The function will return a vector of diagonal blocks corresponding to the RCB partitions.

The function will throw a runtime exception if any of the following conditions are true:

- the sparse matrix ``A`` is not square
- the size of the vector of diagonal blocks ``DiagBlk_v`` is not set or is larger than the number of rows in ``A``
- the size of the vector of diagonal blocks ``DiagBlk_v`` is not a power of 2
- the size of the vector ``DiagBlk_v`` is not equal to the size of the vector ``partition_sizes_rcb``
- the size of the permutation array ``perm_rcb`` is not equal to the number of rows in ``A``
- the size of the reverse permutation array ``reverse_perm_rcb`` is not equal to the number of rows in ``A``

Parameters
==========

:A: the input square CRS matrix (it is expected that column indices are in ascending order).

:perm_rcb: 1-D array describing the mapping from the original order to RCB order.

:reverse_perm_rcb: 1-D array describing the mapping from the RCB order to original order.

:partition_sizes_rcb: the vector containing sizes of RCB partitions

:DiagBlk_v: the vector of the extracted CRS diagonal blocks.

Type Requirements
=================

- ``crsMat_t`` must be a Kokkos Kernels `CrsMatrix <https://kokkos.org/kokkos-kernels/docs/API/sparse/crs_matrix.html>`_

- ``perm_view_type`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1

Example
=======

Will be added later.

