KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_with_rcb_sequential
############################################################################

Defined in header: :code:`KokkosSparse_Utils.hpp`

.. code:: cppkokkos

  template <typename crsMat_t, typename coor_view_type, typename perm_view_type>
  void kk_extract_diagonal_blocks_crsmatrix_with_rcb_sequential(const crsMat_t &A, coor_view_type &coors,
                                                                std::vector<crsMat_t> &DiagBlk_v, perm_view_type &perm_rcb);

Extract diagonal blocks from a square CRS matrix after applying RCB to the coordinates associated with the rows/columns of the CRS matrix.

The number of diagonal blocks is used as the number of partitions for RCB.

The function will return

1. a vector of diagonal blocks corresponding to the RCB partitions, and
2. a permutation array describing the mapping from the original order to RCB order.

The function will throw a runtime exception if any of the following conditions are not met:

- the sparse matrix ``A`` is not square
- the size of the vector of diagonal blocks ``DiagBlk_v`` is not set or is larger than the number of rows in ``A``
- the size of the vector of diagonal blocks ``DiagBlk_v`` is not a power of 2

Parameters
==========

:A: the input square CRS matrix (it is expected that column indices are in ascending order).

:coors: 1/2/3-D coordinates associated with the rows of A.

:DiagBlk_v: the vector of the extracted CRS diagonal blocks.

:perm_rcb: 1-D array describing the mapping from the original order to RCB order.

Type Requirements
=================

- ``crsMat_t`` must be a Kokkos Kernels `CrsMatrix <https://kokkos.org/kokkos-kernels/docs/API/sparse/crs_matrix.html>`_

- ``coor_view_type`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 2

- ``perm_view_type`` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ of rank 1

Example
=======

Will be added later.

