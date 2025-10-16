API: Sparse
###########

.. toctree::
   :maxdepth: 2
   :hidden:

   sparse/bsr_matrix
   sparse/crs_matrix
   sparse/ccs_matrix
   sparse/coo_matrix
   sparse/kokkoskernelshandle
   sparse/sort_crs

   sparse/spmv
   sparse/spadd_symbolic
   sparse/spadd_numeric
   sparse/spgemm_symbolic
   sparse/spgemm_numeric
   sparse/sptrsv_symbolic
   sparse/sptrsv_solve

   sparse/gauss_seidel_symbolic
   sparse/gauss_seidel_numeric
   sparse/gauss_seidel_apply_symmetric
   sparse/gauss_seidel_apply_forward
   sparse/gauss_seidel_apply_backward
   sparse/spiluk_symbolic
   sparse/spiluk_numeric
   sparse/par_ilut

   sparse/extract_diagonal_blocks_rcb.rst

Containers
==========

Sparse containers are data structures that store indices and values associated with linear algebra objects, they allow storage savings all values associated with indices not explicitly stored are assumed to be zeros.

- :doc:`CrsMatrix <sparse/crs_matrix>`
- :doc:`BsrMatrix <sparse/bsr_matrix>`
- :doc:`CcsMatrix <sparse/ccs_matrix>`
- :doc:`CooMatrix <sparse/coo_matrix>`
- :doc:`KokkosKernelsHandle <sparse/kokkoskernelshandle>`

Sorting
=========================

- Sort sparse data structures

  - :doc:`sort_crs_matrix <sparse/sort_crs>`
  - :doc:`sort_brs_matrix <sparse/sort_crs>`
  - :doc:`sort_and_merge_matrix <sparse/sort_crs>`
  - :doc:`sort_crs_graph <sparse/sort_crs>`
  - :doc:`sort_and_merge_graph <sparse/sort_crs>`

Linear Algebra Operations
=========================

- sparse matrix-vector multiplication
   
  - :doc:`spmv <sparse/spmv>`

- sparse matrix-matrix addition

  - :doc:`spadd_symbolic <sparse/spadd_symbolic>`
  - :doc:`spadd_numeric <sparse/spadd_numeric>`

- sparse matrix-matrix multiplication

  - :doc:`spgemm_symbolic <sparse/spgemm_symbolic>`
  - :doc:`spgemm_numeric <sparse/spgemm_numeric>`

- sparse triangular matrix solve

  - :doc:`sptrsv_symbolic <sparse/sptrsv_symbolic>`
  - :doc:`sptrsv_solve <sparse/sptrsv_solve>`

Linear Solver / Preconditioners
===============================

- Gauss-Seidel and SOR

  - :doc:`gauss_seidel_symbolic <sparse/gauss_seidel_symbolic>`
  - :doc:`gauss_seidel_numeric <sparse/gauss_seidel_numeric>`
  - :doc:`symmetric_gauss_seidel_apply <sparse/gauss_seidel_apply_symmetric>`
  - :doc:`forward_sweep_gauss_seidel_apply <sparse/gauss_seidel_apply_forward>`
  - :doc:`backward_sweep_gauss_seidel_apply <sparse/gauss_seidel_apply_backward>`

- incomplete LU factorization

  - :doc:`spiluk_symbolic <sparse/spiluk_symbolic>`
  - :doc:`spiluk_numeric <sparse/spiluk_numeric>`
  - :doc:`par_ilut <sparse/par_ilut>`

Utility Functions
=================

- Apply RCB to the coordinates associated with rows/columns of a crs matrix then extract the diagonal blocks corresponding to the RCB partitions

  - :doc:`extract_diagonal_blocks_rcb <sparse/extract_diagonal_blocks_rcb>`
