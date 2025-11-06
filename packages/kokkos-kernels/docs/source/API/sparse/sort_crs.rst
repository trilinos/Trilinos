KokkosSparse::SortCRS
#####################

Defined in header ``KokkosSparse_sort_crs.hpp``

SortAlgorithm Enum
==================

.. code:: cppkokkos

  enum class SortAlgorithm {
    DEFAULT,
    PARALLEL_THREAD_LEVEL,
    BULK_SORT
  };

The `SortAlgorithm` enum specifies the sorting strategy to use for CRS matrices and graphs. The available options are:

- `DEFAULT`: Automatically selects the best sorting strategy based on the execution space and matrix properties.
- `PARALLEL_THREAD_LEVEL`: Forces parallel thread-level sorting within each row.
- `BULK_SORT`: Order all entries in the matrix/graph using a single sort-by-key. This is the default algorithm for highly imbalanced matrices or graphs.

Functions
=========

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Function
     - Description

   * - :ref:`sort_crs_matrix <sort_crs_matrix>`
     - Sorts the adjacent column list for each row of a CRS matrix into ascending order. Permutes the values accordingly.

   * - :ref:`sort_bsr_matrix <sort_bsr_matrix>`
     - Sorts the rows of a Block Row Storage (BRS) matrix, permuting the values accordingly.

   * - :ref:`sort_crs_graph <sort_crs_graph>`
     - Sorts the adjacent column list for each row of a CRS graph into ascending order.

   * - :ref:`sort_and_merge_matrix <sort_and_merge_matrix>`
     - Produces a new CRS matrix that is sorted and has no duplicate entries. Values for duplicate entries are summed.

   * - :ref:`sort_and_merge_graph <sort_and_merge_graph>`
     - Produces a new CRS graph that is sorted and has no duplicate entries.

.. _sort_crs_matrix:

sort_crs_matrix
^^^^^^^^^^^^^^^

.. code:: cppkokkos

  template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t>
  void sort_crs_matrix(const execution_space& exec, const rowmap_t& rowmap, const entries_t& entries,
                       const values_t& values,
                       typename entries_t::non_const_value_type numCols =
                           KokkosKernels::ArithTraits<typename entries_t::non_const_value_type>::max(),
                       SortAlgorithm option = SortAlgorithm::DEFAULT);

Sorts the adjacent column list for each row of a CRS matrix into ascending order. Permutes the values accordingly.

Template Parameters
===================

:execution_space: The Kokkos execution space to use for parallel operations.
:rowmap_t: Type of the row map view.
:entries_t: Type of the entries view.
:values_t: Type of the values view.

Parameters
==========

:exec: Execution space instance.
:rowmap: Row map view of the CRS matrix.
:entries: Entries view of the CRS matrix.
:values: Values view of the CRS matrix.
:numCols: Number of columns in the matrix (optional).
:option: Sorting strategy (optional, defaults to `SortAlgorithm::DEFAULT`).

.. _sort_bsr_matrix:

sort_bsr_matrix
^^^^^^^^^^^^^^^

.. code:: cppkokkos

  template <typename execution_space, typename rowmap_t, typename entries_t, typename values_t, typename Ordinal>
  void sort_bsr_matrix(const execution_space& exec, Ordinal blockSize, const rowmap_t& rowmap, const entries_t& entries,
                       const values_t& values,
                       typename entries_t::non_const_value_type numCols =
                           KokkosKernels::ArithTraits<typename entries_t::non_const_value_type>::max());

Sorts the rows of a Block Row Storage (BRS) matrix, permuting the values accordingly.

Template Parameters
===================

:execution_space: The Kokkos execution space to use for parallel operations.
:rowmap_t: Type of the row map view.
:entries_t: Type of the entries view.
:values_t: Type of the values view.
:Ordinal: Type of the block size.

Parameters
==========

:exec: Execution space instance.
:blockSize: Size of the blocks in the BRS matrix.
:rowmap: Row map view of the BRS matrix.
:entries: Entries view of the BRS matrix.
:values: Values view of the BRS matrix.
:numCols: Number of columns in the matrix (optional).

.. _sort_crs_graph:

sort_crs_graph
^^^^^^^^^^^^^^

.. code:: cppkokkos

  template <typename execution_space, typename rowmap_t, typename entries_t>
  void sort_crs_graph(const execution_space& exec, const rowmap_t& rowmap, const entries_t& entries,
                      typename entries_t::non_const_value_type numCols =
                          KokkosKernels::ArithTraits<typename entries_t::non_const_value_type>::max(),
                      SortAlgorithm option = SortAlgorithm::DEFAULT);

Sorts the adjacent column list for each row of a CRS graph into ascending order.

Template Parameters
===================

:execution_space: The Kokkos execution space to use for parallel operations.
:rowmap_t: Type of the row map view.
:entries_t: Type of the entries view.

Parameters
==========

:exec: Execution space instance.
:rowmap: Row map view of the CRS graph.
:entries: Entries view of the CRS graph.
:numCols: Number of columns in the graph (optional).
:option: Sorting strategy (optional, defaults to `SortAlgorithm::DEFAULT`).

.. _sort_and_merge_matrix:

sort_and_merge_matrix
^^^^^^^^^^^^^^^^^^^^^

.. code:: cppkokkos

  template <typename exec_space, typename rowmap_t, typename entries_t, typename values_t>
  void sort_and_merge_matrix(const exec_space& exec, const typename rowmap_t::const_type& rowmap_in,
                             const entries_t& entries_in, const values_t& values_in, rowmap_t& rowmap_out,
                             entries_t& entries_out, values_t& values_out,
                             typename entries_t::const_value_type& numCols =
                                 KokkosKernels::ArithTraits<typename entries_t::non_const_value_type>::max(),
                             SortAlgorithm option = SortAlgorithm::DEFAULT);

Produces a new CRS matrix that is sorted and has no duplicate entries. Values for duplicate entries are summed.

Template Parameters
===================

:exec_space: The Kokkos execution space to use for parallel operations.
:rowmap_t: Type of the row map view.
:entries_t: Type of the entries view.
:values_t: Type of the values view.

Parameters
==========

:exec: Execution space instance.
:rowmap_in: Input row map view of the CRS matrix.
:entries_in: Input entries view of the CRS matrix.
:values_in: Input values view of the CRS matrix.
:rowmap_out: Output row map view of the CRS matrix.
:entries_out: Output entries view of the CRS matrix.
:values_out: Output values view of the CRS matrix.
:numCols: Number of columns in the matrix (optional).
:option: Sorting strategy (optional, defaults to `SortAlgorithm::DEFAULT`).

.. _sort_and_merge_graph:

sort_and_merge_graph
^^^^^^^^^^^^^^^^^^^^

.. code:: cppkokkos

  template <typename exec_space, typename rowmap_t, typename entries_t>
  void sort_and_merge_graph(const exec_space& exec, const typename rowmap_t::const_type& rowmap_in,
                            const entries_t& entries_in, rowmap_t& rowmap_out, entries_t& entries_out,
                            typename entries_t::const_value_type& numCols =
                                KokkosKernels::ArithTraits<typename entries_t::non_const_value_type>::max(),
                            SortAlgorithm option = SortAlgorithm::DEFAULT);

Produces a new CRS graph that is sorted and has no duplicate entries.

Template Parameters
===================

:exec_space: The Kokkos execution space to use for parallel operations.
:rowmap_t: Type of the row map view.
:entries_t: Type of the entries view.

Parameters
==========

:exec: Execution space instance.
:rowmap_in: Input row map view of the CRS graph.
:entries_in: Input entries view of the CRS graph.
:rowmap_out: Output row map view of the CRS graph.
:entries_out: Output entries view of the CRS graph.
:numCols: Number of columns in the graph (optional).
:option: Sorting strategy (optional, defaults to `SortAlgorithm::DEFAULT`).