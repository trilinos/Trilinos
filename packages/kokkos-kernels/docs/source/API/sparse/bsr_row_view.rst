KokkosSparse::BsrRowView
########################

Defined in header ``KokkosSparse_BsrMatrix.hpp``.

.. code:: cppkokkos

  template <class MatrixType>
  struct BsrRowView;

View of a row of a block sparse matrix. This class provides a generic view of a block row of a sparse matrix. We intended this class to view a row of a BsrMatrix, but MatrixType need not necessarily be BsrMatrix.

Template Parameters
===================

:MatrixType: type of the matrix that stores the row to be viewed.

Member Types
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Memeber type
     - Definition

   * - value_type
     - The type the values stored in a block.

   * - ordinal_type
     - The type of the column indices.

   * - block_values_type
     - The type of a block in the row (a rank 2 Kokkos View).


Data Members
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Data Member
     - Definition

   * - length
     - Number of entries (i.e. blocks) in the row.

Member Functions
================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - BsrRowView
     - constructor of the block row view, may require the type of the matrix who's row is being viewed.

   * - local_row_in_block
     - Return a pointer to the first value of a given row in a block.

   * - local_block_value
     - Return a pointer to a specific value in a block of the row.

   * - block
     - Gets a rank 2 view representing the data of the requested block.

   * - findRelBlockOffset
     - Finds the offset associated with a given block in the row.

.. _bsrrowview_constructor:

constructor
-----------

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  BsrRowView(value_type* const values, ordinal_type* const colidx,
             const ordinal_type& blockDim, const ordinal_type& count);

  template <class OffsetType>
  KOKKOS_INLINE_FUNCTION
  BsrRowView(const typename MatrixType::values_type& values,
             const typename MatrixType::index_type& colidx, const ordinal_type& blockDim,
	     const ordinal_type& count, const OffsetType& start,
	     const typename std::enable_if<std::is_integral_v<OffsetType>, int>::type& = 0);

Constructs a ``BsrRowBiew`` instance from raw pointers or from views extracted from a block matrix.

1. Pointer based constructor, it assumes that the pointers passed are correctly offset for the desired row.
2. View based constructor, will compute offsets internally.

.. note::

   I am not quite sure why we have the ``std::enable_if<>`` for the second constructor? Would it be possible to replace it with a ``static_assert()`` on the same condition?

Type Requirements
^^^^^^^^^^^^^^^^^

- `OffsetType` must be an integral type.

Parameters
^^^^^^^^^^

:values: Array of the row's values.

:colidx: Array of the row's column indices.

:blockDim: Stride between rows within a block in the above arrays.

:count: The number of blocks in the row.

:start: Offset into values and colidx of the desired block-row start.

.. _local_row_in_block:

local_row_in_block
------------------

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  value_type* local_row_in_block(const ordinal_type& K, const ordinal_type& i) const;

Return a pointer offset to local row i of block K of ``values_`` array; user responsible for indexing into this pointer correctly.

Parameters
^^^^^^^^^^

:K: must be the LOCAL block index within this block-row.

:i: must be the LOCAL row index offset within this block-row.

.. _local_block_value:

local_block_value
-----------------

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  value_type& local_block_value(const ordinal_type& K, const ordinal_type& i,
                                const ordinal_type& j) const;

Return the value at a specified block K of block-row with local row and col offset (i,j).

Parameters
^^^^^^^^^^

:K: must be the LOCAL block index within this block-row.

:i: must be the LOCAL row index offset within this block-row.

:j: must be the LOCAL col index offset within this block-row.

.. _block:

block
-----

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  block_values_type block(const ordinal_type& K) const;

Return unmanaged 2D strided View wrapping local block K from this block-row.

Parameters
^^^^^^^^^^

:K: must be the LOCAL block index within this block-row.

.. _findRelBlockOffset:

findRelBlockOffset
------------------

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  ordinal_type findRelBlockOffset(const ordinal_type idx_to_match, bool /*is_sorted*/ = false) const;

Return offset into ``colidx_`` for the requested block idx. If none found, return ``KokkosKernels::ArithTraits<ordinal_type>::max()``.

Parameters
^^^^^^^^^^

:id_to_math: The index to find in the row, if the index is not found ``KokkosKernels::ArithTraits<ordinal_type>::max()`` is returned.
