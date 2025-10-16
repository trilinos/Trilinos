KokkosSparse::BsrMatrix
#######################

.. toctree::
   :maxdepth: 1
   :hidden:

   BsrMatrix_constructors
   bsr_row_view

Defined in header ``KokkosSparse_BsrMatrix.hpp``

.. code:: cppkokkos

  template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
            class SizeType = KokkosKernels::default_size_type>
  class BsrMatrix;

``KokkosSparse::BsrMatrix`` provides a block compressed sparse row implementation of a sparse matrix, as described, for example, in Saad (2nd ed.).

Template Parameters
===================

:ScalarType: type of the values stored in the matrix.

:OrdinalType: type of the indices storied in the matrix (column indices).

:Device: type of the ``Kokkos::Device`` of the matrix

:MemoryTraits: type of memory traits used for the views of the matrix.

:SizeType: type used to store row offsets.


Member Types
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Memeber type
     - Definition

   * - execution_space
     - Alias for Device::execution_space.

   * - memory_space
     - Alias for Device::memory_space.

   * - device_type
     - Device associated with the matrix.

   * - value_type
     - Type of the values stored in the matrix, alias for ScalarType.

   * - ordinal_type
     - Type of the indices stored in the matrix, alias for OrdinalType.

   * - size_type
     - Type of the offsets stored in the matrix, alias for SizeType.

   * - memory_traits
     - Alias for MemoryTraits.

   * - HostMirror
     - CrsMatrix type templated on ScalarType, OrdinalType, host_mirror_space, MemoryTraits and SizeType.

   * - StaticCrsGraph
     - Type of the underlying static crs graph.

   * - staticcsrgraph
     - Alias for StaticCrsGraph.

   * - index_type
     - Type of the view storing the column indices in the underlying static crs graph.

   * - const_ordinal_type
     - Const variant for the type of the column indices stored in index_type.

   * - non_const_ordinal_type
     - Non-const variant for the type of the column indices stored in index_type.

   * - row_map_type
     - Type of the view storing the row offsets in the underlying static crs graph.

   * - const_size_type
     - Const variant for the type of the offsets stored in row_map_type.

   * - non_const_size_type
     - Non-const variant for the type of the offsets stored in row_map_type.

   * - values_type
     - Type of the view storing the values stored in the matrix.

   * - const_value_type
     - Const variant for the type of the values stored in values_type.

   * - non_const_value_type
     - Non-const variant for the type of the values stored in values_type.

   * - block_layout_type
     - The Kokkos layout used when extracting blocks from the matrix. It is an alias for ``Kokkos::LayoutRight``.

   * - block_type
     - Type used to access of a view of a block as Kokkos::View or rank 2.

   * - const_block_type
     - Type used to access of a const view of a block as Kokkos::View or rank 2.


Data Members
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Data Member
     - Definition

   * - staticcrsgraph_type graph
     - The underlying static crs graph of the matrix storing its sparsity pattern.

   * - values_type values
     - The values stored in the matrix.

   * - const ordinal blockDim\_
     - The size of the blocks stored by the matrix.

Member Functions
================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - `constructor <BsrMatrix_constructors.html>`_
     - Construct a CrsMatrix from inputs, this can perform shallow or deep copies of input parameters depending on inputs and semantic.

   * - :ref:`operator=`
     - Assignment operator, will attempt to assign the data (graph and values) of the right hand side to the left hand side. If the assignment operator of graph or values fail, matrix assignment fails.

   * - :ref:`numRows`
     - Returns the number of rows in the matrix.

   * - :ref:`numCols`
     - Returns the number of columns in the matrix.

   * - :ref:`setNumCols`
     - Modify the number of columns in the matrix.

   * - :ref:`blockDim`
     - Returns the dimension of the blocks stored in the matrix.

   * - :ref:`numPointRows <bsrmatrix_numPointRows>`
     - Returns the number of point rows in the matrix.

   * - :ref:`numPointCols <bsrmatrix_numPointCols>`
     - Returns the number of point columns in the matrix.

   * - :ref:`nnz`
     - Returns the number of structural non-zero values in the matrix (some of these might actually store zero).

   * - :ref:`block_row`
     - Returns a SparseRowView object from a row of the matrix.

   * - :ref:`block_row_Const`
     - Returns a SparseRowViewConst object from a row of the matrix.

   * - :ref:`unmanaged_block`
     - Return a view of the i-th block in the matrix.

   * - :ref:`unmanaged_block_const`
     - Return a const view of the i-th block in the matrix.

.. _operator=:

operator=
^^^^^^^^^

.. code:: cppkokkos

  template <typename aScalarType, typename aOrdinalType, class aDevice, class aMemoryTraits, typename aSizeType>
  BsrMatrix& operator=(const BsrMatrix<aScalarType, aOrdinalType, aDevice, aMemoryTraits, aSizeType>& mtx);

Attempts to assign the underlying ``graph`` and ``values`` of the input matrix ``mtx`` to the matrix.

.. _numRows:

numRows
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const;

Returns the number of rows in the matrix.

.. _numCols:

numCols
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const;

Returns the number of columns in the matrix.

.. _setNumCols:

setNumCols
^^^^^^^^^^

.. code:: cppkokkos

  void setNumCols(ordinal_type c);

Modify the number of columns in the sparse matrix.
This invalidates any algorithm handles which previously used this matrix.

.. _blockDim:

blockDim
^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type blockDim() const;

Returns the dimension of the blocks stored by the matrix.

.. _bsrmatrix_numPointRows:

numPointRows
^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numPointRows() const;

Returns the number of point rows in the matrix. This is the number of (block)
rows times the block size. It is also the dimension of the matrix's range.

.. _bsrmatrix_numPointCols:

numPointCols
^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numPointCols() const;

Returns the number of point columns in the matrix. This is the number of (block)
columns times the block size. It is also the dimension of the matrix's domain.

.. _nnz:

nnz
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type nnz() const;

Returns the number of non-zero entries in the matrix.

.. _block_row:

block_row
^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  BsrRowView<BsrMatrix> block_row(const ordinal_type i) const;

Returns a view of the i-th block row of the matrix as a :doc:`BsrRowView <bsr_row_view>`.

.. _block_row_Const:

block_row_Const
^^^^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  BsrRowViewConst<BsrMatrix> block_row_Const(const ordinal_type i) const;

Returns a view of the i-th block row of the matrix as a :doc:`BsrRowViewConst <bsr_row_view>`.

.. _unmanaged_block:

unmanaged_block
^^^^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  block_type unmanaged_block(const size_type i) const;

Return a view of the i-th block in the matrix.

.. _unmanaged_block_const:

unmanaged_block_const
^^^^^^^^^^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  const_block_type unmanaged_block_const(const size_type i) const;

Return a const view of the i-th block in the matrix.


Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_bsrmatrix_2.cpp
  :language: c++
  :lines: 16-

