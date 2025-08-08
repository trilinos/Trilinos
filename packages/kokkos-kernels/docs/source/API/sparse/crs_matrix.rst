KokkosSparse::CrsMatrix
#######################

.. toctree::
   :maxdepth: 1
   :hidden:

   CrsMatrix_constructors
   sparse_row_view

Defined in header ``KokkosSparse_CrsMatrix.hpp``

.. code:: cppkokkos

  template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
            class SizeType = KokkosKernels::default_size_type>
  class CrsMatrix;

``KokkosSparse::CrsMatrix`` provides a compressed sparse row implementation of a sparse matrix, as described, for example, in Saad (2nd ed.).

Template Parameters
===================

:ScalarType: type of the values stored in the matrix.

:OrdinalType: type of the indices storied in the matrix (column indices).

:Device: type of the ``Kokkos::Device`` of the matrix

:MemoryTraits: type of Kokkos memory traits used for the views of the matrix.

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

   * - const_type
     - Type of an indentical matrix storing const values.


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

Member Functions
================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - `constructor <CrsMatrix_constructors.html>`_
     - Construct a CrsMatrix from inputs, this can perform shallow or deep copies of input parameters depending on inputs and semantic.

   * - :ref:`operator= <crsmatrix_operator=>`
     - Assignment operator, will attempt to assign the data (graph and values) of the right hand side to the left hand side. If the assignment operator of graph or values fail, matrix assignment fails.

   * - :ref:`numRows <crsmatrix_numRows>`
     - Returns the number of rows in the matrix.

   * - :ref:`numCols <crsmatrix_numCols>`
     - Returns the number of columns in the matrix.

   * - :ref:`nnz <crsmatrix_nnz>`
     - Returns the number of structural non-zero values in the matrix (some of these might actually store zero).

   * - :ref:`row <crsmatrix_row>`
     - Returns a SparseRowView object from a row of the matrix.

   * - :ref:`rowConst <crsmatrix_rowConst>`
     - Returns a SparseRowViewConst object from a row of the matrix.

.. _crsmatrix_operator=:

operator=
^^^^^^^^^

.. code:: cppkokkos

  template <typename aScalarType, typename aOrdinalType, class aDevice, class aMemoryTraits, typename aSizeType>
  CrsMatrix& operator=(const CrsMatrix<aScalarType, aOrdinalType, aDevice, aMemoryTraits, aSizeType>& mtx);

Attempts to assign the underlying ``graph`` and ``values`` of the input matrix ``mtx`` to the matrix.

.. _crsmatrix_numRows:

numRows
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const;

Returns the number of rows in the matrix.

.. _crsmatrix_numCols:

numCols
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const;

Returns the number of columns in the matrix.

.. _crsmatrix_nnz:

nnz
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type nnz() const;

Returns the number of non-zero entries in the matrix.

.. _crsmatrix_row:

row
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  SparseRowView<CrsMatrix> row(const ordinal_type i) const;

Returns a view of the i-th row of the matrix as a :doc:`SparseRowView <sparse_row_view>`.

.. _crsmatrix_rowConst:

rowConst
^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst<CrsMatrix> row(const ordinal_type i) const;

Returns a view of the i-th row of the matrix as a :doc:`SparseRowViewConst <sparse_row_view>`.



Example
=======

.. literalinclude:: ../../../../example/wiki/sparse/KokkosSparse_wiki_crsmatrix.cpp
  :language: c++
