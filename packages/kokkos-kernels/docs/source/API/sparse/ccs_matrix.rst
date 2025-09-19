KokkosSparse::CcsMatrix
#######################

.. toctree::
   :maxdepth: 1
   :hidden:

   CcsMatrix_constructors

Defined in header ``KokkosSparse_CcsMatrix.hpp``

.. code:: cppkokkos

  template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
            class SizeType = typename Kokkos::ViewTraits<OrdinalType*, Device, void, void>::size_type>
  class CcsMatrix;

``KokkosSparse::CcsMatrix`` provides a compressed sparse column implementation of a sparse matrix.

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

   * - Member type
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

   * - staticccsgraph_type
     - Type of the underlying graph.

   * - index_type
     - Type of the view storing the column indices in the underlying graph.

   * - col_map_type
     - Type of the view storing the column offsets in the underlying graph.

   * - values_type
     - Type of the view storing the values stored in the matrix.


Data Members
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Data Member
     - Definition

   * - staticccsgraph_type graph
     - The underlying static ccs graph of the matrix storing its sparsity pattern.

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

   * - `constructor <CcsMatrix_constructors.html>`_
     - Construct a CcsMatrix from inputs, this can perform shallow or deep copies of input parameters depending on inputs and semantic.

   * - :ref:`numRows <ccsmatrix_numRows>`
     - Returns the number of rows in the matrix.

   * - :ref:`numCols <ccsmatrix_numCols>`
     - Returns the number of columns in the matrix.

   * - :ref:`setNumRows <ccsmatrix_setNumRows>`
     - Modify the number of columns in the matrix.

   * - :ref:`numPointRows <ccsmatrix_numPointRows>`
     - Equivalent to `numRows()`.

   * - :ref:`numPointCols <ccsmatrix_numPointCols>`
     - Equivalent to `numCols()`.

   * - :ref:`nnz <ccsmatrix_nnz>`
     - Returns the number of structural non-zero values in the matrix (some of these might actually store zero).

.. _ccsmatrix_numRows:

numRows
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const;

Returns the number of rows in the matrix.

.. _ccsmatrix_numCols:

numCols
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const;

Returns the number of columns in the matrix.

.. _ccsmatrix_setNumRows:

setNumRows
^^^^^^^^^^

.. code:: cppkokkos

  void setNumRows(ordinal_type r);

Modify the number of rows in the sparse matrix.
This invalidates any algorithm handles which previously used this matrix.

.. _ccsmatrix_numPointRows:

numPointRows
^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numPointRows() const;

Equivalent to `numRows()`, since this is not a block matrix type.

.. _ccsmatrix_numPointCols:

numPointCols
^^^^^^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numPointCols() const;

Equivalent to `numCols()`, since this is not a block matrix type.

.. _ccsmatrix_nnz:

nnz
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type nnz() const;

Returns the number of non-zero entries in the matrix.

