KokkosSparse::CooMatrix
#######################

.. toctree::
   :maxdepth: 1
   :hidden:

   CooMatrix_constructors

Defined in header ``KokkosSparse_CooMatrix.hpp``

.. code:: cppkokkos

  template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
            class SizeType = typename Kokkos::ViewTraits<OrdinalType*, Device, void, void>::size_type>
  class CooMatrix;

``KokkosSparse::CooMatrix`` provides a COO (coordinate) implementation of a sparse matrix.

Template Parameters
===================

:ScalarType: type of the values stored in the matrix.

:OrdinalType: type of the indices storied in the matrix (row and column indices).

:Device: type of the ``Kokkos::Device`` of the matrix

:MemoryTraits: type of Kokkos memory traits used for the views of the matrix.

:SizeType: unused.


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

   * - ordinal_type
     - Type of the indices stored in the matrix, alias for OrdinalType.

   * - row_type 
     - Alias for ordinal_type.

   * - column_type
     - Alias for ordinal_type.

   * - size_type
     - Alias for SizeType.

   * - memory_traits
     - Alias for MemoryTraits.

   * - row_view
     - Type of the view storing the row indices.

   * - column_view
     - Type of the view storing the column indices.

   * - scalar_view
     - Type of the view storing the numerical values.

Member Functions
================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - `constructor <CooMatrix_constructors.html>`_
     - Construct a CooMatrix from inputs, this can perform shallow or deep copies of input parameters depending on inputs and semantic.

   * - :ref:`numRows <coomatrix_numRows>`
     - Returns the number of rows in the matrix.

   * - :ref:`numCols <coomatrix_numCols>`
     - Returns the number of columns in the matrix.

   * - :ref:`setNumRows <coomatrix_setNumRows>`
     - Modify the number of rows in the matrix.

   * - :ref:`setNumCols <coomatrix_setNumRows>`
     - Modify the number of columns in the matrix.

   * - :ref:`nnz <coomatrix_nnz>`
     - Returns the number of structural non-zero values in the matrix (some of these might actually store zero).

   * - :ref:`row <coomatrix_row>`
     - Returns the row indices view.

   * - :ref:`col <coomatrix_col>`
     - Returns the column indices view.

   * - :ref:`data <coomatrix_data>`
     - Returns the scalar values view.

.. _coomatrix_numRows:

numRows
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const;

Returns the number of rows in the matrix.

.. _coomatrix_numCols:

numCols
^^^^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const;

Returns the number of columns in the matrix.

.. _coomatrix_setNumRows:

setNumRows
^^^^^^^^^^

.. code:: cppkokkos

  void setNumRows(ordinal_type r);

Modify the number of rows in the sparse matrix.
This invalidates any algorithm handles which previously used this matrix.

.. _coomatrix_setNumCols:

setNumCols
^^^^^^^^^^

.. code:: cppkokkos

  void setNumCols(ordinal_type c);

Modify the number of columns in the sparse matrix.
This invalidates any algorithm handles which previously used this matrix.

.. _coomatrix_nnz:

nnz
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type nnz() const;

Returns the number of non-zero entries in the matrix.

.. _coomatrix_row:

row
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type row() const;

Returns the row indices view.

.. _coomatrix_col:

col
^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type col() const;

Returns the column indices view.

.. _coomatrix_data:

data
^^^^

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION size_type data() const;

Returns the scalar values view.

