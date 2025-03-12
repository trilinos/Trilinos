KokkosSparse::SparseRowView
###########################

Defined in header ``KokkosSparse_CrsMatrix.hpp``.

.. code:: cppkokkos

  template <class MatrixType>
  struct SparseRowView;

View of a row of a sparse matrix. This class provides a generic view of a row of a sparse matrix. We intended this class to view a row of a CrsMatrix, but MatrixType need not necessarily be CrsMatrix.

Template Parameters
===================

:MatrixType: type of the matrix that stores the row to be viewed, typically ``KokkosSparse::CrsMatrix``.

Member Types
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Member type
     - Definition

   * - value_type
     - The type of the values in the row.

   * - ordinal_type
     - The type of the column indices in the row.

Data Members
============

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Data Member
     - Definition

   * - const ordinal_type length;
     - Number of entries in the row.

Member Functions
================

.. list-table::
   :widths: 30 70
   :header-rows: 1
   :align: left

   * - Name
     - Definition

   * - :ref:`constructor <sparserowview_constructor>`
     - constructs a ``SparseRowView``.

   * - :ref:`value <sparserowview_value>`
     - (Const) reference to the value of entry i in this row of the sparse matrix.

   * - :ref:`colidx <sparserowview_colidx>`
     - (Const) reference to the column index of entry i in this row of the sparse matrix.

.. _sparserowview_constructor:

constructor
-----------

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst(value_type* const values, ordinal_type* const colidx__,
                     const ordinal_type& stride, const ordinal_type& count);

  template <class OffsetType>
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst(const typename MatrixType::values_type& values,
                     const typename MatrixType::index_type& colidx__,
                     const ordinal_type& stride, const ordinal_type& count,
		     const OffsetType& idx,
                     const typename std::enable_if<std::is_integral_v<OffsetType>, int>::type& = 0);

Constructs a ``SparseRowView`` instance from raw pointers or Kokkos View extracted from a sparse matrix.

1. Pointer based constructor.
2. View based constructor, constructor with offset into ``colidx`` array.

Type Requirements
^^^^^^^^^^^^^^^^^

- `OffsetType` must be an integral type (``std::is_integral_v<OffsetType> == true``).

Parameters
^^^^^^^^^^

:values: Array of the row's values.

:colidx\_\_: Array of the row's column indices.

:stride: (Constant) stride between matrix entries in each of the above arrays.

:count: Number of entries in the row.

:idx: Start offset into ``colidx`` array.

.. _sparserowview_value:

value
-----

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  value_type& value(const ordinal_type& i) const;

(Const) reference to the value of entry i in this row of the sparse matrix.

Parameters
^^^^^^^^^^

:i: Index of the entry to be accessed. Note that this index does not necessarily correspond to the column index or the local row index of the underlying sparse matrix.

.. _sparserowview_colidx:

colidx
------

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx(const ordinal_type& i) const;

(Const) reference to the column index of entry i in this row of the sparse matrix.

Parameters
^^^^^^^^^^

:i: Index of the entry to be accessed. Note that this index does not necessarily correspond to the column index or the local row index of the underlying sparse matrix.
