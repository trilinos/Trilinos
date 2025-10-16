KokkosSparse::CooMatrix<>::CooMatrix
####################################

.. code:: cppkokkos

  KOKKOS_INLINE_FUNCTION CooMatrix();

  CooMatrix(size_type nrows, size_type ncols, row_view row_in, column_view col_in, scalar_view data_in);

Constructs a CooMatrix from specified inputs.

1. Default constructor: creates a 0x0 matrix with empty views for rows, columns, and values.
2. Constructor which shallow-copies the provided views.

Parameters
==========

:nrows: The number of rows in the matrix.

:ncols: The number of columns in the matrix.

:row_in: The row indices of the entries.

:col_in: The column indices of the entries.

:data_in: The numerical values of the entries.

Type Requirements
-----------------

- `OrdinalType` must be an integer type.
