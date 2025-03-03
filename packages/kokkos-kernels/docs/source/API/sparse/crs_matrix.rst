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

.. code:: cppkokkos

  #include <sstream>

  #include "Kokkos_Core.hpp"

  #include "KokkosKernels_default_types.hpp"
  #include "KokkosSparse_CrsMatrix.hpp"
  #include "KokkosSparse_spmv.hpp"

  using Scalar  = default_scalar;
  using Ordinal = default_lno_t;
  using Offset  = default_size_type;
  using Layout  = default_layout;

  int main(int argc, char* argv[]) {
    Kokkos::initialize();

    using device_type  = typename Kokkos::Device<Kokkos::DefaultExecutionSpace,
                                                 typename Kokkos::DefaultExecutionSpace::memory_space>;
    using matrix_type  = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void, Offset>;
    using graph_type   = typename matrix_type::staticcrsgraph_type;
    using row_map_type = typename graph_type::row_map_type;
    using entries_type = typename graph_type::entries_type;
    using values_type  = typename matrix_type::values_type;

    const Scalar SC_ONE = Kokkos::ArithTraits<Scalar>::one();

    Ordinal numRows = 10;

    {
      const Offset numNNZ = 2 + (numRows - 2)*3 + 2;
      typename row_map_type::non_const_type row_map("row pointers", numRows + 1);
      typename entries_type::non_const_type entries("column indices", numNNZ);
      typename values_type::non_const_type values("values", numNNZ);

      {
        // Build the row pointers and store numNNZ                                                                                                                                        
        typename row_map_type::HostMirror row_map_h = Kokkos::create_mirror_view(row_map);
        for(Ordinal rowIdx = 1; rowIdx < numRows + 1; ++rowIdx) {
          if( (rowIdx == 1) || (rowIdx == numRows) ){
            row_map_h(rowIdx) = row_map_h(rowIdx - 1) + 2;
          } else {
            row_map_h(rowIdx) = row_map_h(rowIdx - 1) + 3;
          }
        }
        Kokkos::deep_copy(row_map, row_map_h);
        if(row_map_h(numRows) != numNNZ) {
          std::ostringstream error_msg;
          error_msg << "error: row_map(numRows) != numNNZ, row_map_h(numRows)=" << row_map_h(numRows)
                    << ", numNNZ=" << numNNZ;
          throw std::runtime_error(error_msg.str());
        }

	      typename entries_type::HostMirror entries_h = Kokkos::create_mirror_view(entries);
        typename values_type::HostMirror values_h = Kokkos::create_mirror_view(values);
        for(Ordinal rowIdx = 0; rowIdx < numRows; ++rowIdx) {
          if(rowIdx == 0) {
            entries_h(row_map_h(rowIdx))     = rowIdx;
            entries_h(row_map_h(rowIdx) + 1) = rowIdx + 1;

            values_h(row_map_h(rowIdx))      = SC_ONE;
            values_h(row_map_h(rowIdx) + 1)  = -SC_ONE;
	        } else if(rowIdx == numRows - 1) {
            entries_h(row_map_h(rowIdx))     = rowIdx - 1;
            entries_h(row_map_h(rowIdx) + 1) = rowIdx;

            values_h(row_map_h(rowIdx))      = -SC_ONE;
            values_h(row_map_h(rowIdx) + 1)  = SC_ONE;
	        } else {
            entries_h(row_map_h(rowIdx))     = rowIdx - 1;
            entries_h(row_map_h(rowIdx) + 1) = rowIdx;
            entries_h(row_map_h(rowIdx) + 2) = rowIdx + 1;

            values_h(row_map_h(rowIdx))      = -SC_ONE;
            values_h(row_map_h(rowIdx) + 1)  = SC_ONE + SC_ONE;
            values_h(row_map_h(rowIdx) + 2)  = -SC_ONE;
          }
	      }
	      Kokkos::deep_copy(entries, entries_h);
        Kokkos::deep_copy(values, values_h);
      }

      graph_type myGraph(entries, row_map);
      matrix_type myMatrix("test matrix", numRows, values, myGraph);
    }
    Kokkos::finalize();
    return 0;
  }
