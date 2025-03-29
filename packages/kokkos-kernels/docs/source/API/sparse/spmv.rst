KokkosSparse::spmv
##################

Defined in header ``KokkosSparse_spmv.hpp``

.. code:: cppkokkos

  template <class ExecutionSpace, class Handle, class AlphaType, class AMatrix,
            class XVector, class BetaType, class YVector>
  void spmv(const ExecutionSpace& space, Handle* handle, const char mode[],
            const AlphaType& alpha, const AMatrix& A, const XVector& x,
            const BetaType& beta, const YVector& y);

  template <class ExecutionSpace, class AlphaType, class AMatrix, class XVector,
            class BetaType, class YVector,
	    typename = std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>>>
  void spmv(const ExecutionSpace& space, const char mode[], const AlphaType& alpha,
            const AMatrix& A, const XVector& x,
            const BetaType& beta, const YVector& y);

  template <class Handle, class AlphaType, class AMatrix, class XVector,
            class BetaType, class YVector,
            typename = std::enable_if_t<!Kokkos::is_execution_space<Handle>::value>>
  void spmv(Handle* handle, const char mode[], const AlphaType& alpha,
            const AMatrix& A, const XVector& x,
            const BetaType& beta, const YVector& y);

  template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
  void spmv(const char mode[], const AlphaType& alpha, const AMatrix& A,
            const XVector& x, const BetaType& beta, const YVector& y);

Kokkos sparse matrix-vector multiply. Computes y := alpha*Op(A)*x + beta*y, where Op(A) is controlled by mode (see below).

.. math::

   y = \beta*y + \alpha*A*x


1. Scale the ``y`` vector by ``beta``, and accumulate the result of the sparse matrix-vector product (``A*x``), scaled by ``alpha``, into ``y``.
2. Calls 1. using an spmv handle with ``SPMV_FAST_SETUP`` algorithm for the handle parameter.
3. Calls 1. using an instance of ``Handle::ExecutionSpaceType`` for the execution space parameter.
4. Calls 1. using an spmv handle with ``SPMV_FAST_SETUP`` algorithm and an instance of ``Handle::ExecutionSpaceType`` for the handle and execution space parameters respectively.

Parameters
==========

:space: execution space instance.

:handle: an spmv handle that stores multiple parameters for algorithm and third party libraries choices at run time.

:mode: mode to be applied to ``A``, possible values are "N" (normal), "T" (transpose), "C" (conjugate) and "H" (hermitian or conjugate-transpose).

:alpha, beta: scaling coefficents for the matrix ``A`` and left hand side ``y`` vector respectively.

:A: The matrix used to perform the matrix-vector product.

:x, y: The right and left hand side vectors used as input and output respectively.

Type Requirements
-----------------

- ``AMatrix`` must be either a :doc:`KokkosSparse::CrsMatrix <crs_matrix>` or a :doc:`KokkosSparse::BsrMatrix <bsr_matrix>` and have a memory space compatible with the ``ExecutionSpace`` type:

  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename AMatrix::memory_space>::accessible == true``

- `ExecutionSpace` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- ``XVector`` and ``YVector`` are two Kokkos Views of same rank, either rank 1 or rank 2 with memory spaces accessible from ``ExecutionSpace`` and ``YVector`` must store non-const data:

  - ``Kokkos::is_view_v<XVector> == true && Kokkos::is_view_v<YVector> == true``
  - ``XVector::rank() == YVector::rank() && (XVector::rank() == 1 || XVector::rank() == 2)``
  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename XVector::memory_space>::accessible``
  - ``Kokkos::SpaceAccessibility<ExecutionSpace, typename YVector::memory_space>::accessible``
  - ``!std::is_const_v<typename YVector::value_type> == true``

Example
=======

.. code:: cppkokkos

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

    {
      const Scalar SC_ONE = Kokkos::ArithTraits<Scalar>::one();

      Ordinal numRows = 10;

      // Build the row pointers and store numNNZ                                                                                                                                          
      typename row_map_type::non_const_type row_map("row pointers", numRows + 1);
      typename row_map_type::HostMirror row_map_h = Kokkos::create_mirror_view(row_map);
      for(Ordinal rowIdx = 1; rowIdx < numRows + 1; ++rowIdx) {
        if( (rowIdx == 1) || (rowIdx == numRows) ){
          row_map_h(rowIdx) = row_map_h(rowIdx - 1) + 2;
        } else {
          row_map_h(rowIdx) = row_map_h(rowIdx - 1) + 3;
        }
      }
      const Offset numNNZ = row_map_h(numRows);
      Kokkos::deep_copy(row_map, row_map_h);

      typename entries_type::non_const_type entries("column indices", numNNZ);
      typename entries_type::HostMirror entries_h = Kokkos::create_mirror_view(entries);
      typename values_type::non_const_type values("values", numNNZ);
      typename values_type::HostMirror values_h = Kokkos::create_mirror_view(values);
      for(Ordinal rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        if(rowIdx == 0) {
          entries_h(0) = rowIdx;
          entries_h(1) = rowIdx + 1;

          values_h(0)  = SC_ONE;
          values_h(1)  = -SC_ONE;
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

      graph_type myGraph(entries, row_map);
      matrix_type myMatrix("test matrix", numRows, values, myGraph);

      const Scalar alpha = SC_ONE;
      const Scalar beta  = SC_ONE;

      typename values_type::non_const_type x("lhs", numRows);
      typename values_type::non_const_type y("rhs", numRows);

      KokkosSparse::spmv("N", alpha, myMatrix, x, beta, y);
    }

    Kokkos::finalize();

    return 0;
  }

