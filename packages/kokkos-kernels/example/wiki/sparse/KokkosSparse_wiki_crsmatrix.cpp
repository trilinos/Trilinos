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

  using device_type = typename Kokkos::Device<
      Kokkos::DefaultExecutionSpace,
      typename Kokkos::DefaultExecutionSpace::memory_space>;
  using matrix_type =
      typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void,
                                       Offset>;
  using graph_type   = typename matrix_type::staticcrsgraph_type;
  using row_map_type = typename graph_type::row_map_type;
  using entries_type = typename graph_type::entries_type;
  using values_type  = typename matrix_type::values_type;

  const Scalar SC_ONE = Kokkos::ArithTraits<Scalar>::one();

  Ordinal numRows = 10;

  {
    const Offset numNNZ = 2 + (numRows - 2) * 3 + 2;
    typename row_map_type::non_const_type row_map("row pointers", numRows + 1);
    typename entries_type::non_const_type entries("column indices", numNNZ);
    typename values_type::non_const_type values("values", numNNZ);

    {
      // Build the row pointers and store numNNZ
      typename row_map_type::HostMirror row_map_h =
          Kokkos::create_mirror_view(row_map);
      for (Ordinal rowIdx = 1; rowIdx < numRows + 1; ++rowIdx) {
        if ((rowIdx == 1) || (rowIdx == numRows)) {
          row_map_h(rowIdx) = row_map_h(rowIdx - 1) + 2;
        } else {
          row_map_h(rowIdx) = row_map_h(rowIdx - 1) + 3;
        }
      }
      Kokkos::deep_copy(row_map, row_map_h);
      if (row_map_h(numRows) != numNNZ) {
        std::ostringstream error_msg;
        error_msg << "error: row_map(numRows) != numNNZ, row_map_h(numRows)="
                  << row_map_h(numRows) << ", numNNZ=" << numNNZ;
        throw std::runtime_error(error_msg.str());
      }

      typename entries_type::HostMirror entries_h =
          Kokkos::create_mirror_view(entries);
      typename values_type::HostMirror values_h =
          Kokkos::create_mirror_view(values);
      for (Ordinal rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        if (rowIdx == 0) {
          entries_h(row_map_h(rowIdx))     = rowIdx;
          entries_h(row_map_h(rowIdx) + 1) = rowIdx + 1;

          values_h(row_map_h(rowIdx))     = SC_ONE;
          values_h(row_map_h(rowIdx) + 1) = -SC_ONE;
        } else if (rowIdx == numRows - 1) {
          entries_h(row_map_h(rowIdx))     = rowIdx - 1;
          entries_h(row_map_h(rowIdx) + 1) = rowIdx;

          values_h(row_map_h(rowIdx))     = -SC_ONE;
          values_h(row_map_h(rowIdx) + 1) = SC_ONE;
        } else {
          entries_h(row_map_h(rowIdx))     = rowIdx - 1;
          entries_h(row_map_h(rowIdx) + 1) = rowIdx;
          entries_h(row_map_h(rowIdx) + 2) = rowIdx + 1;

          values_h(row_map_h(rowIdx))     = -SC_ONE;
          values_h(row_map_h(rowIdx) + 1) = SC_ONE + SC_ONE;
          values_h(row_map_h(rowIdx) + 2) = -SC_ONE;
        }
      }
      Kokkos::deep_copy(entries, entries_h);
      Kokkos::deep_copy(values, values_h);
    }

    graph_type myGraph(entries, row_map);
    matrix_type myMatrix("test matrix", numRows, values, myGraph);
    std::cout << "myMatrix has been created successfully:" << std::endl
              << "  - numRows=" << myMatrix.numRows() << std::endl
              << "  - numCols=" << myMatrix.numCols() << std::endl
              << "  - numNNZ= " << myMatrix.nnz() << std::endl;
  }

  Kokkos::finalize();

  return 0;
}
