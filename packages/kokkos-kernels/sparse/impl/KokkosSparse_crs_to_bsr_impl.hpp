namespace KokkosSparse {

#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

namespace Impl {

/*! \brief Expand each entry of a crs matrix to a block in a bsr matrix
    The scalar, ordinal, and device types of the two matrices do not need
    to be compatible
*/
template <typename Bsr, typename Crs>
Bsr expand_crs_to_bsr(const Crs &crs, size_t blockSize) {
  using bsr_device_type     = typename Bsr::device_type;
  using bsr_execution_space = typename Bsr::execution_space;

  using crs_values_type = typename Crs::values_type;
  using bsr_values_type = typename Bsr::values_type;

  using crs_index_type = typename Crs::index_type;
  using bsr_index_type = typename Bsr::index_type;

  using crs_row_map_type = typename Crs::row_map_type;
  using bsr_row_map_type = Kokkos::View<typename Bsr::row_map_type::non_const_data_type,
                                        bsr_device_type>;  // need non-const version

  using bsr_size_type = typename Bsr::non_const_size_type;

  {
    size_t nnz = crs.nnz() * blockSize * blockSize;
    if (nnz > size_t(Kokkos::ArithTraits<bsr_size_type>::max())) {
      std::stringstream ss;
      ss << "expanding " << crs.nnz() << " non-zeros of CrsMatrix into blocks of " << blockSize
         << " would overflow size_type of requested BsrMatrix " << Kokkos::ArithTraits<bsr_size_type>::name();
      throw std::runtime_error(ss.str());
    }
  }

  // construct the Bsr row map
  bsr_row_map_type bsrRowMap("bsrRowMap", crs.graph.row_map.size());
  {
    // clone Crs row map in Bsr memory space
    Kokkos::View<typename crs_row_map_type::non_const_data_type, bsr_device_type> crows("crows",
                                                                                        crs.graph.row_map.size());
    Kokkos::deep_copy(crows, crs.graph.row_map);

    // copy to actual row map
    Kokkos::RangePolicy<bsr_execution_space> policy(0, crs.graph.row_map.size());
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(size_t i) { bsrRowMap(i) = crows(i); });
  }

  // construct the BSR col indices
  bsr_index_type bsrIndices("bsrIndices", crs.graph.entries.size());
  {
    // clone Crs row map in Bsr memory space
    Kokkos::View<typename crs_index_type::non_const_data_type, bsr_device_type> cinds("cinds",
                                                                                      crs.graph.entries.size());
    Kokkos::deep_copy(cinds, crs.graph.entries);

    // copy to actual row map
    Kokkos::RangePolicy<bsr_execution_space> policy(0, crs.graph.entries.size());
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(size_t i) { bsrIndices(i) = cinds(i); });
  }

  // construct BSR values
  bsr_values_type bsrVals("bsrVals", crs.nnz() * blockSize * blockSize);
  {
    // clone Crs values in Bsr memory space
    Kokkos::View<typename crs_values_type::non_const_data_type, bsr_device_type> cvals("cvals", crs.values.size());
    Kokkos::deep_copy(cvals, crs.values);

    // copy to actual values
    Kokkos::RangePolicy<bsr_execution_space> policy(0, crs.values.size());
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(size_t i) {
          for (size_t ii = i; ii < i + blockSize * blockSize; ++ii) {
            bsrVals(ii) = cvals(i);
          }
        });
  }

  Bsr bsr("", crs.numRows(), crs.numCols(), crs.nnz(), bsrVals, bsrRowMap, bsrIndices, blockSize);
  return bsr;
}  // expand_crs_to_bsr

/*! \brief convert a crs already in block format to a Bsr matrix
 */
template <typename Bsr, typename Crs>
Bsr blocked_crs_to_bsr(const Crs &crs, size_t blockSize) {
  using bsr_value_type   = typename Bsr::value_type;
  using bsr_ordinal_type = typename Bsr::ordinal_type;
  using crs_size_type    = typename Crs::non_const_size_type;

  // copy matrix data to host
  auto hRowMap  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crs.graph.row_map);
  auto hColInds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crs.graph.entries);
  auto hVals    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), crs.values);
  Kokkos::fence();

  // construct COO data on host
  std::vector<bsr_value_type> vals;
  std::vector<bsr_ordinal_type> rows, cols;

  vals.reserve(crs.nnz());
  rows.reserve(crs.nnz());
  cols.reserve(crs.nnz());

  for (bsr_ordinal_type row = 0; row < bsr_ordinal_type(hRowMap.size()) - 1; ++row) {
    for (crs_size_type ci = hRowMap(row); ci < hRowMap(row + 1); ++ci) {
      bsr_ordinal_type col = hColInds(ci);
      bsr_value_type val   = hVals(ci);

      rows.push_back(row);
      cols.push_back(col);
      vals.push_back(val);
    }
  }

  Bsr bsr("", crs.numRows(), crs.numCols(), crs.nnz(), vals.data(), rows.data(), cols.data(), blockSize);
  return bsr;
}  // expand_crs_to_bsr

}  // namespace Impl
}  // namespace KokkosSparse
