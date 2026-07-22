// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_CcsMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifndef KOKKOSSPARSE_CRS2CCS_HPP
#define KOKKOSSPARSE_CRS2CCS_HPP
namespace KokkosSparse {
namespace Impl {
template <class OrdinalType, class SizeType, class ValViewType, class RowMapViewType, class ColIdViewType>
class Crs2Ccs {
 private:
  using CcsST             = typename ValViewType::value_type;
  using CcsOT             = OrdinalType;
  using CcsET             = typename ValViewType::execution_space;
  using CcsMT             = void;
  using CcsSzT            = SizeType;
  using CcsType           = CcsMatrix<CcsST, CcsOT, CcsET, CcsMT, CcsSzT>;
  using CcsValsViewType   = typename CcsType::values_type;
  using CcsColMapViewType = typename CcsType::col_map_type::non_const_type;
  using CcsRowIdViewType  = typename CcsType::index_type;

  OrdinalType nrows_;
  OrdinalType ncols_;
  SizeType nnz_;
  ValViewType vals_;
  RowMapViewType row_map_;
  ColIdViewType col_ids_;

  CcsValsViewType ccs_vals_;
  CcsColMapViewType ccs_col_map_;
  CcsRowIdViewType ccs_row_ids_;

 public:
  Crs2Ccs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map,
          ColIdViewType col_ids)
      : nrows_(nrows), ncols_(ncols), nnz_(nnz), vals_(vals), row_map_(row_map), col_ids_(col_ids) {
    ccs_vals_    = CcsValsViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ccs_vals_"), nnz);
    ccs_col_map_ = CcsColMapViewType(Kokkos::view_alloc("ccs_col_map_"), ncols + 1);
    ccs_row_ids_ = CcsRowIdViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "ccs_row_ids_"), nnz);

    KokkosSparse::Impl::transpose_matrix<RowMapViewType, ColIdViewType, ValViewType, CcsColMapViewType,
                                         CcsRowIdViewType, CcsValsViewType, CcsColMapViewType, CcsET>(
        nrows_, ncols_, row_map_, col_ids_, vals_, ccs_col_map_, ccs_row_ids_, ccs_vals_);
  }

  CcsType get_ccsMat() { return CcsType("crs2ccs", nrows_, ncols_, nnz_, ccs_vals_, ccs_col_map_, ccs_row_ids_); }
};
}  // namespace Impl
// clang-format off
///
/// \brief Blocking function that converts a CrsMatrix to a CcsMatrix.
/// Crs values are copied from row-contiguous layout into column-contiguous layout.
/// \tparam OrdinalType The view value type associated with the RowIdViewType
/// \tparam SizeType The type of nnz
/// \tparam ValViewType    The values view type
/// \tparam RowMapViewType The column map view type
/// \tparam ColIdViewType  The row ids view type
/// \param nrows   The number of rows in the crs matrix
/// \param ncols   The number of columns in the crs matrix
/// \param nnz     The number of non-zeros in the crs matrix
/// \param vals    The values view of the crs matrix
/// \param row_map The row map view of the crs matrix
/// \param col_ids The col ids view of the crs matrix
/// \return A KokkosSparse::CcsMatrix.
///
/// \note In KokkosKernels sparse code, adj stands for adjacency list
///   and here we're passing in a crs matrix with xadj=row_map and adj=col_ids.
// clang-format on
template <class OrdinalType, class SizeType, class ValViewType, class RowMapViewType, class ColIdViewType>
auto crs2ccs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map,
             ColIdViewType col_ids) {
  static_assert(std::is_same_v<SizeType, typename RowMapViewType::non_const_value_type>,
                "crs2ccs: SizeType (type of nnz) must match the element type of "
                "RowMapViewType");
  static_assert(std::is_same_v<OrdinalType, typename ColIdViewType::non_const_value_type>,
                "crs2ccs: OrdinalType (type of nrows, ncols) must match the element type "
                "of ColIdViewType");
  using Crs2ccsType = Impl::Crs2Ccs<OrdinalType, SizeType, ValViewType, RowMapViewType, ColIdViewType>;
  Crs2ccsType crs2Ccs(nrows, ncols, nnz, vals, row_map, col_ids);
  return crs2Ccs.get_ccsMat();
}

///
/// @brief Blocking function that converts a crs matrix to a CcsMatrix.
/// Crs values are copied from row-contiguous layout into column-contiguous
/// layout.
///
/// \tparam ScalarType   The crsMatrix::scalar_type
/// \tparam OrdinalType  The crsMatrix::ordinal_type
/// \tparam DeviceType   The crsMatrix::device_type
/// \tparam MemoryTraits The crsMatrix::memory_traits
/// \tparam SizeType     The crsMatrix::size_type
/// \param crsMatrix The KokkosSparse::CrsMatrix.
/// \return A KokkosSparse::CcsMatrix.
template <typename ScalarType, typename OrdinalType, class DeviceType, class MemoryTraitsType, typename SizeType>
auto crs2ccs(KokkosSparse::CrsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &crsMatrix) {
  return crs2ccs(crsMatrix.numRows(), crsMatrix.numCols(), crsMatrix.nnz(), crsMatrix.values, crsMatrix.graph.row_map,
                 crsMatrix.graph.entries);
}
}  // namespace KokkosSparse

#endif  //  KOKKOSSPARSE_CRS2CCS_HPP
