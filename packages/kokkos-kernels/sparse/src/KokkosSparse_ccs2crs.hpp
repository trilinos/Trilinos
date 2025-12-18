// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_CcsMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifndef KOKKOSSPARSE_CCS2CRS_HPP
#define KOKKOSSPARSE_CCS2CRS_HPP
namespace KokkosSparse {
namespace Impl {
template <class OrdinalType, class SizeType, class ValViewType, class ColMapViewType, class RowIdViewType>
class Ccs2Crs {
 private:
  using CrsST             = typename ValViewType::value_type;
  using CrsOT             = OrdinalType;
  using CrsET             = typename ValViewType::execution_space;
  using CrsMT             = void;
  using CrsSzT            = SizeType;
  using CrsType           = CrsMatrix<CrsST, CrsOT, CrsET, CrsMT, CrsSzT>;
  using CrsValsViewType   = typename CrsType::values_type;
  using CrsRowMapViewType = typename CrsType::row_map_type::non_const_type;
  using CrsColIdViewType  = typename CrsType::index_type;

  OrdinalType nrows_;
  OrdinalType ncols_;
  SizeType nnz_;
  ValViewType vals_;
  ColMapViewType col_map_;

  CrsValsViewType crs_vals_;
  CrsRowMapViewType crs_row_map_;
  CrsColIdViewType crs_col_ids_;

 public:
  Ccs2Crs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, ColMapViewType col_map,
          RowIdViewType row_ids)
      : nrows_(nrows), ncols_(ncols), nnz_(nnz), vals_(vals), col_map_(col_map) {
    crs_vals_    = CrsValsViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "crs_vals_"), nnz);
    crs_row_map_ = CrsRowMapViewType(Kokkos::view_alloc("crs_row_map_"), nrows + 1);
    crs_col_ids_ = CrsColIdViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "crs_col_ids_"), nnz);

    KokkosSparse::Impl::transpose_matrix<ColMapViewType, RowIdViewType, ValViewType, CrsRowMapViewType,
                                         CrsColIdViewType, CrsValsViewType, CrsRowMapViewType, CrsET>(
        ncols_, nrows_, col_map_, row_ids, vals_, crs_row_map_, crs_col_ids_, crs_vals_);
  }

  CrsType get_crsMat() { return CrsType("ccs2crs", nrows_, ncols_, nnz_, crs_vals_, crs_row_map_, crs_col_ids_); }
};
}  // namespace Impl
// clang-format off
///
/// \brief Blocking function that converts a ccs matrix to a CrsMatrix.
/// Ccs values are copied from column-contiguous layout into row-contiguous layout.
/// \tparam OrdinalType    The view value type associated with the RowIdViewType
/// \tparam SizeType       The type of nnz
/// \tparam ValViewType    The values view type
/// \tparam ColMapViewType The column map view type
/// \tparam RowIdViewType  The row ids view type
/// \param nrows   The number of rows in the ccs matrix
/// \param ncols   The number of columns in the ccs matrix
/// \param nnz     The number of non-zeros in the ccs matrix
/// \param vals    The values view of the ccs matrix
/// \param col_map The column map view of the ccs matrix
/// \param row_ids The row ids view of the ccs matrix
/// \return A KokkosSparse::CrsMatrix.
///
/// \note In KokkosKernels sparse code, adj stands for adjacency list
///   and here we're passing in a ccs matrix with xadj=col_map and adj=row_ids.
// clang-format on
template <class OrdinalType, class SizeType, class ValViewType, class ColMapViewType, class RowIdViewType>
auto ccs2crs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, ColMapViewType col_map,
             RowIdViewType row_ids) {
  static_assert(std::is_same_v<SizeType, typename ColMapViewType::non_const_value_type>,
                "ccs2crs: SizeType (type of nnz) must match the element type of "
                "ColMapViewType");
  static_assert(std::is_same_v<OrdinalType, typename RowIdViewType::non_const_value_type>,
                "ccs2crs: OrdinalType (type of nrows, ncols) must match the element type "
                "of RowIdViewType");
  using Ccs2crsType = Impl::Ccs2Crs<OrdinalType, SizeType, ValViewType, ColMapViewType, RowIdViewType>;
  Ccs2crsType ccs2Crs(nrows, ncols, nnz, vals, col_map, row_ids);
  return ccs2Crs.get_crsMat();
}

///
/// @brief Blocking function that converts a crs matrix to a CcsMatrix.
/// Ccs values are copied from column-contiguous layout into row-contiguous
/// layout.
///
/// \tparam ScalarType   The ccsMatrix::scalar_type
/// \tparam OrdinalType  The ccsMatrix::ordinal_type
/// \tparam DeviceType   The ccsMatrix::device_type
/// \tparam MemoryTraits The ccsMatrix::memory_traits
/// \tparam SizeType     The ccsMatrix::size_type
/// \param ccsMatrix The KokkosSparse::CcsMatrix.
/// \return A KokkosSparse::CrsMatrix.
template <typename ScalarType, typename OrdinalType, class DeviceType, class MemoryTraitsType, typename SizeType>
auto ccs2crs(KokkosSparse::CcsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &ccsMatrix) {
  return ccs2crs(ccsMatrix.numRows(), ccsMatrix.numCols(), ccsMatrix.nnz(), ccsMatrix.values, ccsMatrix.graph.col_map,
                 ccsMatrix.graph.entries);
}
}  // namespace KokkosSparse
#endif  //  KOKKOSSPARSE_CCS2CRS_HPP
