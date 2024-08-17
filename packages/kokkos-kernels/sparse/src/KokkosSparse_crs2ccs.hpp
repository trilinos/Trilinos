//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_CcsMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifndef _KOKKOSSPARSE_CRS2CCS_HPP
#define _KOKKOSSPARSE_CRS2CCS_HPP
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

  OrdinalType __nrows;
  OrdinalType __ncols;
  SizeType __nnz;
  ValViewType __vals;
  RowMapViewType __row_map;
  ColIdViewType __col_ids;

  CcsValsViewType __ccs_vals;
  CcsColMapViewType __ccs_col_map;
  CcsRowIdViewType __ccs_row_ids;

 public:
  Crs2Ccs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map,
          ColIdViewType col_ids)
      : __nrows(nrows), __ncols(ncols), __nnz(nnz), __vals(vals), __row_map(row_map), __col_ids(col_ids) {
    __ccs_vals    = CcsValsViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "__ccs_vals"), nnz);
    __ccs_col_map = CcsColMapViewType(Kokkos::view_alloc("__ccs_col_map"), ncols + 1);
    __ccs_row_ids = CcsRowIdViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "__ccs_row_ids"), nnz);

    KokkosSparse::Impl::transpose_matrix<RowMapViewType, ColIdViewType, ValViewType, CcsColMapViewType,
                                         CcsRowIdViewType, CcsValsViewType, CcsColMapViewType, CcsET>(
        __nrows, __ncols, __row_map, __col_ids, __vals, __ccs_col_map, __ccs_row_ids, __ccs_vals);
  }

  CcsType get_ccsMat() { return CcsType("crs2ccs", __nrows, __ncols, __nnz, __ccs_vals, __ccs_col_map, __ccs_row_ids); }
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

#endif  //  _KOKKOSSPARSE_CRS2CCS_HPP
