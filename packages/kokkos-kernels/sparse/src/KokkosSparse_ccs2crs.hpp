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

#ifndef _KOKKOSSPARSE_CCS2CRS_HPP
#define _KOKKOSSPARSE_CCS2CRS_HPP
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

  OrdinalType __nrows;
  OrdinalType __ncols;
  SizeType __nnz;
  ValViewType __vals;
  ColMapViewType __col_map;
  RowIdViewType __row_ids;

  RowIdViewType __crs_row_cnt;

  CrsValsViewType __crs_vals;
  CrsRowMapViewType __crs_row_map;
  CrsColIdViewType __crs_col_ids;

 public:
  Ccs2Crs(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, ColMapViewType col_map,
          RowIdViewType row_ids)
      : __nrows(nrows), __ncols(ncols), __nnz(nnz), __vals(vals), __col_map(col_map), __row_ids(row_ids) {
    __crs_vals    = CrsValsViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "__crs_vals"), nnz);
    __crs_row_map = CrsRowMapViewType(Kokkos::view_alloc("__crs_row_map"), nrows + 1);
    __crs_col_ids = CrsColIdViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "__crs_col_ids"), nnz);

    KokkosSparse::Impl::transpose_matrix<ColMapViewType, RowIdViewType, ValViewType, CrsRowMapViewType,
                                         CrsColIdViewType, CrsValsViewType, CrsRowMapViewType, CrsET>(
        __ncols, __nrows, __col_map, __row_ids, __vals, __crs_row_map, __crs_col_ids, __crs_vals);
  }

  CrsType get_crsMat() { return CrsType("ccs2crs", __nrows, __ncols, __nnz, __crs_vals, __crs_row_map, __crs_col_ids); }
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
#endif  //  _KOKKOSSPARSE_CCS2CRS_HPP
