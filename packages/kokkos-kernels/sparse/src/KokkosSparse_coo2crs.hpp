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

#ifndef _KOKKOSSPARSE_COO2CRS_HPP
#define _KOKKOSSPARSE_COO2CRS_HPP

#include "KokkosSparse_CooMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Utils.hpp"

#include "KokkosSparse_coo2crs_impl.hpp"

namespace KokkosSparse {
// clang-format off
///
/// \brief Blocking function that converts a CooMatrix into a CrsMatrix. Values are summed.
/// \tparam DimType the dimension type
/// \tparam RowViewType The row array view type
/// \tparam ColViewType The column array view type
/// \tparam DataViewType The data array view type
/// \param m the number of rows
/// \param n the number of columns
/// \param row the array of row ids
/// \param col the array of col ids
/// \param data the array of data
/// \return A KokkosSparse::CrsMatrix.
// clang-format on
template <class DimType, class RowViewType, class ColViewType, class DataViewType>
auto coo2crs(DimType m, DimType n, RowViewType row, ColViewType col, DataViewType data) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<RowViewType>::value, "RowViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<ColViewType>::value, "CalViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<DataViewType>::value, "DataViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(RowViewType::rank) == 1, "RowViewType must have rank 1.");
  static_assert(static_cast<int>(ColViewType::rank) == 1, "ColViewType must have rank 1.");
  static_assert(static_cast<int>(DataViewType::rank) == 1, "DataViewType must have rank 1.");
#endif

  static_assert(std::is_integral<typename RowViewType::value_type>::value,
                "RowViewType::value_type must be an integral.");
  static_assert(std::is_integral<typename ColViewType::value_type>::value,
                "ColViewType::value_type must be an integral.");

  if (row.extent(0) != col.extent(0) || row.extent(0) != data.extent(0))
    Kokkos::abort("row.extent(0) = col.extent(0) = data.extent(0) required.");

  if constexpr (std::is_signed_v<DimType>) {
    if (m < 0 || n < 0) Kokkos::abort("m >= 0 and n >= 0 required.");
  }

  using Coo2crsType = Impl::Coo2Crs<DimType, RowViewType, ColViewType, DataViewType, true>;
  Coo2crsType Coo2Crs(m, n, row, col, data);
  return Coo2Crs.get_crsMat();
}

// clang-format off
///
/// \brief Blocking function that converts a CooMatrix into a CrsMatrix. Values are summed.
/// \tparam ScalarType   The `KokkosSparse::CooMatrix::scalar_type`
/// \tparam OrdinalType  The KokkosSparse::CooMatrix::ordinal_type
/// \tparam DeviceType   The KokkosSparse::CooMatrix::device_type
/// \tparam MemoryTraits The KokkosSparse::CooMatrix::memory_traits
/// \tparam SizeType     The KokkosSparse::CooMatrix::size_type
/// \param cooMatrix     The sparse matrix stored in coordinate ("Coo") format.
/// \return A KokkosSparse::CrsMatrix.
// clang-format on
template <typename ScalarType, typename OrdinalType, class DeviceType, class MemoryTraitsType, typename SizeType>
auto coo2crs(KokkosSparse::CooMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &cooMatrix) {
  return coo2crs(cooMatrix.numRows(), cooMatrix.numCols(), cooMatrix.row(), cooMatrix.col(), cooMatrix.data());
}
}  // namespace KokkosSparse
#endif  //  _KOKKOSSPARSE_COO2CRS_HPP
