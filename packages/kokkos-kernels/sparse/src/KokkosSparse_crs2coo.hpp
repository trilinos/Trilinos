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
#include "KokkosSparse_CooMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifndef _KOKKOSSPARSE_CRS2COO_HPP
#define _KOKKOSSPARSE_CRS2COO_HPP
namespace KokkosSparse {
namespace Impl {
template <class OrdinalType, class SizeType, class ValViewType, class RowMapViewType, class ColIdViewType,
          class DeviceType = typename ValViewType::execution_space>
class Crs2Coo {
 private:
  using scalar_type           = typename ValViewType::value_type;
  using const_scalar_type     = const std::remove_const_t<scalar_type>;
  using non_const_scalar_type = std::remove_const_t<scalar_type>;

  using ordinal_type           = OrdinalType;
  using const_ordinal_type     = const std::remove_const_t<ordinal_type>;
  using non_const_ordinal_type = std::remove_const_t<ordinal_type>;

  using size_type           = SizeType;
  using const_size_type     = const std::remove_const_t<size_type>;
  using non_const_size_type = std::remove_const_t<size_type>;

  using device_type = DeviceType;

  using row_view                = typename Kokkos::View<ordinal_type *, device_type>;
  using col_view                = row_view;
  using non_const_coo_data_view = typename ValViewType::non_const_type;
  using coo_type                = CooMatrix<scalar_type, ordinal_type, device_type>;

  non_const_ordinal_type m_nrows;
  non_const_ordinal_type m_ncols;
  non_const_size_type m_nnz;

  non_const_coo_data_view m_data;
  col_view m_col;
  row_view m_row;

  ValViewType m_vals;
  RowMapViewType m_row_map;
  ColIdViewType m_col_ids;

  using copy_tp1_pt          = Kokkos::TeamPolicy<DeviceType>;
  using copy_tp1_member_type = typename copy_tp1_pt::member_type;

 public:
  Crs2Coo(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map,
          ColIdViewType col_ids)
      : m_nrows(nrows), m_ncols(ncols), m_nnz(nnz), m_vals(vals), m_row_map(row_map), m_col_ids(col_ids) {
    m_data = non_const_coo_data_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_data"), nnz);
    m_col  = col_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_col"), nnz);
    m_row  = row_view(Kokkos::view_alloc(Kokkos::WithoutInitializing, "m_row"), nnz);

    copy_tp1_pt policy(m_nrows, 1, 1);
    {
      auto vec_len_max = policy.vector_length_max();
      copy_tp1_pt query_policy(m_nrows, 1, vec_len_max);
      policy = copy_tp1_pt(m_nrows, query_policy.team_size_recommended(*this, Kokkos::ParallelForTag()), vec_len_max);
    }

    Kokkos::parallel_for("Crs2Coo", policy, *this);
    DeviceType().fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const copy_tp1_member_type &member) const {
    auto i         = member.league_rank();
    auto row_start = m_row_map(i);
    auto row_len   = m_row_map(i + 1) - row_start;
    auto row_end   = row_start + row_len;

    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, row_start, row_end), [&](const size_type &id) {
      m_data(id) = m_vals(id);
      m_col(id)  = m_col_ids(id);
      m_row(id)  = i;
    });
  }

  coo_type get_cooMat() { return coo_type(m_nrows, m_ncols, m_row, m_col, m_data); }
};
}  // namespace Impl
// clang-format off
///
/// \brief Blocking function that converts a CrsMatrix to a CooMatrix.
/// Crs values are copied into the CooMatrix in the order they appear
/// within the CrsMatrix, starting from row 0 to row nrows - 1.
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
/// \return A KokkosSparse::CooMatrix.
///
// clang-format on
template <class OrdinalType, class SizeType, class ValViewType, class RowMapViewType, class ColIdViewType>
auto crs2coo(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals, RowMapViewType row_map,
             ColIdViewType col_ids) {
  static_assert(std::is_same_v<SizeType, typename RowMapViewType::non_const_value_type>,
                "crs2coo: SizeType (type of nnz) must match the element type of "
                "RowMapViewType");
  static_assert(std::is_same_v<OrdinalType, typename ColIdViewType::non_const_value_type>,
                "crs2coo: OrdinalType (type of nrows, ncols) must match the element type "
                "of ColIdViewType");
  using Crs2cooType = Impl::Crs2Coo<OrdinalType, SizeType, ValViewType, RowMapViewType, ColIdViewType>;
  Crs2cooType crs2Coo(nrows, ncols, nnz, vals, row_map, col_ids);
  return crs2Coo.get_cooMat();
}

///
/// @brief Blocking function that converts a CrsMatrix to a CooMatrix.
/// Crs values are copied into the CooMatrix in the order they appear
/// within the CrsMatrix, starting from row 0 to row nrows - 1.
///
/// \tparam ScalarType   The crsMatrix::scalar_type
/// \tparam OrdinalType  The crsMatrix::ordinal_type
/// \tparam DeviceType   The crsMatrix::device_type
/// \tparam MemoryTraits The crsMatrix::memory_traits
/// \tparam SizeType     The crsMatrix::size_type
/// \param crsMatrix The KokkosSparse::CrsMatrix.
/// \return A KokkosSparse::CooMatrix.
template <typename ScalarType, typename OrdinalType, class DeviceType, class MemoryTraitsType, typename SizeType>
auto crs2coo(KokkosSparse::CrsMatrix<ScalarType, OrdinalType, DeviceType, MemoryTraitsType, SizeType> &crsMatrix) {
  return crs2coo(crsMatrix.numRows(), crsMatrix.numCols(), crsMatrix.nnz(), crsMatrix.values, crsMatrix.graph.row_map,
                 crsMatrix.graph.entries);
}
}  // namespace KokkosSparse
#endif  //  _KOKKOSSPARSE_CRS2COO_HPP
