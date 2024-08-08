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

/// \file KokkosSparse_CooMatrix.hpp
/// \brief Local sparse matrix interface
///
/// This file provides KokkosSparse::CooMatrix.  This implements a
/// local (no MPI) sparse matrix stored in coordinate ("Coo") format
/// which is also known as ivj or triplet format.

#ifndef KOKKOS_SPARSE_COOMATRIX_HPP_
#define KOKKOS_SPARSE_COOMATRIX_HPP_

#include "Kokkos_Core.hpp"
#include "KokkosKernels_Error.hpp"
#include <sstream>
#include <stdexcept>

namespace KokkosSparse {
/// \class CooMatrix
///
/// \brief Coordinate format implementation of a sparse matrix.
///
/// \tparam ScalarType The type of scalar entries in the sparse matrix.
/// \tparam OrdinalType The type of index entries in the sparse matrix.
/// \tparam Device The Kokkos Device type.
/// \tparam MemoryTraits Traits describing how Kokkos manages and
///   accesses data.  The default parameter suffices for most users.
/// "Coo" stands for "coordinate format".
template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
          class SizeType = typename Kokkos::ViewTraits<OrdinalType*, Device, void, void>::size_type>
class CooMatrix {
 public:
  //! Type of each value in the matrix
  using scalar_type = ScalarType;
  //! Type of each value in the const matrix
  using const_scalar_type = const std::remove_const_t<scalar_type>;
  //! Non constant scalar type
  using non_const_scalar_type = std::remove_const_t<scalar_type>;
  //! Type of each index in the matrix
  using ordinal_type = OrdinalType;
  //! Type of each value in the const matrix
  using const_ordinal_type = const std::remove_const_t<ordinal_type>;
  //! Non constant ordinal type
  using non_const_ordinal_type = std::remove_const_t<ordinal_type>;
  //! Type of each row index in the matrix
  using row_type = ordinal_type;
  //! Type of each column index in the matrix
  using column_type = ordinal_type;
  //! Type of the Kokkos::Device
  using device_type = Device;
  //! Type of the Kokkos::Device::execution_space
  using execution_space = typename device_type::execution_space;
  //! Type of the Kokkos::Device::memory_space
  using memory_space = typename device_type::memory_space;
  //! Type of the Kokkos::MemoryTraits
  using memory_traits = MemoryTraits;
  //! Type of all integral class members
  using size_type = SizeType;

  static_assert(std::is_integral_v<OrdinalType>, "OrdinalType must be an integral.");

  //! The type of the row index view in the matrix
  using row_view = Kokkos::View<row_type*, Kokkos::LayoutRight, device_type, memory_traits>;
  //! The type of the column index view in the matrix
  using column_view = Kokkos::View<column_type*, Kokkos::LayoutRight, device_type, memory_traits>;
  //! The type of the scalar values view in the matrix
  using scalar_view = Kokkos::View<scalar_type*, Kokkos::LayoutRight, device_type, memory_traits>;
  //! The type of a constant CooMatrix
  using const_type = CooMatrix<const_scalar_type, const_ordinal_type, device_type, memory_traits, size_type>;

 private:
  size_type m_num_rows, m_num_cols;
  row_view m_row;
  column_view m_col;
  scalar_view m_data;

 public:
  /// \brief Default constructor; constructs an empty sparse matrix.
  KOKKOS_INLINE_FUNCTION
  CooMatrix() : m_num_rows(0), m_num_cols(0) {}

  // clang-format off
  /// \brief Constructor that accepts a column indicies view, row indices view, and
  ///        values view.
  ///
  /// The matrix will store and use the column indices, rows indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param nrows   [in] The number of rows.
  /// \param ncols   [in] The number of columns.
  /// \param row_in  [in] The row indexes.
  /// \param col_in  [in] The column indexes.
  /// \param data_in [in] The values.
  // clang-format on
  CooMatrix(size_type nrows, size_type ncols, row_view row_in, column_view col_in, scalar_view data_in)
      : m_num_rows(nrows), m_num_cols(ncols), m_row(row_in), m_col(col_in), m_data(data_in) {
    if (m_data.extent(0) != m_row.extent(0) || m_row.extent(0) != m_col.extent(0)) {
      std::ostringstream os;
      os << "data.extent(0): " << m_data.extent(0) << " != "
         << "row.extent(0): " << m_row.extent(0) << " != "
         << "col.extent(0): " << m_col.extent(0) << ".";
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type numCols() const { return m_num_cols; }

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type numRows() const { return m_num_rows; }

  //! The number of stored entries in the sparse matrix, including zeros.
  KOKKOS_INLINE_FUNCTION size_type nnz() const { return m_data.extent(0); }

  //! The row indexes of the matrix
  KOKKOS_INLINE_FUNCTION row_view row() const { return m_row; }

  //! The column indexes of the matrix
  KOKKOS_INLINE_FUNCTION column_view col() const { return m_col; }

  //! The scalar values of the matrix
  KOKKOS_INLINE_FUNCTION scalar_view data() const { return m_data; }
};

/// \class is_coo_matrix
/// \brief is_coo_matrix<T>::value is true if T is a CooMatrix<...>, false
/// otherwise
template <typename>
struct is_coo_matrix : public std::false_type {};
template <typename... P>
struct is_coo_matrix<CooMatrix<P...>> : public std::true_type {};
template <typename... P>
struct is_coo_matrix<const CooMatrix<P...>> : public std::true_type {};

}  // namespace KokkosSparse
#endif
