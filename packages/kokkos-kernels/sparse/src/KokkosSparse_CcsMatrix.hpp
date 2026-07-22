// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosSparse_CcsMatrix.hpp
/// \brief Local sparse matrix interface
///
/// This file provides KokkosSparse::CcsMatrix.  This implements a
/// local (no MPI) sparse matrix stored in compressed column sparse
/// ("Ccs") format.

#ifndef KOKKOSSPARSE_CCSMATRIX_HPP_
#define KOKKOSSPARSE_CCSMATRIX_HPP_

#include "Kokkos_Core.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "KokkosSparse_findRelOffset.hpp"
#include "KokkosSparse_StaticCcsGraph.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Macros.hpp"

namespace KokkosSparse {
/// \class CcsMatrix
/// \brief Compressed sparse column implementation of a sparse matrix.
/// \tparam ScalarType The type of entries in the sparse matrix.
/// \tparam OrdinalType The type of column indices in the sparse matrix.
/// \tparam Device The Kokkos Device type.
/// \tparam MemoryTraits Traits describing how Kokkos manages and
///   accesses data.  The default parameter suffices for most users.
///
/// "Ccs" stands for "compressed column sparse."
template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
          class SizeType = typename Kokkos::ViewTraits<OrdinalType*, Device, void, void>::size_type>
class CcsMatrix {
  static_assert(std::is_signed<OrdinalType>::value, "CcsMatrix requires that OrdinalType is a signed integer type.");

 public:
  //! Type of the matrix's execution space.
  typedef typename Device::execution_space execution_space;
  //! Type of the matrix's memory space.
  typedef typename Device::memory_space memory_space;
  //! Canonical device type
  typedef Kokkos::Device<execution_space, memory_space> device_type;
  typedef MemoryTraits memory_traits;

  /// \brief Type of each entry of the "column map."
  ///
  /// The "column map" corresponds to the \c ptr array of column offsets in
  /// compressed sparse column (CCS) storage.
  typedef SizeType size_type;
  //! Type of each value in the matrix.
  typedef ScalarType value_type;
  //! Type of each (column) index in the matrix.
  typedef OrdinalType ordinal_type;
  //! Type of the graph structure of the sparse matrix - consistent with Kokkos.
  typedef StaticCcsGraph<ordinal_type, KokkosKernels::default_layout, device_type, memory_traits, size_type>
      staticccsgraph_type;
  //! Type of the "column map" (which contains the offset for each column's
  //! data).
  typedef typename staticccsgraph_type::col_map_type col_map_type;
  typedef Kokkos::View<value_type*, Kokkos::LayoutRight, device_type, MemoryTraits> values_type;
  //! Type of column indices in the sparse matrix.
  typedef typename staticccsgraph_type::entries_type index_type;

  /// \name Storage of the actual sparsity structure and values.
  ///
  /// CcsMatrix uses the compressed sparse column (CCS) storage format to
  /// store the sparse matrix.
  //@{
  //! The graph (sparsity structure) of the sparse matrix.
  staticccsgraph_type graph;
  //! The 1-D array of values of the sparse matrix.
  values_type values;
  //@}

 private:
  /// \brief The number of rows in the CCS matrix
  ordinal_type numRows_;

 public:
  /// \brief Default constructor; constructs an empty sparse matrix.
  KOKKOS_INLINE_FUNCTION
  CcsMatrix() : numRows_(0) {}

  // clang-format off
  /// \brief Constructor that accepts a column map, row indices, and
  ///   values.
  ///
  /// The matrix will store and use the column map, indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param nrows [in] The number of rows.
  /// \param ncols [in] The number of columns.
  /// \param annz [in] The number of entries.
  /// \param vals [in] The entries.
  /// \param colmap [in] The column map (containing the offsets to the data in
  /// each column).
  /// \param rows [in] The row indices.
  // clang-format on
  CcsMatrix(const std::string& /* label */, const OrdinalType nrows, const OrdinalType ncols, const size_type annz,
            const values_type& vals, const col_map_type& colmap, const index_type& rows)
      : graph(rows, colmap), values(vals), numRows_(nrows) {
    const ordinal_type actualNumRows = (colmap.extent(0) != 0)
                                           ? static_cast<ordinal_type>(colmap.extent(0) - static_cast<size_type>(1))
                                           : static_cast<ordinal_type>(0);
    if (ncols != actualNumRows) {
      std::ostringstream os;
      os << "Input argument ncols = " << ncols
         << " != the actual number of "
            "rows "
         << actualNumRows << " according to the 'rows' input argument.";
      throw std::invalid_argument(os.str());
    }
    if (annz != nnz()) {
      std::ostringstream os;
      os << "Input argument annz = " << annz << " != this->nnz () = " << nnz() << ".";
      throw std::invalid_argument(os.str());
    }
  }

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const { return graph.numCols(); }

  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const { return numRows_; }

  /// \brief Modify the number of rows in the sparse matrix.
  ///
  /// This invalidates any algorithm handles which previously used this matrix.
  void setNumRows(ordinal_type r) { numRows_ = r; }

  //! The number of "point" (non-block) rows in the matrix. Since Ccs is not
  //! blocked, this is just the number of regular rows.
  KOKKOS_INLINE_FUNCTION ordinal_type numPointRows() const { return numRows(); }

  //! The number of "point" (non-block) columns in the matrix. Since Ccs is not
  //! blocked, this is just the number of regular columns.
  KOKKOS_INLINE_FUNCTION ordinal_type numPointCols() const { return numCols(); }

  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type nnz() const { return graph.entries.extent(0); }
};

/// \class is_ccs_matrix
/// \brief is_ccs_matrix<T>::value is true if T is a CcsMatrix<...>, false
/// otherwise
template <typename>
struct is_ccs_matrix : public std::false_type {};
template <typename... P>
struct is_ccs_matrix<CcsMatrix<P...>> : public std::true_type {};
template <typename... P>
struct is_ccs_matrix<const CcsMatrix<P...>> : public std::true_type {};

}  // namespace KokkosSparse
#endif  // KOKKOSSPARSE_CCSMATRIX_HPP_
