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

/// \file KokkosSparse_CcsMatrix.hpp
/// \brief Local sparse matrix interface
///
/// This file provides KokkosSparse::CcsMatrix.  This implements a
/// local (no MPI) sparse matrix stored in compressed column sparse
/// ("Ccs") format.

#ifndef KOKKOS_SPARSE_CCSMATRIX_HPP_
#define KOKKOS_SPARSE_CCSMATRIX_HPP_

#include "Kokkos_Core.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "KokkosSparse_findRelOffset.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Macros.hpp"

namespace Kokkos {
/// \class StaticCcsGraph
/// \brief Compressed column storage array copied from Kokkos::StaticCrsGraph.
///
/// \tparam DataType The type of stored entries.  If a StaticCcsGraph is
///   used as the graph of a sparse matrix, then this is usually an
///   integer type, the type of the column indices in the sparse
///   matrix.
///
/// \tparam Arg1Type The second template parameter, corresponding
///   either to the Device type (if there are no more template
///   parameters) or to the Layout type (if there is at least one more
///   template parameter).
///
/// \tparam Arg2Type The third template parameter, which if provided
///   corresponds to the Device type.
///
/// \tparam Arg3Type The third template parameter, which if provided
///   corresponds to the MemoryTraits.
///
/// \tparam SizeType The type of col offsets.  Usually the default
///   parameter suffices.  However, setting a nondefault value is
///   necessary in some cases, for example, if you want to have a
///   sparse matrices with dimensions (and therefore column indices)
///   that fit in \c int, but want to store more than <tt>INT_MAX</tt>
///   entries in the sparse matrix.
///
/// A col has a range of entries:
/// <ul>
/// <li> <tt> col_map[i0] <= entry < col_map[i0+1] </tt> </li>
/// <li> <tt> 0 <= i1 < col_map[i0+1] - col_map[i0] </tt> </li>
/// <li> <tt> entries( entry ,            i2 , i3 , ... ); </tt> </li>
/// <li> <tt> entries( col_map[i0] + i1 , i2 , i3 , ... ); </tt> </li>
/// </ul>
template <class DataType, class Arg1Type, class Arg2Type = void, class Arg3Type = void,
          typename SizeType = typename ViewTraits<DataType*, Arg1Type, Arg2Type, Arg3Type>::size_type>
class StaticCcsGraph {
 private:
  using traits = ViewTraits<DataType*, Arg1Type, Arg2Type, Arg3Type>;

 public:
  using data_type       = DataType;
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  using size_type       = SizeType;

  using col_map_type   = View<const size_type*, array_layout, device_type, memory_traits>;
  using entries_type   = View<data_type*, array_layout, device_type, memory_traits>;
  using row_block_type = View<const size_type*, array_layout, device_type, memory_traits>;

  entries_type entries;
  col_map_type col_map;

  //! Construct an empty view.
  KOKKOS_INLINE_FUNCTION
  StaticCcsGraph() : entries(), col_map() {}

  //! Copy constructor (shallow copy).
  KOKKOS_INLINE_FUNCTION
  StaticCcsGraph(const StaticCcsGraph& rhs) : entries(rhs.entries), col_map(rhs.col_map) {}

  template <class EntriesType, class ColMapType>
  KOKKOS_INLINE_FUNCTION StaticCcsGraph(const EntriesType& entries_, const ColMapType& col_map_)
      : entries(entries_), col_map(col_map_) {}

  /**  \brief  Return number of columns in the graph
   */
  KOKKOS_INLINE_FUNCTION
  size_type numCols() const {
    return (col_map.extent(0) != 0) ? col_map.extent(0) - static_cast<size_type>(1) : static_cast<size_type>(0);
  }
};
}  // namespace Kokkos

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
  typedef Kokkos::StaticCcsGraph<ordinal_type, default_layout, device_type, memory_traits, size_type>
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
#endif
