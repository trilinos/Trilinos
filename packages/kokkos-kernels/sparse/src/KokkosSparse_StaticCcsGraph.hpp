// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSSPARSE_STATICCCSGRAPH_HPP_
#define KOKKOSSPARSE_STATICCCSGRAPH_HPP_

#include <Kokkos_Core.hpp>

namespace KokkosSparse {

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
/// \tparam Arg3Type The fourth template parameter, which if provided
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
          typename SizeType = typename Kokkos::ViewTraits<DataType*, Arg1Type, Arg2Type, Arg3Type>::size_type>
class StaticCcsGraph {
 private:
  using traits = Kokkos::ViewTraits<DataType*, Arg1Type, Arg2Type, Arg3Type>;

 public:
  using data_type       = DataType;
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  using size_type       = SizeType;

  using col_map_type   = Kokkos::View<const size_type*, array_layout, device_type, memory_traits>;
  using entries_type   = Kokkos::View<data_type*, array_layout, device_type, memory_traits>;
  using row_block_type = Kokkos::View<const size_type*, array_layout, device_type, memory_traits>;

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

}  // namespace KokkosSparse

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Kokkos {
template <class DataType, class Arg1Type, class Arg2Type = void, class Arg3Type = void,
          typename SizeType = typename Kokkos::ViewTraits<DataType*, Arg1Type, Arg2Type, Arg3Type>::size_type>
using StaticCcsGraph = KokkosSparse::StaticCcsGraph<DataType, Arg1Type, Arg2Type, Arg3Type, SizeType>;
}
#endif

#endif
