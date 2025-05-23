// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GETGRAPHDIAGOFFSETS_DECL_HPP
#define TPETRA_DETAILS_GETGRAPHDIAGOFFSETS_DECL_HPP

/// \file Tpetra_Details_getGraphDiagOffsets_decl.hpp
/// \brief Declare and define the function
///   Tpetra::Details::getGraphDiagOffsets, an implementation detail
///   of Tpetra::CrsGraph.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "KokkosSparse_StaticCrsGraph.hpp"
#include "Tpetra_Details_LocalMap.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Impl {

/// \brief Implementation detail of
///   Tpetra::Details::getGraphDiagOffsets, which in turn is an
///   implementation detail of Tpetra::CrsGraph.
///
/// FIXME (mfh 12 Mar 2016) There's currently no way to make a
/// MemoryUnmanaged KokkosSparse::StaticCrsGraph.  Thus, we have to do this
/// separately for its column indices.  We want the column indices to
/// be unmanaged because we need to take subviews in this kernel.
/// Taking a subview of a managed View updates the reference count,
/// which is a thread scalability bottleneck.
///
/// mfh 12 Mar 2016: Tpetra::CrsGraph::getLocalDiagOffsets returns
/// offsets as size_t.  However, see Github Issue #213.
template<class LO,
         class GO,
         class DeviceType,
         class DiagOffsetType = size_t>
class GetGraphDiagOffsets {
public:
  typedef typename DeviceType::device_type device_type;
  typedef DiagOffsetType diag_offset_type;
  typedef ::Kokkos::View<diag_offset_type*,
                         device_type,
                         ::Kokkos::MemoryUnmanaged> diag_offsets_type;
  using local_graph_device_type = ::KokkosSparse::StaticCrsGraph<LO,
                                  ::Kokkos::LayoutLeft,
                                  device_type, void, size_t>;
  typedef ::Tpetra::Details::LocalMap<LO, GO, device_type> local_map_type;
  typedef ::Kokkos::View<const typename local_graph_device_type::size_type*,
                         ::Kokkos::LayoutLeft,
                         device_type,
                         ::Kokkos::MemoryUnmanaged> row_offsets_type;
  // This is unmanaged for performance, because we need to take
  // subviews inside the functor.
  typedef ::Kokkos::View<const LO*,
                         ::Kokkos::LayoutLeft,
                         device_type,
                         ::Kokkos::MemoryUnmanaged> lcl_col_inds_type;

  //! Constructor; also runs the functor.
  GetGraphDiagOffsets (const diag_offsets_type& diagOffsets,
                       const local_map_type& lclRowMap,
                       const local_map_type& lclColMap,
                       const row_offsets_type& ptr,
                       const lcl_col_inds_type& ind,
                       const bool isSorted);

  //! Kokkos::parallel_for loop body.
  KOKKOS_FUNCTION void operator() (const LO& lclRowInd) const;

private:
  diag_offsets_type diagOffsets_;
  local_map_type lclRowMap_;
  local_map_type lclColMap_;
  row_offsets_type ptr_;
  lcl_col_inds_type ind_;
  bool isSorted_;
};

} // namespace Impl

template<class DiagOffsetsType,
         class LclMapType,
         class RowOffsetsType,
         class LclColIndsType>
void
getGraphDiagOffsets (const DiagOffsetsType& diagOffsets,
                     const LclMapType& lclRowMap,
                     const LclMapType& lclColMap,
                     const RowOffsetsType& ptr,
                     const LclColIndsType& ind,
                     const bool isSorted)
{
  static_assert (Kokkos::is_view<DiagOffsetsType>::value,
                 "DiagOffsetsType (the type of diagOffsets) must be a Kokkos::View.");
  static_assert (Kokkos::is_view<RowOffsetsType>::value,
                 "RowOffsetsType (the type of ptr) must be a Kokkos::View.");
  static_assert (Kokkos::is_view<LclColIndsType>::value,
                 "LclColIndsType (the type of ind) must be a Kokkos::View.");
  static_assert (static_cast<int> (DiagOffsetsType::rank) == 1,
                 "DiagOffsetsType (the type of diagOffsets) must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (RowOffsetsType::rank) == 1,
                 "RowOffsetsType (the type of ptr) must be a rank-1 Kokkos::View.");
  static_assert (static_cast<int> (LclColIndsType::rank) == 1,
                 "LclColIndsType (the type of ind) must be a rank-1 Kokkos::View.");
  typedef typename DiagOffsetsType::non_const_value_type diag_offset_type;
  static_assert (std::is_same<typename DiagOffsetsType::value_type, diag_offset_type>::value,
                 "DiagOffsetsType (the type of diagOffsets) must be a nonconst "
                 "Kokkos::View, since it is the output argument of this function.");
  static_assert (std::is_integral<diag_offset_type>::value,
                 "The type of each entry of diagOffsets must be an integer.");
  typedef typename LclColIndsType::non_const_value_type local_ordinal_type;
  static_assert (std::is_integral<local_ordinal_type>::value,
                 "The type of each entry of ind (the array of column indices) "
                 "must be an integer.");
  static_assert (sizeof (diag_offset_type) >= sizeof (local_ordinal_type),
                 "Diagonal offset type must be big enough to count the number "
                 "of column indices, since the diagonal entry (if it exists) "
                 "may be anywhere in a row.");
  typedef typename RowOffsetsType::non_const_value_type row_offset_type;
  static_assert (std::is_integral<row_offset_type>::value,
                 "The type of each entry of ptr must be an integer.");

  typedef typename LclMapType::local_ordinal_type LO;
  typedef typename LclMapType::global_ordinal_type GO;
  typedef typename LclMapType::device_type DT;

  typedef Impl::GetGraphDiagOffsets<LO, GO, DT, diag_offset_type> impl_type;
  // The functor's constructor runs the functor.
  impl_type impl (diagOffsets, lclRowMap, lclColMap, ptr, ind, isSorted);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETGRAPHDIAGOFFSETS_DECL_HPP
