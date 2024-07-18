// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_GETGRAPHOFFRANKOFFSETS_DECL_HPP
#define TPETRA_DETAILS_GETGRAPHOFFRANKOFFSETS_DECL_HPP

/// \file Tpetra_Details_getGraphOffRankOffsets_decl.hpp
/// \brief Declare and define the function
///   Tpetra::Details::getGraphOffRankOffsets, an implementation detail
///   of Tpetra::CrsGraph.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_StaticCrsGraph.hpp"
#include "Tpetra_Details_LocalMap.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Impl {

/// \brief Implementation detail of
///   Tpetra::Details::getGraphOffRankOffsets, which in turn is an
///   implementation detail of Tpetra::CrsGraph.
///
/// FIXME (mfh 12 Mar 2016) There's currently no way to make a
/// MemoryUnmanaged Kokkos::StaticCrsGraph.  Thus, we have to do this
/// separately for its column indices.  We want the column indices to
/// be unmanaged because we need to take subviews in this kernel.
/// Taking a subview of a managed View updates the reference count,
/// which is a thread scalability bottleneck.
///
/// mfh 12 Mar 2016: Tpetra::CrsGraph::getLocalOffRankOffsets returns
/// offsets as size_t.  However, see Github Issue #213.
template<class LO,
         class GO,
         class DeviceType,
         class OffsetType = size_t>
class GetGraphOffRankOffsets {
public:
  typedef typename DeviceType::device_type device_type;
  typedef OffsetType offset_type;
  typedef ::Kokkos::View<offset_type*,
                         device_type,
                         ::Kokkos::MemoryUnmanaged> offsets_type;
  typedef ::Kokkos::StaticCrsGraph<LO,
                                   ::Kokkos::LayoutLeft,
                                   device_type,
                                   void, size_t> local_graph_type;
  typedef ::Tpetra::Details::LocalMap<LO, GO, device_type> local_map_type;
  typedef ::Kokkos::View<const typename local_graph_type::size_type*,
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
  GetGraphOffRankOffsets (const offsets_type& OffRankOffsets,
                          const local_map_type& lclColMap,
                          const local_map_type& lclDomMap,
                          const row_offsets_type& ptr,
                          const lcl_col_inds_type& ind);

  //! Kokkos::parallel_for loop body.
  KOKKOS_FUNCTION void operator() (const LO& lclRowInd) const;

private:
  offsets_type OffRankOffsets_;
  local_map_type lclColMap_;
  local_map_type lclDomMap_;
  row_offsets_type ptr_;
  lcl_col_inds_type ind_;
  LO lclNumRows_;
};

} // namespace Impl

template<class OffsetsType,
         class LclMapType,
         class LclGraphType>
void
getGraphOffRankOffsets (const OffsetsType& OffRankOffsets,
                        const LclMapType& lclColMap,
                        const LclMapType& lclDomMap,
                        const LclGraphType& lclGraph)
{
  typedef typename OffsetsType::non_const_value_type offset_type;
  typedef typename LclMapType::local_ordinal_type LO;
  typedef typename LclMapType::global_ordinal_type GO;
  typedef typename LclMapType::device_type DT;

  typedef Impl::GetGraphOffRankOffsets<LO, GO, DT, offset_type> impl_type;

  // The functor's constructor runs the functor.
  impl_type impl (OffRankOffsets, lclColMap, lclDomMap, lclGraph.row_map, lclGraph.entries);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETGRAPHOFFRANKOFFSETS_DECL_HPP
