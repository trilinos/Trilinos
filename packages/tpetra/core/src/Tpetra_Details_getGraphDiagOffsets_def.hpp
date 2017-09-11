/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_DETAILS_GETGRAPHDIAGOFFSETS_DEF_HPP
#define TPETRA_DETAILS_GETGRAPHDIAGOFFSETS_DEF_HPP

/// \file Tpetra_Details_getGraphDiagOffsets_def.hpp
/// \brief Define the implementation of the function
///   Tpetra::Details::getGraphDiagOffsets, an implementation detail
///   of Tpetra::CrsGraph.

#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Map.hpp"
#include "KokkosSparse_findRelOffset.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

/// \brief Implementation detail of
///   Tpetra::Details::getGraphDiagOffsets, which in turn is an
///   implementation detail of Tpetra::CrsGraph.
///
/// FIXME (mfh 12 Mar 2016) There's currently no way to make a
/// MemoryUnmanaged Kokkos::StaticCrsGraph.  Thus, we have to do this
/// separately for its column indices.  We want the column indices to
/// be unmanaged because we need to take subviews in this kernel.
/// Taking a subview of a managed View updates the reference count,
/// which is a thread scalability bottleneck.
///
/// mfh 12 Mar 2016: Tpetra::CrsGraph::getLocalDiagOffsets returns
/// offsets as size_t.  However, see Github Issue #213.
template<class LO,
         class GO,
         class Node,
         class DiagOffsetType>
GetGraphDiagOffsets<LO, GO, Node, DiagOffsetType>::
GetGraphDiagOffsets (const diag_offsets_type& diagOffsets,
                     const local_map_type& lclRowMap,
                     const local_map_type& lclColMap,
                     const row_offsets_type& ptr,
                     const lcl_col_inds_type& ind,
                     const bool isSorted) :
  diagOffsets_ (diagOffsets),
  lclRowMap_ (lclRowMap),
  lclColMap_ (lclColMap),
  ptr_ (ptr),
  ind_ (ind),
  isSorted_ (isSorted)
{
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> policy_type;

  const LO lclNumRows = lclRowMap.getNodeNumElements ();
  policy_type range (0, lclNumRows);
  Kokkos::parallel_for (range, *this);
}

template<class LO,
         class GO,
         class Node,
         class DiagOffsetType>
KOKKOS_FUNCTION void
GetGraphDiagOffsets<LO, GO, Node, DiagOffsetType>::
operator() (const LO& lclRowInd) const
{
  const size_t STINV =
    Tpetra::Details::OrdinalTraits<diag_offset_type>::invalid ();
  const GO gblRowInd = lclRowMap_.getGlobalElement (lclRowInd);
  const GO gblColInd = gblRowInd;
  const LO lclColInd = lclColMap_.getLocalElement (gblColInd);

  if (lclColInd == Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
    diagOffsets_[lclRowInd] = STINV;
  }
  else {
    // Could be empty, but that's OK.
    const LO numEnt = ptr_[lclRowInd+1] - ptr_[lclRowInd];
    // std::pair doesn't have its methods marked as device
    // functions, so we have to use Kokkos::pair.
    auto lclColInds =
      Kokkos::subview (ind_, Kokkos::make_pair (ptr_[lclRowInd],
                                                ptr_[lclRowInd+1]));
    using ::KokkosSparse::findRelOffset;
    const LO diagOffset =
      findRelOffset<LO, lcl_col_inds_type> (lclColInds, numEnt,
                                            lclColInd, 0, isSorted_);
    diagOffsets_[lclRowInd] = (diagOffset == numEnt) ? STINV :
      static_cast<diag_offset_type> (diagOffset);
  }
}

} // namespace Impl
} // namespace Details
} // namespace Tpetra

// Explicit template instantiation macro for
// Tpetra::Details::Impl::GetGraphDiagOffsets.  NOT FOR USERS!!!  Must
// be used inside the Tpetra namespace.
#define TPETRA_DETAILS_IMPL_GETGRAPHDIAGOFFSETS_INSTANT( LO, GO, NODE ) \
  template class Details::Impl::GetGraphDiagOffsets< LO, GO, NODE::device_type >;

#endif // TPETRA_DETAILS_GETGRAPHDIAGOFFSETS_DEF_HPP
