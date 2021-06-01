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

#ifndef TPETRA_DETAILS_GETGRAPHOFFRANKOFFSETS_DEF_HPP
#define TPETRA_DETAILS_GETGRAPHOFFRANKOFFSETS_DEF_HPP

/// \file Tpetra_Details_getGraphOffRankOffsets_def.hpp
/// \brief Define the implementation of the function
///   Tpetra::Details::getGraphOffRankOffsets, an implementation detail
///   of Tpetra::CrsGraph.

#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Map.hpp"
#include "KokkosSparse_findRelOffset.hpp"

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
         class Node,
         class OffsetType>
GetGraphOffRankOffsets<LO, GO, Node, OffsetType>::
GetGraphOffRankOffsets (const offsets_type& OffRankOffsets,
                        const local_map_type& lclColMap,
                        const local_map_type& lclDomMap,
                        const row_offsets_type& ptr,
                        const lcl_col_inds_type& ind) :
  OffRankOffsets_ (OffRankOffsets),
  lclColMap_ (lclColMap),
  lclDomMap_ (lclDomMap),
  ptr_ (ptr),
  ind_ (ind)
{
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> policy_type;

  lclNumRows_ = ptr.extent(0)-1;
  policy_type range (0, ptr.extent(0));
  Kokkos::parallel_for (range, *this);
}

template<class LO,
         class GO,
         class Node,
         class OffsetType>
KOKKOS_FUNCTION void
GetGraphOffRankOffsets<LO, GO, Node, OffsetType>::
operator() (const LO& lclRowInd) const
{
  const LO INVALID =
    Tpetra::Details::OrdinalTraits<LO>::invalid ();

  if (lclRowInd == lclNumRows_)
    OffRankOffsets_[lclRowInd] = ptr_[lclRowInd];
  else {
    // TODO: use parallel reduce
    size_t offset = ptr_[lclRowInd+1];
    for (size_t j = ptr_[lclRowInd]; j < ptr_[lclRowInd+1]; j++) {
      const LO lclColInd = ind_[j];
      const GO gblColInd = lclColMap_.getGlobalElement (lclColInd);
      const LO lclDomInd = lclDomMap_.getLocalElement (gblColInd);
      if ((lclDomInd == INVALID) && (j < offset))
        offset = j;
    }
    OffRankOffsets_[lclRowInd] = offset;
  }
}

} // namespace Impl
} // namespace Details
} // namespace Tpetra

// Explicit template instantiation macro for
// Tpetra::Details::Impl::GetGraphOffRankOffsets.  NOT FOR USERS!!!  Must
// be used inside the Tpetra namespace.
#define TPETRA_DETAILS_IMPL_GETGRAPHOFFRANKOFFSETS_INSTANT( LO, GO, NODE ) \
  template class Details::Impl::GetGraphOffRankOffsets< LO, GO, NODE::device_type >;

#endif // TPETRA_DETAILS_GETGRAPHOFFRANKOFFSETS_DEF_HPP
