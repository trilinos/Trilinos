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
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LOCALROWOFFSETS_DEF_HPP
#define TPETRA_DETAILS_LOCALROWOFFSETS_DEF_HPP

/// \file Tpetra_Details_localRowOffsets_def.hpp
/// \brief Definition of function for getting local row offsets from a
///   Tpetra::RowGraph.

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {
namespace Impl {

template <class LO, class GO, class NT>
std::pair<typename LocalRowOffsetsResult<NT>::offsets_type, size_t>
localRowCounts (const RowGraph<LO, GO, NT>& G)
{
  using result_type = LocalRowOffsetsResult<NT>;
  using offsets_type = typename result_type::offsets_type;
  using offset_type = typename result_type::offset_type;

  const LO lclNumRows (G.getNodeNumRows ());
  offsets_type entPerRow;
  if (lclNumRows != 0) {
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    entPerRow =
      offsets_type (view_alloc ("entPerRow", WithoutInitializing),
                    lclNumRows);
  }
  using host = Kokkos::DefaultHostExecutionSpace;
  auto entPerRow_h = Kokkos::create_mirror_view (host (), entPerRow);

  // Don't trust G.getNodeMaxNumRowEntries() unless G is fillComplete.
  // Even then, I would rather this method didn't exist (since it adds
  // state and imposes overhead on fillComplete), and it's easy to
  // compute ourselves here.
  size_t maxNumEnt = 0;
  for (LO i = 0; i < lclNumRows; ++i) {
    const size_t lclNumEnt = G.getNumEntriesInLocalRow (i);
    entPerRow_h[i] = offset_type (lclNumEnt);
    maxNumEnt = maxNumEnt < lclNumEnt ? lclNumEnt : maxNumEnt;
  }
  Kokkos::deep_copy (entPerRow, entPerRow_h);
  return {entPerRow, maxNumEnt};
}

template <class LO, class GO, class NT>
LocalRowOffsetsResult<NT>
localRowOffsetsFromRowGraph (const RowGraph<LO, GO, NT>& G)
{
  using result_type = LocalRowOffsetsResult<NT>;
  using offsets_type = typename result_type::offsets_type;
  using offset_type = typename result_type::offset_type;

  offsets_type entPerRow;
  size_t maxNumEnt = 0;
  {
    auto result = localRowCounts (G);
    entPerRow = result.first;
    maxNumEnt = result.second;
  }

  const LO lclNumRows (G.getNodeNumRows ());
  offsets_type ptr;
  offset_type nnz = 0;
  if (lclNumRows != 0) {
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    ptr = offsets_type (view_alloc ("ptr", WithoutInitializing),
                        lclNumRows + 1);
    using ::Tpetra::Details::computeOffsetsFromCounts;
    nnz = computeOffsetsFromCounts (ptr, entPerRow);
  }
  return {ptr, nnz, maxNumEnt};
}

template <class LO, class GO, class NT>
LocalRowOffsetsResult<NT>
localRowOffsetsFromFillCompleteCrsGraph (const CrsGraph<LO, GO, NT>& G)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using result_type = LocalRowOffsetsResult<NT>;
  using offsets_type = typename result_type::offsets_type;
  using offset_type = typename result_type::offset_type;

  auto G_lcl = G.getLocalGraph ();
  offsets_type ptr (view_alloc ("ptr", WithoutInitializing),
                    G_lcl.row_map.extent (0));
  Kokkos::deep_copy (ptr, G_lcl.row_map);

  const offset_type nnz = G.getNodeNumEntries ();
  const size_t maxNumEnt = G.getNodeMaxNumRowEntries ();
  return {ptr, nnz, maxNumEnt};
}

} // namespace Impl

template <class LO, class GO, class NT>
LocalRowOffsetsResult<NT>
localRowOffsets (const RowGraph<LO, GO, NT>& G)
{
  if (G.isFillComplete ()) {
    using crs_graph_type = CrsGraph<LO, GO, NT>;
    const crs_graph_type* G_crs =
      dynamic_cast<const crs_graph_type*> (&G);
    if (G_crs != nullptr) {
      return Impl::localRowOffsetsFromFillCompleteCrsGraph (*G_crs);
    }
  }
  return Impl::localRowOffsetsFromRowGraph (G);
}

} // namespace Details
} // namespace Tpetra

//
// Explicit instantiation macros
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_DETAILS_LOCALROWOFFSETS_INSTANT(LO, GO, NT) \
namespace Details { \
namespace Impl { \
  \
template std::pair<LocalRowOffsetsResult<NT>::offsets_type, size_t> \
localRowCounts (const RowGraph<LO, GO, NT>& G); \
  \
template LocalRowOffsetsResult<NT> \
localRowOffsetsFromRowGraph (const RowGraph<LO, GO, NT>& G); \
  \
template LocalRowOffsetsResult<NT> \
localRowOffsetsFromFillCompleteCrsGraph (const CrsGraph<LO, GO, NT>& G); \
  \
} \
  \
template LocalRowOffsetsResult<NT> \
localRowOffsets (const RowGraph<LO, GO, NT>& A); \
}

#endif // TPETRA_DETAILS_LOCALROWOFFSETS_DEF_HPP
