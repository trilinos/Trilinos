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

#ifndef TPETRA_DETAILS_LOCALROWOFFSETS_DECL_HPP
#define TPETRA_DETAILS_LOCALROWOFFSETS_DECL_HPP

/// \file Tpetra_Details_localRowOffsets_decl.hpp
/// \brief Declaration of function for getting local row offsets from
///   a Tpetra::RowGraph.

#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_RowGraph_fwd.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include <utility> // pair

namespace Tpetra {
namespace Details {

//! Result returned by localRowOffsets (see below).
template <class NT>
struct LocalRowOffsetsResult {
private:
  using local_graph_type =
    typename KokkosSparse::CrsMatrix<
      double, int, typename NT::device_type, void, size_t>::
        staticcrsgraph_type;
public:
  using offsets_type =
    typename local_graph_type::row_map_type::non_const_type;
  using offset_type = typename offsets_type::non_const_value_type;

  offsets_type ptr; //!< Local row offsets (Kokkos::View)
  offset_type nnz;  //!< Local number of graph / matrix entries
  size_t maxNumEnt; //!< Max number of entries over all local rows
};

namespace Impl {

template <class LO, class GO, class NT>
std::pair<typename LocalRowOffsetsResult<NT>::offsets_type, size_t>
localRowCounts (const RowGraph<LO, GO, NT>& G);

template <class LO, class GO, class NT>
LocalRowOffsetsResult<NT>
localRowOffsetsFromRowGraph (const RowGraph<LO, GO, NT>& G);

template <class LO, class GO, class NT>
LocalRowOffsetsResult<NT>
localRowOffsetsFromFillCompleteCrsGraph (const CrsGraph<LO, GO, NT>& G);

} // namespace Impl

/// \brief Get local row offsets ("ptr", in compressed sparse row
///   terms) for the given graph.
template <class LO, class GO, class NT>
LocalRowOffsetsResult<NT>
localRowOffsets (const RowGraph<LO, GO, NT>& G);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LOCALROWOFFSETS_DECL_HPP
