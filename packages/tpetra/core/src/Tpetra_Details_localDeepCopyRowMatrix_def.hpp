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

#ifndef TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DEF_HPP
#define TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DEF_HPP

/// \file Tpetra_Details_localDeepCopyRowMatrix_def.hpp
/// \brief Definition of function for making a deep copy of a
///   Tpetra::RowMatrix's local matrix.

#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Details_localRowOffsets.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Kokkos_Core.hpp"
#include <algorithm>
#include <stdexcept>

namespace Tpetra {
namespace Details {

template <class SC, class LO, class GO, class NT>
KokkosSparse::CrsMatrix<
  typename Kokkos::ArithTraits<SC>::val_type,
    LO,
    typename NT::execution_space,
    void,
    size_t>
localDeepCopyLocallyIndexedRowMatrix
(const RowMatrix<SC, LO, GO, NT>& A,
 const char label[])
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A.isGloballyIndexed (), std::logic_error,
     "Tpetra::Details::localDeepCopyLocallyIndexedRowMatrix: "
     "Input RowMatrix is globally indexed.");

  using lro_result_type = LocalRowOffsetsResult<NT>;
  using offsets_type = typename lro_result_type::offsets_type;
  using offset_type = typename lro_result_type::offset_type;
  offsets_type ptr;
  offset_type nnz = 0;
  size_t maxNumEnt = 0;
  {
    const auto& G = * (A.getGraph ());
    const lro_result_type result = localRowOffsets (G);
    ptr = result.ptr;
    nnz = result.nnz;
    maxNumEnt = result.maxNumEnt;
  }

  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using IST = typename Kokkos::ArithTraits<SC>::val_type;
  using local_matrix_type = KokkosSparse::CrsMatrix<
    IST, LO, typename NT::device_type, void, size_t>;
  using local_graph_type =
    typename local_matrix_type::staticcrsgraph_type;
  using inds_type = typename local_graph_type::entries_type;
  inds_type ind (view_alloc ("ind", WithoutInitializing), nnz);
  auto ind_h = Kokkos::create_mirror_view (ind);

  using values_type = typename local_matrix_type::values_type;
  values_type val (view_alloc ("val", WithoutInitializing), nnz);
  auto val_h = Kokkos::create_mirror_view (val);

  const bool hasViews = A.supportsRowViews ();

  Teuchos::Array<LO> inputIndsBuf;
  Teuchos::Array<SC> inputValsBuf;
  if (! hasViews) {
    inputIndsBuf.resize (maxNumEnt);
    inputValsBuf.resize (maxNumEnt);
  }

  const LO lclNumRows (A.getNodeNumRows ());
  offset_type curPos = 0;
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    Teuchos::ArrayView<const LO> inputInds_av;
    Teuchos::ArrayView<const SC> inputVals_av;
    size_t numEnt = 0;
    if (hasViews) {
      A.getLocalRowView (lclRow, inputInds_av,
                         inputVals_av);
      numEnt = static_cast<size_t> (inputInds_av.size ());
    }
    else {
      A.getLocalRowCopy (lclRow, inputIndsBuf (),
                         inputValsBuf (), numEnt);
      inputInds_av = inputIndsBuf.view (0, numEnt);
      inputVals_av = inputValsBuf.view (0, numEnt);
    }
    const IST* inVals =
      reinterpret_cast<const IST*> (inputVals_av.getRawPtr ());
    const LO* inInds = inputInds_av.getRawPtr ();
    std::copy (inInds, inInds + numEnt, ind_h.data () + curPos);
    std::copy (inVals, inVals + numEnt, val_h.data () + curPos);
    curPos += offset_type (numEnt);
  }
  Kokkos::deep_copy (ind, ind_h);
  Kokkos::deep_copy (val, val_h);

  local_graph_type lclGraph (ind, ptr);
  const size_t numCols = A.getColMap ()->getNodeNumElements ();
  return local_matrix_type (label, numCols, val, lclGraph);
}

} // namespace Details
} // namespace Tpetra

//
// Explicit instantiation macros
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_INSTANT(SC, LO, GO, NT) \
namespace Details { \
  template KokkosSparse::CrsMatrix< \
    Kokkos::ArithTraits<SC>::val_type, \
    LO, NT::execution_space, void, size_t> \
  localDeepCopyLocallyIndexedRowMatrix<SC, LO, GO, NT> \
    (const RowMatrix<SC, LO, GO, NT>& A, \
     const char label[]); \
}

#endif // TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DEF_HPP
