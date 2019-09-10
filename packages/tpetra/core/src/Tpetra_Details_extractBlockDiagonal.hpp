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

#ifndef TPETRA_DETAILS_EXTRACTBLOCKDIAGONAL_HPP
#define TPETRA_DETAILS_EXTRACTBLOCKDIAGONAL_HPP

#include "TpetraCore_config.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Details_Behavior.hpp"

/// \file Tpetra_Details_extractBlockDiagonal.hpp
/// \brief Functions that allow for the extraction of a block diagonal
/// from a Tpetra::CrsMatrix.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {


template<class SparseMatrixType,
         class MultiVectorType>
void extractBlockDiagonal(const SparseMatrixType& A, MultiVectorType & diagonal) {
  using local_map_type = typename SparseMatrixType::map_type::local_map_type;
  using SC             = typename MultiVectorType::scalar_type;
  using LO             = typename SparseMatrixType::local_ordinal_type;
  using KCRS           = typename SparseMatrixType::local_matrix_type;
  using lno_view_t     = typename KCRS::StaticCrsGraphType::row_map_type::const_type;
  using lno_nnz_view_t = typename KCRS::StaticCrsGraphType::entries_type::const_type;
  using scalar_view_t  = typename KCRS::values_type::const_type;
  using local_mv_type  = typename MultiVectorType::dual_view_type::t_dev;
  using range_type     = Kokkos::RangePolicy<typename SparseMatrixType::node_type::execution_space, LO>;

  // Sanity checking: Map Compatibility (A's rowmap matches diagonal's map)
  if (Tpetra::Details::Behavior::debug() == true) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.getRowMap()->isSameAs(*diagonal.getMap()),
       std::runtime_error, "Tpetra::Details::extractBlockDiagonal was given incompatible maps");
  }

  LO numrows   = diagonal.getLocalLength();
  LO blocksize = diagonal.getNumVectors();
  SC ZERO = Teuchos::ScalarTraits<typename MultiVectorType::scalar_type>::zero();

  // Get Kokkos versions of objects
  local_map_type rowmap  = A.getRowMap()->getLocalMap();
  local_map_type colmap  = A.getRowMap()->getLocalMap();
  local_mv_type diag     = diagonal.getLocalViewDevice();
  const KCRS   Amat      = A.getLocalMatrix();
  lno_view_t Arowptr     = Amat.graph.row_map;
  lno_nnz_view_t Acolind = Amat.graph.entries;
  scalar_view_t Avals    = Amat.values;

  Kokkos::parallel_for("Tpetra::extractBlockDiagonal",range_type(0,numrows),KOKKOS_LAMBDA(const LO i){
      LO diag_col   = colmap.getLocalElement(rowmap.getGlobalElement(i));
      LO blockStart = diag_col - (diag_col % blocksize);
      LO blockStop  = blockStart + blocksize;
      for(LO k=0; k<blocksize; k++)
        diag(i,k)=ZERO;

      for (size_t k = Arowptr(i); k < Arowptr(i+1); k++) {
        LO col = Acolind(k);
        if (blockStart <= col && col < blockStop) {
          diag(i,col-blockStart) = Avals(k);
        }
      }
    });
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_EXTRACTBLOCKDIAGONAL_HPP
