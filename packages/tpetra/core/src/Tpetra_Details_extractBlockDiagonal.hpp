// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  using KCRS           = typename SparseMatrixType::local_matrix_device_type;
  using lno_view_t     = typename KCRS::StaticCrsGraphType::row_map_type::const_type;
  using lno_nnz_view_t = typename KCRS::StaticCrsGraphType::entries_type::const_type;
  using scalar_view_t  = typename KCRS::values_type::const_type;
  using local_mv_type  = typename MultiVectorType::dual_view_type::t_dev;
  using range_type     = Kokkos::RangePolicy<typename SparseMatrixType::node_type::execution_space, LO>;
  using ATS        = Kokkos::ArithTraits<SC>;
  using impl_ATS = Kokkos::ArithTraits<typename ATS::val_type>;

  // Sanity checking: Map Compatibility (A's rowmap matches diagonal's map)
  if (Tpetra::Details::Behavior::debug() == true) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.getRowMap()->isSameAs(*diagonal.getMap()),
       std::runtime_error, "Tpetra::Details::extractBlockDiagonal was given incompatible maps");
  }

  LO numrows   = diagonal.getLocalLength();
  LO blocksize = diagonal.getNumVectors();

  // Get Kokkos versions of objects
  local_map_type rowmap  = A.getRowMap()->getLocalMap();
  local_map_type colmap  = A.getRowMap()->getLocalMap();
  local_mv_type diag     = diagonal.getLocalViewDevice(Access::OverwriteAll);
  const KCRS   Amat      = A.getLocalMatrixDevice();
  lno_view_t Arowptr     = Amat.graph.row_map;
  lno_nnz_view_t Acolind = Amat.graph.entries;
  scalar_view_t Avals    = Amat.values;

  Kokkos::parallel_for("Tpetra::extractBlockDiagonal",range_type(0,numrows),KOKKOS_LAMBDA(const LO i){
      LO diag_col   = colmap.getLocalElement(rowmap.getGlobalElement(i));
      LO blockStart = diag_col - (diag_col % blocksize);
      LO blockStop  = blockStart + blocksize;
      for(LO k=0; k<blocksize; k++)
        diag(i,k)=impl_ATS::zero();

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
