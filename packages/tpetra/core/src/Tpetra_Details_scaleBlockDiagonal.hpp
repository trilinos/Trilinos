// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_SCALEBLOCKDIAGONAL_HPP
#define TPETRA_DETAILS_SCALEBLOCKDIAGONAL_HPP

#include "TpetraCore_config.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"


/// \file Tpetra_Details_scaleBlockDiagonal.hpp
/// \brief Functions that rescales a multivector (in-place) by a inverse of a block-diagonal matrix 
/// (as implied by the given MultiVector) or its transpose.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {


template<class MultiVectorType>        
void inverseScaleBlockDiagonal(MultiVectorType & blockDiagonal, bool doTranspose, MultiVectorType & multiVectorToBeScaled) {
  using LO             = typename MultiVectorType::local_ordinal_type;
  using range_type     = Kokkos::RangePolicy<typename MultiVectorType::node_type::execution_space, LO>;
  using namespace KokkosBatched;
  typename MultiVectorType::impl_scalar_type SC_one = Teuchos::ScalarTraits<typename MultiVectorType::impl_scalar_type>::one();

  // Sanity checking: Map Compatibility (A's rowmap matches diagonal's map)
  if (Tpetra::Details::Behavior::debug() == true) {
    TEUCHOS_TEST_FOR_EXCEPTION(!blockDiagonal.getMap()->isSameAs(*multiVectorToBeScaled.getMap()),
       std::runtime_error, "Tpetra::Details::scaledBlockDiagonal was given incompatible maps");
  }

  LO numrows   = blockDiagonal.getLocalLength();
  LO blocksize = blockDiagonal.getNumVectors();
  LO numblocks = numrows / blocksize;

  // Get Kokkos versions of objects

  auto blockDiag = blockDiagonal.getLocalViewDevice(Access::OverwriteAll);
  auto toScale   = multiVectorToBeScaled.getLocalViewDevice(Access::ReadWrite);
  
  typedef Algo::Level3::Unblocked algo_type;
  Kokkos::parallel_for("scaleBlockDiagonal",range_type(0,numblocks),KOKKOS_LAMBDA(const LO i){
      Kokkos::pair<LO,LO> row_range(i*blocksize,(i+1)*blocksize);
      auto A = Kokkos::subview(blockDiag,row_range, Kokkos::ALL());
      auto B = Kokkos::subview(toScale,  row_range, Kokkos::ALL());

      // Factor
      SerialLU<algo_type>::invoke(A);

      if(doTranspose) {
	// Solve U^T
	SerialTrsm<Side::Left,Uplo::Upper,Trans::Transpose,Diag::NonUnit,algo_type>::invoke(SC_one,A,B);
	// Solver L^T
	SerialTrsm<Side::Left,Uplo::Lower,Trans::Transpose,Diag::Unit,algo_type>::invoke(SC_one,A,B);
      }
      else {
	// Solve L
	SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,algo_type>::invoke(SC_one,A,B);	
	// Solve U
	SerialTrsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,algo_type>::invoke(SC_one,A,B);
      }
    });  
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_SCALEBLOCKDIAGONAL_HPP
