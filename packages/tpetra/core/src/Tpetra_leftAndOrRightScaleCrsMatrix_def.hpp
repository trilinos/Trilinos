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

#ifndef TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DEF_HPP
#define TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DEF_HPP

/// \file Tpetra_leftAndOrRightScaleCrsMatrix_def.hpp
/// \brief Definition of Tpetra::leftAndOrRightScaleCrsMatrix
///
/// For the declaration of this function and its public Doxygen
/// documentation, please see
/// Tpetra_leftAndOrRightScaleCrsMatrix_decl.hpp in this directory.

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Details_EquilibrationInfo.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

/// \brief Kokkos::parallel_for functor that left-scales a
///   KokkosSparse::CrsMatrix.
///
/// \tparam LocalSparseMatrixType KokkosSparse::CrsMatrix specialization.
template<class LocalSparseMatrixType>
class LeftScaleLocalCrsMatrix {
public:
  using local_matrix_type = LocalSparseMatrixType;
  using val_type = typename local_matrix_type::value_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using device_type = typename local_matrix_type::device_type;

  LeftScaleLocalCrsMatrix (const local_matrix_type& A_lcl,
                           const EquilibrationInfo<val_type, device_type>& result) :
    A_lcl_ (A_lcl),
    rowNorms_ (result.rowNorms),
    assumeSymmetric_ (result.assumeSymmetric)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const typename local_matrix_type::ordinal_type lclRow) const
  {
    using LO = typename local_matrix_type::ordinal_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    const mag_type curRowNorm = rowNorms_(lclRow);
    // Dividing by zero gives unpleasant results, so don't do it.
    if (curRowNorm > KAM::zero ()) {
      const mag_type scalingFactor = assumeSymmetric_ ?
        KAM::sqrt (curRowNorm) : curRowNorm;
      auto curRow = A_lcl_.row (lclRow);
      const LO numEnt = curRow.length;
      for (LO k = 0; k < numEnt; ++k) {
        curRow.value (k) = curRow.value(k) / scalingFactor;
      }
    }
  }

  /// \brief leftAndOrRightScaleCrsMatrix should call this.
  ///
  /// Create a functor instance and invoke it in a Kokkos::parallel_for.
  static void
  run (const local_matrix_type& A_lcl,
       const EquilibrationInfo<val_type, device_type>& result)
  {
    using execution_space = typename device_type::execution_space;
    using LO = typename local_matrix_type::ordinal_type;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;

    LeftScaleLocalCrsMatrix<local_matrix_type> functor (A_lcl, result);
    const LO lclNumRows = A_lcl.numRows ();
    Kokkos::parallel_for ("leftScaleLocalCrsMatrix",
                          range_type (0, lclNumRows), functor);
  }

private:
  local_matrix_type A_lcl_;
  Kokkos::View<const mag_type*, device_type> rowNorms_;
  bool assumeSymmetric_;
};

/// \brief Kokkos::parallel_for functor that right-scales a
///   KokkosSparse::CrsMatrix.
///
/// \tparam LocalSparseMatrixType KokkosSparse::CrsMatrix specialization.
template<class LocalSparseMatrixType>
class RightScaleLocalCrsMatrix {
public:
  using local_matrix_type = LocalSparseMatrixType;
  using val_type = typename local_matrix_type::value_type;
  using device_type = typename local_matrix_type::device_type;

  RightScaleLocalCrsMatrix (const local_matrix_type& A_lcl,
                            const EquilibrationInfo<val_type, device_type>& result) :
    A_lcl_ (A_lcl),
    scalingFactors_ (result.assumeSymmetric ? result.colNorms : result.rowScaledColNorms),
    assumeSymmetric_ (result.assumeSymmetric)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const typename local_matrix_type::ordinal_type lclRow) const
  {
    using LO = typename local_matrix_type::ordinal_type;
    using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
    using KAM = Kokkos::ArithTraits<mag_type>;

    auto curRow = A_lcl_.row (lclRow);
    const LO numEnt = curRow.length;
    for (LO k = 0; k < numEnt; ++k) {
      const LO lclColInd = curRow.colidx(k);
      const mag_type curColNorm = scalingFactors_(lclColInd);
      if (curColNorm != KAM::zero ()) {
        const mag_type scalingFactor = assumeSymmetric_ ?
          KAM::sqrt (curColNorm) : curColNorm;
        curRow.value(k) = curRow.value(k) / scalingFactor;
      }
    }
  }

  /// \brief leftAndOrRightScaleCrsMatrix should call this.
  ///
  /// Create a functor instance and invoke it in a Kokkos::parallel_for.
  static void
  run (const local_matrix_type& A_lcl,
       const EquilibrationInfo<val_type, device_type>& result)
  {
    using execution_space = typename device_type::execution_space;
    using LO = typename local_matrix_type::ordinal_type;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;

    RightScaleLocalCrsMatrix<local_matrix_type> functor (A_lcl, result);
    Kokkos::parallel_for ("rightScaleLocalCrsMatrix",
                          range_type (0, A_lcl.numRows ()), functor);
  }

private:
  local_matrix_type A_lcl_;

  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  Kokkos::View<const mag_type*, device_type> scalingFactors_;
  bool assumeSymmetric_;
};

} // namespace Details

template<class SC, class LO, class GO, class NT>
void
leftAndOrRightScaleCrsMatrix (Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                              const Details::EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type, typename NT::device_type>& equib,
                              const bool leftScale,
                              const bool rightScale)
{
  if (! leftScale && ! rightScale) {
    return;
  }

  const bool A_fillComplete_on_input = A.isFillComplete ();
  if (! A_fillComplete_on_input) {
    // Make sure that A has a valid local matrix.  It might not if it
    // was not created with a local matrix, and if fillComplete has
    // never been called on it before.  It's invalid if the local
    // number of rows is zero, but shouldn't be.
    auto A_lcl = A.getLocalMatrix ();
    const LO lclNumRows =
      static_cast<LO> (A.getRowMap ()->getNodeNumElements ());
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_lcl.numRows () != lclNumRows, std::invalid_argument,
       "leftAndOrRightScaleCrsMatrix: Local matrix is not valid.  "
       "This means that A was not created with a local matrix, "
       "and that fillComplete has never yet been called on A before.  "
       "Please call fillComplete on A at least once first "
       "before calling this method.");
  }
  else {
    A.resumeFill ();
  }

  if (leftScale) {
    using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using impl_type = Details::LeftScaleLocalCrsMatrix<local_matrix_type>;
    impl_type::run (A.getLocalMatrix (), equib);
  }
  if (rightScale) {
    using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using impl_type = Details::RightScaleLocalCrsMatrix<local_matrix_type>;
    impl_type::run (A.getLocalMatrix (), equib);
  }

  if (A_fillComplete_on_input) { // put A back how we found it
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList ();
    params->set ("No Nonlocal Changes", true);
    A.fillComplete (A.getDomainMap (), A.getRangeMap (), params);
  }
}

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_INSTANT(SC,LO,GO,NT) \
  template void \
  leftAndOrRightScaleCrsMatrix ( \
    Tpetra::CrsMatrix<SC, LO, GO, NT>& A, \
    const Details::EquilibrationInfo<Kokkos::ArithTraits<SC>::val_type, NT::device_type>& equib, \
    const bool leftScale, \
    const bool rightScale);

#endif // TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DEF_HPP
