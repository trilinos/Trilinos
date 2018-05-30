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
#include "Tpetra_Details_leftScaleLocalCrsMatrix.hpp"
#include "Tpetra_Details_rightScaleLocalCrsMatrix.hpp"

namespace Tpetra {

template<class SC, class LO, class GO, class NT>
void
leftAndOrRightScaleCrsMatrix (Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                              const Kokkos::View<
                                const typename Kokkos::ArithTraits<SC>::mag_type*,
                                typename NT::device_type>& rowScalingFactors,
                              const Kokkos::View<
                                const typename Kokkos::ArithTraits<SC>::mag_type*,
                                typename NT::device_type>& colScalingFactors,
                              const bool leftScale,
                              const bool rightScale,
                              const bool assumeSymmetric,
                              const EScaling scaling)
{
  if (! leftScale && ! rightScale) {
    return;
  }

  const bool A_fillComplete_on_input = A.isFillComplete ();
  if (! A_fillComplete_on_input) {
    // Make sure that A has a valid local matrix.  It might not if it
    // was not created with a local matrix, and if fillComplete has
    // never been called on it before.  A never-initialized (and thus
    // invalid) local matrix has zero rows, because it was default
    // constructed.
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

  const bool divide = scaling == SCALING_DIVIDE;
  if (leftScale) {
    Details::leftScaleLocalCrsMatrix (A.getLocalMatrix (),
                                      rowScalingFactors,
                                      assumeSymmetric,
                                      divide);
  }
  if (rightScale) {
    Details::rightScaleLocalCrsMatrix (A.getLocalMatrix (),
                                       colScalingFactors,
                                       assumeSymmetric,
                                       divide);
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
    const Kokkos::View< \
      const Kokkos::ArithTraits<SC>::mag_type*, \
      NT::device_type>& rowScalingFactors, \
    const Kokkos::View< \
      const Kokkos::ArithTraits<SC>::mag_type*, \
      NT::device_type>& colScalingFactors, \
    const bool leftScale, \
    const bool rightScale, \
    const bool assumeSymmetric, \
    const EScaling scaling);

#endif // TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DEF_HPP
