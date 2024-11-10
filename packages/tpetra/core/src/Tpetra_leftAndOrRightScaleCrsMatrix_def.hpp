// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_leftScaleLocalCrsMatrix.hpp"
#include "Tpetra_Details_rightScaleLocalCrsMatrix.hpp"
#include "Teuchos_TestForException.hpp"

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
    auto A_lcl = A.getLocalMatrixDevice ();
    const LO lclNumRows =
      static_cast<LO> (A.getRowMap ()->getLocalNumElements ());
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
    Details::leftScaleLocalCrsMatrix (A.getLocalMatrixDevice (),
                                      rowScalingFactors,
                                      assumeSymmetric,
                                      divide);
  }
  if (rightScale) {
    Details::rightScaleLocalCrsMatrix (A.getLocalMatrixDevice (),
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

template<class SC, class LO, class GO, class NT>
void
leftAndOrRightScaleCrsMatrix (Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                              const Tpetra::Vector<
                                typename Kokkos::ArithTraits<SC>::mag_type,
                                LO, GO, NT>& rowScalingFactors,
                              const Tpetra::Vector<
                                typename Kokkos::ArithTraits<SC>::mag_type,
                                LO, GO, NT>& colScalingFactors,
                              const bool leftScale,
                              const bool rightScale,
                              const bool assumeSymmetric,
                              const EScaling scaling)
{
  using device_type = typename NT::device_type;
  using dev_memory_space = typename device_type::memory_space;
  using mag_type = typename Kokkos::ArithTraits<SC>::mag_type;
  const char prefix[] = "leftAndOrRightScaleCrsMatrix: ";
  const bool debug = ::Tpetra::Details::Behavior::debug ();

  Kokkos::View<const mag_type*, device_type> row_lcl;
  Kokkos::View<const mag_type*, device_type> col_lcl;
  if (leftScale) {
    if (debug) {
      const bool same = rowScalingFactors.getMap ()->isSameAs (* (A.getRowMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION
        (! same, std::invalid_argument, prefix << "rowScalingFactors's Map "
         "must be the same as the CrsMatrix's row Map.  If you see this "
         "message, it's likely that you are using a range Map Vector and that "
         "the CrsMatrix's row Map is overlapping.");
    }
    // if (rowScalingFactors.template need_sync<dev_memory_space> ()) {
    //   const_cast<vec_type&> (rowScalingFactors).template sync<dev_memory_space> ();
    // }
    auto row_lcl_2d = rowScalingFactors.template getLocalView<dev_memory_space> (Access::ReadOnly);
    row_lcl = Kokkos::subview (row_lcl_2d, Kokkos::ALL (), 0);
  }
  if (rightScale) {
    if (debug) {
      const bool same = colScalingFactors.getMap ()->isSameAs (* (A.getColMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION
        (! same, std::invalid_argument, prefix << "colScalingFactors's Map "
         "must be the same as the CrsMatrix's column Map.  If you see this "
         "message, it's likely that you are using a domain Map Vector.");
    }
    // if (colScalingFactors.template need_sync<dev_memory_space> ()) {
    //   const_cast<vec_type&> (colScalingFactors).template sync<dev_memory_space> ();
    // }
    auto col_lcl_2d = colScalingFactors.template getLocalView<dev_memory_space> (Access::ReadOnly);
    col_lcl = Kokkos::subview (col_lcl_2d, Kokkos::ALL (), 0);
  }

  leftAndOrRightScaleCrsMatrix (A, row_lcl, col_lcl, leftScale, rightScale,
                                assumeSymmetric, scaling);
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
    const EScaling scaling); \
  \
  template void \
  leftAndOrRightScaleCrsMatrix ( \
    Tpetra::CrsMatrix<SC, LO, GO, NT>& A, \
    const Tpetra::Vector<Kokkos::ArithTraits<SC>::mag_type, LO, GO, NT>& rowScalingFactors, \
    const Tpetra::Vector<Kokkos::ArithTraits<SC>::mag_type, LO, GO, NT>& colScalingFactors, \
    const bool leftScale, \
    const bool rightScale, \
    const bool assumeSymmetric, \
    const EScaling scaling);

#endif // TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DEF_HPP
