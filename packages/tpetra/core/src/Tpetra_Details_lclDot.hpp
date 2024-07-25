// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LCLDOT_HPP
#define TPETRA_DETAILS_LCLDOT_HPP

/// \file Tpetra_Details_lclDot.hpp
/// \brief Declaration and definition of Tpetra::Details::lclDot, an
///   implementation detail of Tpetra::MultiVector.
///
/// \warning This file, and its contents, are implementation details
///   of Tpetra, and may change or disappear at any time.

#include "Kokkos_DualView.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosBlas1_dot.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_TestForException.hpp"

namespace Tpetra {
namespace Details {

template<class RV, class XMV>
void
lclDot (const RV& dotsOut,
        const XMV& X_lcl,
        const XMV& Y_lcl,
        const size_t lclNumRows,
        const size_t numVecs,
        const size_t whichVecsX[],
        const size_t whichVecsY[],
        const bool constantStrideX,
        const bool constantStrideY)
{
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef typename RV::non_const_value_type dot_type;
#ifdef HAVE_TPETRA_DEBUG
  const char prefix[] = "Tpetra::MultiVector::lclDotImpl: ";
#endif // HAVE_TPETRA_DEBUG

  static_assert (Kokkos::is_view<RV>::value,
                 "Tpetra::MultiVector::lclDotImpl: "
                 "The first argument dotsOut is not a Kokkos::View.");
  static_assert (RV::rank == 1, "Tpetra::MultiVector::lclDotImpl: "
                 "The first argument dotsOut must have rank 1.");
  static_assert (Kokkos::is_view<XMV>::value,
                 "Tpetra::MultiVector::lclDotImpl: The type of the 2nd and "
                 "3rd arguments (X_lcl and Y_lcl) is not a Kokkos::View.");
  static_assert (XMV::rank == 2, "Tpetra::MultiVector::lclDotImpl: "
                 "X_lcl and Y_lcl must have rank 2.");

  // In case the input dimensions don't match, make sure that we
  // don't overwrite memory that doesn't belong to us, by using
  // subset views with the minimum dimensions over all input.
  const std::pair<size_t, size_t> rowRng (0, lclNumRows);
  const std::pair<size_t, size_t> colRng (0, numVecs);
  RV theDots = subview (dotsOut, colRng);
  XMV X = subview (X_lcl, rowRng, ALL ());
  XMV Y = subview (Y_lcl, rowRng, ALL ());

#ifdef HAVE_TPETRA_DEBUG
  if (lclNumRows != 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.extent (0) != lclNumRows, std::logic_error, prefix <<
       "X.extent(0) = " << X.extent (0) << " != lclNumRows "
       "= " << lclNumRows << ".  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (Y.extent (0) != lclNumRows, std::logic_error, prefix <<
       "Y.extent(0) = " << Y.extent (0) << " != lclNumRows "
       "= " << lclNumRows << ".  "
       "Please report this bug to the Tpetra developers.");
    // If a MultiVector is constant stride, then numVecs should
    // equal its View's number of columns.  Otherwise, numVecs
    // should be less than its View's number of columns.
    TEUCHOS_TEST_FOR_EXCEPTION
      (constantStrideX &&
       (X.extent (0) != lclNumRows || X.extent (1) != numVecs),
       std::logic_error, prefix << "X is " << X.extent (0) << " x " <<
       X.extent (1) << " (constant stride), which differs from the "
       "local dimensions " << lclNumRows << " x " << numVecs << ".  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! constantStrideX &&
       (X.extent (0) != lclNumRows || X.extent (1) < numVecs),
       std::logic_error, prefix << "X is " << X.extent (0) << " x " <<
       X.extent (1) << " (NOT constant stride), but the local "
       "dimensions are " << lclNumRows << " x " << numVecs << ".  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (constantStrideY &&
       (Y.extent (0) != lclNumRows || Y.extent (1) != numVecs),
       std::logic_error, prefix << "Y is " << Y.extent (0) << " x " <<
       Y.extent (1) << " (constant stride), which differs from the "
       "local dimensions " << lclNumRows << " x " << numVecs << ".  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (! constantStrideY &&
       (Y.extent (0) != lclNumRows || Y.extent (1) < numVecs),
       std::logic_error, prefix << "Y is " << Y.extent (0) << " x " <<
       Y.extent (1) << " (NOT constant stride), but the local "
       "dimensions are " << lclNumRows << " x " << numVecs << ".  "
       "Please report this bug to the Tpetra developers.");
  }
#endif // HAVE_TPETRA_DEBUG

  if (lclNumRows == 0) {
    const dot_type zero = Kokkos::ArithTraits<dot_type>::zero ();
    // DEEP_COPY REVIEW - NOT TESTED
    Kokkos::deep_copy (theDots, zero);
  }
  else { // lclNumRows != 0
    if (constantStrideX && constantStrideY) {
      if (X.extent (1) == 1) {
        typename RV::non_const_value_type result =
          KokkosBlas::dot (subview (X, ALL (), 0), subview (Y, ALL (), 0));
        // DEEP_COPY REVIEW - NOT TESTED
        Kokkos::deep_copy (theDots, result);
      }
      else {
        KokkosBlas::dot (theDots, X, Y);
      }
    }
    else { // not constant stride
      // NOTE (mfh 15 Jul 2014) This does a kernel launch for
      // every column.  It might be better to have a kernel that
      // does the work all at once.  On the other hand, we don't
      // prioritize performance of MultiVector views of
      // noncontiguous columns.
      for (size_t k = 0; k < numVecs; ++k) {
        const size_t X_col = constantStrideX ? k : whichVecsX[k];
        const size_t Y_col = constantStrideY ? k : whichVecsY[k];
        KokkosBlas::dot (subview (theDots, k), subview (X, ALL (), X_col),
                         subview (Y, ALL (), Y_col));
      } // for each column
    } // constantStride
  } // lclNumRows != 0
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LCLDOT_HPP
