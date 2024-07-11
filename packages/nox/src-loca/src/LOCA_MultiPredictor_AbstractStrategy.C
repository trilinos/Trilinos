// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

void
LOCA::MultiPredictor::AbstractStrategy::setPredictorOrientation(
          bool baseOnSecant,
          const std::vector<double>& stepSize,
          const LOCA::MultiContinuation::ExtendedGroup& grp,
          const LOCA::MultiContinuation::ExtendedVector& prevXVec,
          const LOCA::MultiContinuation::ExtendedVector& xVec,
          LOCA::MultiContinuation::ExtendedVector& secant,
          LOCA::MultiContinuation::ExtendedMultiVector& tangent)
{
  // If orientation is not based on a secant vector (i.e., first or last
  // steps in a continuation run) make parameter component of predictor
  // positive
  if (!baseOnSecant) {
    for (int i=0; i<tangent.numVectors(); i++)
      if (tangent.getScalar(i,i) < 0.0)
    tangent[i].scale(-1.0);
    return;
  }

  // compute secant
  secant.update(1.0, xVec, -1.0, prevXVec, 0.0);

  for (int i=0; i<tangent.numVectors(); i++)
    if (grp.computeScaledDotProduct(secant, tangent[i])*stepSize[i] < 0.0)
      tangent[i].scale(-1.0);
}
