// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_MultiPredictor_AbstractStrategy.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"

void
LOCA::MultiPredictor::AbstractStrategy::setPredictorOrientation(
	      bool baseOnSecant, 
	      const vector<double>& stepSize,
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
