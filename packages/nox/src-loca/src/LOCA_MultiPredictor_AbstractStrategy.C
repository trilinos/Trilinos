// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
