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

#include "Teuchos_ParameterList.hpp"

#include "LOCA_MultiPredictor_Secant.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"

LOCA::MultiPredictor::Secant::Secant(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& predParams) :
  globalData(global_data),
  firstStepPredictor(),
  isFirstStep(true),
  isFirstStepComputed(false),
  predictor(),
  secant(),
  initialized(false)
{
  Teuchos::RCP<Teuchos::ParameterList> firstStepList = 
    topParams->getSublist("First Step Predictor");
  // change default method to constant to avoid infinite stack recursion
  firstStepList->get("Method", "Constant");
  firstStepPredictor = globalData->locaFactory->createPredictorStrategy(
							       topParams,
							       firstStepList);
}

LOCA::MultiPredictor::Secant::~Secant()
{
}

LOCA::MultiPredictor::Secant::Secant(
				 const LOCA::MultiPredictor::Secant& source,
				 NOX::CopyType type) :
  globalData(source.globalData),
  firstStepPredictor(source.firstStepPredictor->clone(type)),
  isFirstStep(source.isFirstStep),
  isFirstStepComputed(source.isFirstStepComputed),
  predictor(),
  secant(),
  initialized(source.initialized)
{
  if (source.initialized) {
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.predictor->clone(type));

    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(type));
  }
}

LOCA::MultiPredictor::AbstractStrategy&
LOCA::MultiPredictor::Secant::operator=(
			  const LOCA::MultiPredictor::AbstractStrategy& s)
{
  const LOCA::MultiPredictor::Secant& source = 
    dynamic_cast<const LOCA::MultiPredictor::Secant&>(s);

  if (this != &source) {
    globalData = source.globalData;
    firstStepPredictor = source.firstStepPredictor->clone(NOX::DeepCopy);
    isFirstStep = source.isFirstStep;
    isFirstStepComputed = source.isFirstStepComputed;
    initialized = source.initialized;

    if (source.initialized) {
      predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.predictor->clone(NOX::DeepCopy));

      secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(NOX::DeepCopy));
    }
  }

  return *this;
}

Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Secant::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Secant(*this, type));
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Secant::compute(
	      bool baseOnSecant, const std::vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails))
    globalData->locaUtils->out() << 
      "\n\tCalling Predictor with method: Secant" << std::endl;

  // Number of continuation parameters
  int numParams = stepSize.size();

  if (!initialized) {

    // Allocate predictor vector
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(xVec.createMultiVector(numParams, NOX::ShapeCopy));
    
    // Allocate secant
    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(xVec.clone(NOX::ShapeCopy));
    
    initialized = true;
  }

  // Use first step predictor if this is the first step
  if (isFirstStep && !isFirstStepComputed) {
    isFirstStepComputed = true;
    return firstStepPredictor->compute(baseOnSecant, stepSize, grp, 
				       prevXVec, xVec);
  }

  if (isFirstStep && isFirstStepComputed)
    isFirstStep = false;
   
  // Compute x - xold
  (*predictor)[0].update(1.0, xVec, -1.0, prevXVec, 0.0);
  
  for (int i=0; i<numParams; i++) {
    
    (*predictor)[i] = (*predictor)[0];

    // Rescale so parameter component = 1
    (*predictor)[i].scale(1.0/fabs(predictor->getScalar(i,i)));

    // Set off-diagonal elements to 0
    for (int j=0; j<numParams; j++)
      if (i != j)
	predictor->getScalar(i,j) = 0.0;
  }
  
  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXVec, 
			  xVec, *secant, *predictor);
  
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Secant::evaluate(
	      const std::vector<double>& stepSize,
	      const LOCA::MultiContinuation::ExtendedVector& xVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result) const
{
  // Use first-step predictor if this is the first step
  if (isFirstStep) {
    return firstStepPredictor->evaluate(stepSize, xVec, result);
  }

  // Number of continuation parameters
  int numParams = stepSize.size();

  for (int i=0; i<numParams; i++)
    result[i].update(1.0, xVec, stepSize[i], (*predictor)[i], 0.0);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Secant::computeTangent(
			LOCA::MultiContinuation::ExtendedMultiVector& v)
{
  // Use first-step predictor if this is the first step
  if (isFirstStep) {
    return firstStepPredictor->computeTangent(v);
  }

  v = *predictor;

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiPredictor::Secant::isTangentScalable() const
{
  // Use first-step predictor if this is the first step
  if (isFirstStep) {
    return firstStepPredictor->isTangentScalable();
  }

  return true;
}
