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

#include "LOCA_MultiPredictor_Tangent.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"
#include "Teuchos_ParameterList.hpp"

LOCA::MultiPredictor::Tangent::Tangent(
	      const Teuchos::RCP<LOCA::GlobalData>& global_data,
	      const Teuchos::RCP<Teuchos::ParameterList>& predParams,
	      const Teuchos::RCP<Teuchos::ParameterList>& solverParams) :
  globalData(global_data),
  linSolverParams(solverParams),
  fdfdp(),
  tangent(),
  secant(),
  initialized(false)
{
}

LOCA::MultiPredictor::Tangent::~Tangent()
{
}

LOCA::MultiPredictor::Tangent::Tangent(
				 const LOCA::MultiPredictor::Tangent& source,
				 NOX::CopyType type) :
  globalData(source.globalData),
  linSolverParams(source.linSolverParams),
  fdfdp(),
  tangent(),
  secant(),
  initialized(source.initialized)
{
  if (source.initialized) {
    fdfdp = source.fdfdp->clone(type);
  
    tangent = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.tangent->clone(type));

    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(type));
  }
}

LOCA::MultiPredictor::AbstractStrategy&
LOCA::MultiPredictor::Tangent::operator=(
			  const LOCA::MultiPredictor::AbstractStrategy& s)
{
  const LOCA::MultiPredictor::Tangent& source = 
    dynamic_cast<const LOCA::MultiPredictor::Tangent&>(s);

  if (this != &source) {
    globalData = source.globalData;
    linSolverParams = source.linSolverParams;
    initialized = source.initialized;

    if (source.initialized) {
      fdfdp = source.fdfdp->clone(NOX::DeepCopy);
  
      tangent = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.tangent->clone(NOX::DeepCopy));

      secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(NOX::DeepCopy));
    }
  }

  return *this;
}

Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Tangent::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Tangent(*this, type));
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Tangent::compute(
	      bool baseOnSecant, const std::vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec)
{
  std::string callingFunction = "LOCA::MultiPredictor::Tangent::compute()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails))
    globalData->locaUtils->out() << 
      "\n\tCalling Predictor with method: Tangent" << std::endl;

  // Number of continuation parameters
  int numParams = stepSize.size();

  // Get underlying group
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGroup 
    = grp.getUnderlyingGroup();

  if (!initialized) {

    // Allocate dfdp
    fdfdp = underlyingGroup->getX().createMultiVector(numParams+1, 
						      NOX::ShapeCopy);

    // Allocate tangent
    tangent = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(xVec.createMultiVector(numParams, NOX::ShapeCopy));
    
    // Allocate secant
    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(xVec.clone(NOX::ShapeCopy));

    initialized = true;
  }

  // Get references to x, parameter components of predictor
  Teuchos::RCP<NOX::Abstract::MultiVector> tanX = 
    tangent->getXMultiVec();
  Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> tanP = 
    tangent->getScalars();

  // Get continuation parameter IDs
  const std::vector<int>& conParamIDs = grp.getContinuationParameterIDs();

  // Compute derivative of residual w.r.t. parameter
  finalStatus = underlyingGroup->computeDfDpMulti(conParamIDs, *fdfdp, false);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  std::vector<int> index_dfdp(conParamIDs.size());
  for (unsigned int i=0; i<conParamIDs.size(); i++)
    index_dfdp[i] = i+1;
  Teuchos::RCP<NOX::Abstract::MultiVector>dfdp = 
    fdfdp->subView(index_dfdp);

  // Scale dfdp by -1.0
  for (unsigned int i=0; i<conParamIDs.size(); i++)
    (*dfdp)[i].scale(-1.0);

  // Compute Jacobian
  status = underlyingGroup->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  // Solve J*tanX = -df/dp
  status = underlyingGroup->applyJacobianInverseMultiVector(*linSolverParams, 
							    *dfdp, 
							    *tanX);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Set parameter component equal to identity
  tanP->putScalar(0.0);
  for (unsigned int i=0; i<conParamIDs.size(); i++)
    (*tanP)(i,i) = 1.0;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXVec, 
			  xVec, *secant, *tangent);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Tangent::evaluate(
	      const std::vector<double>& stepSize,
	      const LOCA::MultiContinuation::ExtendedVector& xVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result) const
{
  // Number of continuation parameters
  int numParams = stepSize.size();

  for (int i=0; i<numParams; i++)
    result[i].update(1.0, xVec, stepSize[i], (*tangent)[i], 0.0);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Tangent::computeTangent(
			LOCA::MultiContinuation::ExtendedMultiVector& v)
{
  v = *tangent;

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiPredictor::Tangent::isTangentScalable() const
{
  return true;
}
