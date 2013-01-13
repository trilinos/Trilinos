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

#include "LOCA_Homotopy_Group.H"
#include "LOCA_Homotopy_AbstractGroup.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"

LOCA::Homotopy::Group::Group(
	 Teuchos::ParameterList& locaSublist,
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Homotopy::AbstractGroup>& g,
	 double scalarRandom,
	 double scalarInitialGuess) :
  globalData(global_data),
  grpPtr(g),
  gVecPtr(g->getX().clone(NOX::ShapeCopy)),
  randomVecPtr(gVecPtr->clone(NOX::ShapeCopy)),
  newtonVecPtr(),
  gradVecPtr(),
  paramVec(grpPtr->getParams()),
  conParam(0.0),
  conParamID(-1),
  conParamLabel("Homotopy Continuation Parameter"),
  augmentJacForHomotopyNotImplemented(false)
{
  // construct a random vector for the problem 
  randomVecPtr->random();
  randomVecPtr->abs(*randomVecPtr);
  randomVecPtr->update(scalarInitialGuess, grpPtr->getX(), scalarRandom);

  // Set the isValid flags to false
  resetIsValidFlags();

  // Set the homotopy parameter as a parameter in the ParameterVector.
  // This will allow us to do an invasive homotopy since the app can now 
  // get access to the homotopy continuation parameter.
  paramVec.addParameter(conParamLabel, conParam);
  grpPtr->setParams(paramVec);

  // Set up paramID
  conParamID = paramVec.getIndex(conParamLabel);

  setStepperParameters(locaSublist);

}

LOCA::Homotopy::Group::Group(
	 Teuchos::ParameterList& locaSublist,
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Homotopy::AbstractGroup>& g,
	 const NOX::Abstract::Vector& randomVector) :
  globalData(global_data),
  grpPtr(g),
  gVecPtr(g->getX().clone(NOX::ShapeCopy)),
  randomVecPtr(gVecPtr->clone(NOX::ShapeCopy)),
  newtonVecPtr(),
  gradVecPtr(),
  paramVec(grpPtr->getParams()),
  conParam(0.0),
  conParamID(-1),
  conParamLabel("Homotopy Continuation Parameter"),
  augmentJacForHomotopyNotImplemented(false)
{
  // construct a random vector for the problem 
  *randomVecPtr = randomVector;

  // Set the isValid flags to false
  resetIsValidFlags();

  // Set the homotopy parameter as a parameter in the ParameterVector.
  // This will allow us to do an invasive homotopy since the app can now 
  // get access to the homotopy continuation parameter.
  paramVec.addParameter(conParamLabel, conParam);
  grpPtr->setParams(paramVec);

  // Set up paramID
  conParamID = paramVec.getIndex(conParamLabel);

  setStepperParameters(locaSublist);

}

LOCA::Homotopy::Group::Group(const LOCA::Homotopy::Group& source, 
			     NOX::CopyType type) : 
  globalData(source.globalData),
  grpPtr(Teuchos::rcp_dynamic_cast<LOCA::Homotopy::AbstractGroup>(source.grpPtr->clone(type))),
  gVecPtr(source.gVecPtr->clone(type)),
  randomVecPtr(source.randomVecPtr->clone(NOX::DeepCopy)), //deep copy to keep the set of equations consistent
  newtonVecPtr(),
  gradVecPtr(),
  paramVec(source.paramVec),
  conParam(source.conParam),
  conParamID(source.conParamID),
  conParamLabel(source.conParamLabel),
  augmentJacForHomotopyNotImplemented(source.augmentJacForHomotopyNotImplemented)
{
  if (source.newtonVecPtr != Teuchos::null)
    newtonVecPtr = source.newtonVecPtr->clone(type);

  if (source.gradVecPtr != Teuchos::null)
    newtonVecPtr = source.gradVecPtr->clone(type);

  if (type == NOX::DeepCopy) {
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;
  }
  else if (type == NOX::ShapeCopy) {
    resetIsValidFlags();
  }
  else {
    globalData->locaErrorCheck->throwError(
			       "LOCA::Homotopy::Group::Group(copy ctor)",
			       "CopyType is invalid!");
  }

}


LOCA::Homotopy::Group::~Group() 
{
}

NOX::Abstract::Group&
LOCA::Homotopy::Group::operator=(const NOX::Abstract::Group& source) 
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Homotopy::Group::clone(NOX::CopyType type) const 
{
  return Teuchos::rcp(new LOCA::Homotopy::Group(*this, type));
}

void
LOCA::Homotopy::Group::setX(const NOX::Abstract::Vector& y) 
{
  resetIsValidFlags();
  grpPtr->setX(y);
}

void
LOCA::Homotopy::Group::computeX(const NOX::Abstract::Group& grp, 
				const NOX::Abstract::Vector& d,
				double step) 
{
  resetIsValidFlags();
  const LOCA::Homotopy::Group& g = 
    dynamic_cast<const LOCA::Homotopy::Group&>(grp);
  grpPtr->computeX(*(g.grpPtr), d, step);
}

NOX::Abstract::Group::ReturnType 
LOCA::Homotopy::Group::computeF()
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  // compute f(x)
  NOX::Abstract::Group::ReturnType status = grpPtr->computeF();
  
  // g = conParam * f(x) + ((1.0 - conParam) * (x - randomVec))
  *gVecPtr = grpPtr->getX();
  gVecPtr->update(-1.0, *randomVecPtr, 1.0); 
  gVecPtr->scale( 1.0 - conParam);
  gVecPtr->update(conParam, grpPtr->getF(), 1.0);

  isValidF = true;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  // Compute J
  NOX::Abstract::Group::ReturnType status = grpPtr->computeJacobian();

  // Compute p*J + (1-p)*I
  NOX::Abstract::Group::ReturnType augHomTest = 
    grpPtr->augmentJacobianForHomotopy(conParam, 1.0-conParam);
  // If it is not implemented, augment the Jacobian during the 
  // applyJacobian() call.
  if (augHomTest == NOX::Abstract::Group::NotDefined)
    augmentJacForHomotopyNotImplemented = true;

  isValidJacobian = true;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeGradient() 
{
  if (isValidGradient)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Homotopy::Group::computeGradient()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  if (gradVecPtr == Teuchos::null)
    gradVecPtr = gVecPtr->clone(NOX::ShapeCopy);

  finalStatus = computeF();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);
  
  status = computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  status = applyJacobianTranspose(*gVecPtr, *gradVecPtr);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  return finalStatus;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeNewton(Teuchos::ParameterList& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  std::string callingFunction = 
    "LOCA::Homotopy::Group::computeNewton()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  
  if (newtonVecPtr == Teuchos::null)
    newtonVecPtr = gVecPtr->clone(NOX::ShapeCopy);
  
  finalStatus = computeF();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);
  
  status = computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  status = applyJacobianInverse(params, *gVecPtr, *newtonVecPtr);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  newtonVecPtr->scale(-1.0);
  
  isValidNewton = true;
  
  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobian(const NOX::Abstract::Vector& input,
				     NOX::Abstract::Vector& result) const 
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;
  
  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobian(input, result);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented)
    result.update(1.0-conParam, input, conParam);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianTranspose(
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;

  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianTranspose(input, result);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented)
    result.update(1.0-conParam, input, conParam);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianInverse(Teuchos::ParameterList& params,
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  return grpPtr->applyJacobianInverse(params, input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianMultiVector(
				      const NOX::Abstract::MultiVector& input,
				      NOX::Abstract::MultiVector& result) const
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;
  
  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianMultiVector(input, result);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented)
    result.update(1.0-conParam, input, conParam);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianTransposeMultiVector(
			              const NOX::Abstract::MultiVector& input,
			              NOX::Abstract::MultiVector& result) const
{
  if (!isValidJacobian)
    return NOX::Abstract::Group::BadDependency;

  NOX::Abstract::Group::ReturnType status = 
    grpPtr->applyJacobianTransposeMultiVector(input, result);

  // If the Jacobian is not augmented for homotopy (i.e. using MFNK)
  // then lets augment it here.
  if (augmentJacForHomotopyNotImplemented)
    result.update(1.0-conParam, input, conParam);

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianInverseMultiVector(
				      Teuchos::ParameterList& params,
				      const NOX::Abstract::MultiVector& input,
				      NOX::Abstract::MultiVector& result) const
{
  return grpPtr->applyJacobianInverseMultiVector(params, input, result);
}

bool
LOCA::Homotopy::Group::isF() const 
{
  return isValidF;
}

bool
LOCA::Homotopy::Group::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Homotopy::Group::isNewton() const 
{
  return isValidNewton;
}
  
bool
LOCA::Homotopy::Group::isGradient() const 
{
  return isValidGradient;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getX() const 
{
  return grpPtr->getX();
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getF() const 
{
  return *gVecPtr;
}

double
LOCA::Homotopy::Group::getNormF() const 
{
  return gVecPtr->norm();
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getGradient() const 
{
  if (gradVecPtr == Teuchos::null) {
    globalData->locaErrorCheck->throwError(
				      "LOCA::Homotopy::Group::getGradient", 
				      "gradVecPtr is NULL!");
  }
  return *gradVecPtr;
}

const NOX::Abstract::Vector&
LOCA::Homotopy::Group::getNewton() const 
{
  if (newtonVecPtr == Teuchos::null) {
    globalData->locaErrorCheck->throwError("LOCA::Homotopy::Group::getNewton", 
					   "newtonVecPtr is NULL!");
  }
  return *newtonVecPtr;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::Group::getXPtr() const 
{
  return grpPtr->getXPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::Group::getFPtr() const 
{
  return gVecPtr;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::Group::getGradientPtr() const 
{
  if (gradVecPtr == Teuchos::null) {
    globalData->locaErrorCheck->throwError(
				      "LOCA::Homotopy::Group::getGradientPtr", 
				      "gradVecPtr is NULL!");
  }
  return gradVecPtr;
}

Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::Group::getNewtonPtr() const 
{
  if (newtonVecPtr == Teuchos::null) {
    globalData->locaErrorCheck->throwError("LOCA::Homotopy::Group::getNewtonPtr", 
					   "newtonVecPtr is NULL!");
  }
  return newtonVecPtr;
}

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup> 
LOCA::Homotopy::Group::getUnderlyingGroup() const
{
  return grpPtr;
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> 
LOCA::Homotopy::Group::getUnderlyingGroup()
{
  return grpPtr;
}


void
LOCA::Homotopy::Group::copy(const NOX::Abstract::Group& src) 
{

  const LOCA::Homotopy::Group& source = 
    dynamic_cast<const LOCA::Homotopy::Group&>(src);

  // Protect against A = A
  if (this != &source) {

    globalData = source.globalData;
     grpPtr->copy(*(source.grpPtr));
    *gVecPtr = *(source.gVecPtr);
    *randomVecPtr = *(source.randomVecPtr);
    if (newtonVecPtr != Teuchos::null)
      *newtonVecPtr = *(source.newtonVecPtr);
    if (gradVecPtr != Teuchos::null)
      *gradVecPtr = *(source.gradVecPtr);
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
    isValidGradient = source.isValidGradient;
    paramVec = source.paramVec;
    conParam = source.conParam;
    conParamID = source.conParamID;
    augmentJacForHomotopyNotImplemented = 
      source.augmentJacForHomotopyNotImplemented;
  }

}

void
LOCA::Homotopy::Group::setParamsMulti(
			  const std::vector<int>& paramIDs, 
			  const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  resetIsValidFlags();
  grpPtr->setParamsMulti(paramIDs, vals);
  for (unsigned int i=0; i<paramIDs.size(); i++)
    if (paramIDs[i] == conParamID)
      conParam = vals(i,0);
}

void
LOCA::Homotopy::Group::setParams(const LOCA::ParameterVector& p) 
{
  resetIsValidFlags();
  grpPtr->setParams(p);
  conParam = p.getValue(conParamLabel);
}

void
LOCA::Homotopy::Group::setParam(int paramID, double val)
{
  resetIsValidFlags();
  grpPtr->setParam(paramID, val);
  if (paramID == conParamID)
    conParam = val;
}

void
LOCA::Homotopy::Group::setParam(std::string paramID, double val)
{
  resetIsValidFlags();
  grpPtr->setParam(paramID, val);
  if (paramID == conParamLabel)
    conParam = val;
}

const LOCA::ParameterVector&
LOCA::Homotopy::Group::getParams() const 
{
  return grpPtr->getParams();
}

double
LOCA::Homotopy::Group::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

double
LOCA::Homotopy::Group::getParam(std::string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeDfDpMulti(const std::vector<int>& paramIDs, 
					NOX::Abstract::MultiVector& dfdp, 
					bool isValidF)
{
  // g = conParam * f(x) + ((1.0 - conParam) * (x - randomVec))
  // dg/dp = f(x) - (x - randomVec) when p = conParam
  // dg/dp = conParam * df/dp when p != conParam

  // For simplicity, we always recompute g, even if isValidF is true
  // Things get kind of messy otherwise

  // Extract parameter IDs that are not the continuation parameter
  std::vector<int> pIDs;
  std::vector<int> idx(1);
  idx[0] = 0; // index 0 corrsponds to f in dfdp
  for (unsigned int i=0; i<paramIDs.size(); i++)
    if (paramIDs[i] != conParamID) {
      pIDs.push_back(paramIDs[i]);
      idx.push_back(i+1);
    }

  // Create view of dfdp for parameters that aren't the continuation parameter
  Teuchos::RCP<NOX::Abstract::MultiVector> fp = 
    dfdp.subView(idx);

  // Compute df/dp for non-continuation parameter parameters
  // We force recomputation of f for simplicity
  NOX::Abstract::Group::ReturnType status = 
    grpPtr->computeDfDpMulti(pIDs, *fp, false);

  // Compute conParam * df/dp
  fp->scale(conParam);

  // Compute g
  double v = 1.0-conParam;
  dfdp[0].update(v, grpPtr->getX(), -v, *randomVecPtr, 1.0);

  // Compute dg/dp for p = conParam
  grpPtr->computeF();
  for (unsigned int i=0; i<paramIDs.size(); i++)
    if (paramIDs[i] == conParamID) {
      dfdp[i+1] = grpPtr->getF();
      dfdp[i+1].update(-1.0, grpPtr->getX(), 1.0, *randomVecPtr, 1.0);
    }

  return status;
}

void
LOCA::Homotopy::Group::preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->preProcessContinuationStep(stepStatus);
}

void
LOCA::Homotopy::Group::postProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  grpPtr->postProcessContinuationStep(stepStatus);
}

void
LOCA::Homotopy::Group::projectToDraw(const NOX::Abstract::Vector& x,
				     double *px) const
{
  grpPtr->projectToDraw(x, px);
  px[this->projectToDrawDimension()] = conParam;
}

int
LOCA::Homotopy::Group::projectToDrawDimension() const
{
  return grpPtr->projectToDrawDimension()+1;
}

void
LOCA::Homotopy::Group::printSolution(const double conParm) const 
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for homotopy parameter = " << 
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(conParam);
  return;
}

void
LOCA::Homotopy::Group::printSolution(const NOX::Abstract::Vector& x_,
				     const double conParm) const 
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for homotopy parameter = " << 
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(conParam);
  return;
}

void
LOCA::Homotopy::Group::resetIsValidFlags()
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
  isValidGradient = false;
  return;
}

void
LOCA::Homotopy::Group::setStepperParameters(Teuchos::ParameterList& params)
{
  
  // Create the stepper sublist
  Teuchos::ParameterList& stepperList = params.sublist("Stepper");
  stepperList.set("Continuation Method", "Natural");
  stepperList.set("Continuation Parameter", conParamLabel);
  stepperList.set("Initial Value", 0.0);
  stepperList.set("Max Value", 1.0);
  stepperList.set("Min Value", -1.0);
  stepperList.set("Max Steps", 50);

  // Create predictor sublist
  Teuchos::ParameterList& predictorList = params.sublist("Predictor");
  predictorList.set("Method", "Constant");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = params.sublist("Step Size");
  //stepSizeList.set("Method", "Constant");
  stepSizeList.set("Method", "Adaptive");
  stepSizeList.set("Initial Step Size", 0.1);
  stepSizeList.set("Min Step Size", 1.0e-2);
  stepSizeList.set("Max Step Size", 1.0);
  stepSizeList.set("Aggressiveness", 0.5);
  return;
}
