// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_Constraint.H"
#include "LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
ExtendedGroup(
      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& tpParams,
      const Teuchos::RefCountPtr<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& grp)
  : LOCA::Extended::MultiAbstractGroup(),
    LOCA::MultiContinuation::AbstractGroup(),
    globalData(global_data),
    parsedParams(topParams),
    turningPointParams(tpParams),
    grpPtr(grp),
    constraint(),
    conGroup(),
    bifParamID(0)
{
  const char *func = "LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup()";

  // Get bifurcation parameter name
  if (!turningPointParams->isParameter("Bifurcation Parameter")) {
    globalData->locaErrorCheck->throwError(func,
				 "\"Bifurcation Parameter\" name is not set!");
  }
  string bifParamName = turningPointParams->get(
						  "Bifurcation Parameter",
						  "None");
  const ParameterVector& p = grpPtr->getParams();
  bifParamID = p.getIndex(bifParamName);

  // Get symmetric flag
  bool isSymmetric = 
    turningPointParams->get("Symmetric Jacobian", false);

  // Compute/get initial "a" & "b" vectors
  Teuchos::RefCountPtr<NOX::Abstract::Vector> aVecPtr;
  Teuchos::RefCountPtr<NOX::Abstract::Vector> bVecPtr;
  getInitialVectors(aVecPtr, bVecPtr, isSymmetric);

  // Create constraint equation
  string constraintMethod = turningPointParams->get(
						  "Constraint Method",
						  "Default");
  if (constraintMethod == "Default")
    constraint = 
      Teuchos::rcp(new LOCA::TurningPoint::MinimallyAugmented::Constraint(
							       globalData,
							       parsedParams,
							       tpParams,
							       grpPtr,
							       isSymmetric,
							       *aVecPtr,
							       bVecPtr.get(),
							       bifParamID));
  else if (constraintMethod == "Modified")
    constraint = 
      Teuchos::rcp(new LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint(
							       globalData,
							       parsedParams,
							       tpParams,
							       grpPtr,
							       isSymmetric,
							       *aVecPtr,
							       bVecPtr.get(),
							       bifParamID));
  else 
    globalData->locaErrorCheck->throwError(
		    func,
		    string("Unknown constraint method:  ") + constraintMethod);

  // Create constrained group
  std::vector<int> bifParamIDs(1);
  bifParamIDs[0] = bifParamID;
  conGroup = 
    Teuchos::rcp(new LOCA::MultiContinuation::ConstrainedGroup(
							    globalData, 
							    parsedParams,
							    turningPointParams,
							    grpPtr, 
							    constraint,
							    bifParamIDs));
}

LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
ExtendedGroup(
	  const LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup& source,
	  NOX::CopyType type)
  : globalData(source.globalData),
    parsedParams(source.parsedParams),
    turningPointParams(source.turningPointParams),
    grpPtr(),
    constraint(),
    conGroup(),
    bifParamID(source.bifParamID)
{
  conGroup = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ConstrainedGroup>(source.conGroup->clone(type));
  grpPtr = Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>(conGroup->getGroup());
  constraint = Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::Constraint>(conGroup->getConstraints());
  constraint->setGroup(grpPtr);
}


LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
~ExtendedGroup() 
{
}

NOX::Abstract::Group&
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
operator=(const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::Group>
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new ExtendedGroup(*this, type));
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
setX(const NOX::Abstract::Vector& y)  
{
  conGroup->setX(y);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
computeX(const NOX::Abstract::Group& g, 
	 const NOX::Abstract::Vector& d,
	 double step) 
{
  const LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup& mg = 
    dynamic_cast<const LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup&>(g);

  // set newton update in constraint
  Teuchos::RefCountPtr<LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint> mod_constraint = 
    Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint>(constraint);
  if (mod_constraint != Teuchos::null) {
    const LOCA::MultiContinuation::ExtendedVector& emv_d = 
      dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(d);
    mod_constraint->setNewtonUpdates(*(emv_d.getXVec()), emv_d.getScalar(0), 
				     step);
  }

  conGroup->computeX(*(mg.conGroup), d, step);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
computeF() 
{
  return conGroup->computeF();
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
computeJacobian() 
{
  return conGroup->computeJacobian();
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
computeGradient() 
{
  return conGroup->computeGradient();
}
   
NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
computeNewton(Teuchos::ParameterList& params) 
{
  return conGroup->computeNewton(params);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
applyJacobian(const NOX::Abstract::Vector& input,
	      NOX::Abstract::Vector& result) const 
{
  return conGroup->applyJacobian(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
applyJacobianTranspose(const NOX::Abstract::Vector& input,
		       NOX::Abstract::Vector& result) const 
{
  return conGroup->applyJacobianTranspose(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
applyJacobianInverse(Teuchos::ParameterList& params, 
		     const NOX::Abstract::Vector& input,
		     NOX::Abstract::Vector& result) const 
{
  return conGroup->applyJacobianInverse(params, input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
applyJacobianMultiVector(const NOX::Abstract::MultiVector& input,
			 NOX::Abstract::MultiVector& result) const 
{
  return conGroup->applyJacobianMultiVector(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
applyJacobianTransposeMultiVector(const NOX::Abstract::MultiVector& input,
				  NOX::Abstract::MultiVector& result) const 
{
  return conGroup->applyJacobianTransposeMultiVector(input, result);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
applyJacobianInverseMultiVector(Teuchos::ParameterList& params,
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const 
{
  return conGroup->applyJacobianInverseMultiVector(params, input, result);
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isF() const 
{
  return conGroup->isF();
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isJacobian() const 
{
  return conGroup->isJacobian();
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isGradient() const 
{
  return conGroup->isGradient();
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isNewton() const 
{
  return conGroup->isNewton();
}
  
const NOX::Abstract::Vector&
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getX() const 
{
  return conGroup->getX();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getF() const 
{
  return conGroup->getF();
}

double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getNormF() const 
{
  return conGroup->getNormF();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getGradient() const 
{
  return conGroup->getGradient();
}

const NOX::Abstract::Vector&
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getNewton() const 
{
  return conGroup->getNewton();
}

double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getNormNewtonSolveResidual() const 
{
  return conGroup->getNormNewtonSolveResidual();
}

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getUnderlyingGroup() const
{
  return conGroup->getUnderlyingGroup();
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getUnderlyingGroup()
{
  return conGroup->getUnderlyingGroup();
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
copy(const NOX::Abstract::Group& src) 
{

  const LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup& source = 
    dynamic_cast<const LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup&>(src);

  // Protect against A = A
  if (this != &source) {
    globalData = source.globalData;
    parsedParams = source.parsedParams;
    turningPointParams = source.turningPointParams;
    conGroup->copy(*source.conGroup);
    grpPtr = Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>(conGroup->getGroup());
    constraint = Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::Constraint>(conGroup->getConstraints());
    constraint->setGroup(grpPtr);
    bifParamID = source.bifParamID;
  }
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
setParamsMulti(const vector<int>& paramIDs, 
	       const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  conGroup->setParamsMulti(paramIDs, vals);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
setParams(const LOCA::ParameterVector& p)
{
  conGroup->setParams(p);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
setParam(int paramID, double val)
{
  conGroup->setParam(paramID, val);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
setParam(string paramID, double val)
{
  conGroup->setParam(paramID, val);
}

const LOCA::ParameterVector&
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getParams() const
{
  return conGroup->getParams();
}

double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getParam(int paramID) const
{
  return conGroup->getParam(paramID);
}

double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getParam(string paramID) const
{
  return conGroup->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
computeDfDpMulti(const vector<int>& paramIDs, 
		 NOX::Abstract::MultiVector& dfdp, 
		 bool isValidF)
{
  return conGroup->computeDfDpMulti(paramIDs, dfdp, isValidF);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  conGroup->preProcessContinuationStep(stepStatus);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  conGroup->postProcessContinuationStep(stepStatus);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
projectToDraw(const NOX::Abstract::Vector& x,
	      double *px) const
{
  conGroup->projectToDraw(x, px);
}

int
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
projectToDrawDimension() const
{
  return conGroup->projectToDrawDimension();
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
printSolution(const double conParam) const
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Turning Point located at: " << 
      globalData->locaUtils->sciformat(conParam) << "   " << 
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;

    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for conParam = " << 
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Right Null Vector for bif param = " << 
      globalData->locaUtils->sciformat(getBifParam()) << std::endl;
  }
  grpPtr->printSolution(*(constraint->getRightNullVec()), getBifParam());
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Left Null Vector for sigma = " << 
      globalData->locaUtils->sciformat(constraint->getSigma()) << std::endl;
  }
  grpPtr->printSolution(*(constraint->getLeftNullVec()), 
			constraint->getSigma());
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
printSolution(const NOX::Abstract::Vector& x,
	      const double conParam) const
{
  const LOCA::MultiContinuation::ExtendedVector& tp_x = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(x);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::printSolution\n";

    globalData->locaUtils->out() << "Turning Point located at: " << 
      globalData->locaUtils->sciformat(conParam) << "   " << 
      globalData->locaUtils->sciformat(tp_x.getScalar(0)) << std::endl;

    globalData->locaUtils->out() << 
      "\tPrinting Solution Vector for conParam = " << 
      globalData->locaUtils->sciformat(conParam) << std::endl;
  }
  grpPtr->printSolution(*(tp_x.getXVec()), conParam);
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Right Null Vector for bif param = " << 
      globalData->locaUtils->sciformat(tp_x.getScalar(0)) << std::endl;
  }
  grpPtr->printSolution(*(constraint->getRightNullVec()), tp_x.getScalar(0));
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
    globalData->locaUtils->out() << 
      "\tPrinting Left Null Vector for sigma = " << 
      globalData->locaUtils->sciformat(constraint->getSigma()) << std::endl;
  }
  grpPtr->printSolution(*(constraint->getLeftNullVec()), 
			constraint->getSigma());
}

int
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getBorderedWidth() const
{
  return conGroup->getBorderedWidth();
}

Teuchos::RefCountPtr<const NOX::Abstract::Group>
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getUnborderedGroup() const
{
  return conGroup->getUnborderedGroup();
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isCombinedAZero() const
{
  return conGroup->isCombinedAZero();
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isCombinedBZero() const
{
  return conGroup->isCombinedBZero();
}

bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
isCombinedCZero() const
{
  return conGroup->isCombinedCZero();
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
extractSolutionComponent(const NOX::Abstract::MultiVector& v,
			 NOX::Abstract::MultiVector& v_x) const
{
  conGroup->extractSolutionComponent(v, v_x);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
extractParameterComponent(bool use_transpose,
			  const NOX::Abstract::MultiVector& v,
			  NOX::Abstract::MultiVector::DenseMatrix& v_p) const
{
  conGroup->extractParameterComponent(use_transpose, v, v_p);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
loadNestedComponents(const NOX::Abstract::MultiVector& v_x,
		     const NOX::Abstract::MultiVector::DenseMatrix& v_p,
		     NOX::Abstract::MultiVector& v) const
{
  conGroup->loadNestedComponents(v_x, v_p, v);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
fillA(NOX::Abstract::MultiVector& A) const
{
  conGroup->fillA(A);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
fillB(NOX::Abstract::MultiVector& B) const
{
  conGroup->fillB(B);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
fillC(NOX::Abstract::MultiVector::DenseMatrix& C) const
{
  conGroup->fillC(C);
}

double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getBifParam() const
{
  return grpPtr->getParam(bifParamID);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
setBifParam(double param)
{
  conGroup->setParam(bifParamID, param);
}

void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::
getInitialVectors(Teuchos::RefCountPtr<NOX::Abstract::Vector>& aVecPtr,
		  Teuchos::RefCountPtr<NOX::Abstract::Vector>& bVecPtr,
		  bool isSymmetric)
{
  string callingFunction = 
    "LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getIntitialVectors()";

  // Get method
  string method = 
    turningPointParams->get("Initial Null Vector Computation",
				     "User Provided");
  if (method == "Solve df/dp") {
    NOX::Abstract::Group::ReturnType status;
    NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
    std::vector<int> paramID(1);
    paramID[0] = bifParamID;
    Teuchos::RefCountPtr<NOX::Abstract::MultiVector> fdfdp = 
      grpPtr->getX().createMultiVector(2);
    aVecPtr = grpPtr->getX().clone(NOX::ShapeCopy);
    bVecPtr = grpPtr->getX().clone(NOX::ShapeCopy);
    aVecPtr->init(0.0);
    bVecPtr->init(0.0);

    // Compute df/dp
    status = grpPtr->computeDfDpMulti(paramID, *fdfdp, false);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    // Compute J
    status = grpPtr->computeJacobian();
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    // Compute b = J^-1*dfdp
    Teuchos::RefCountPtr<Teuchos::ParameterList> lsParams =
      parsedParams->getSublist("Linear Solver");
    status = grpPtr->applyJacobianInverse(*lsParams, (*fdfdp)[1], *bVecPtr);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    // Compute a = J^-T*dfdp if necessary
    if (!isSymmetric) {
      // Cast group to one that can solve J^T
      Teuchos::RefCountPtr<LOCA::Abstract::TransposeSolveGroup> ts_grp = 
	Teuchos::rcp_dynamic_cast<LOCA::Abstract::TransposeSolveGroup>(grpPtr);
      if (ts_grp == Teuchos::null)
	globalData->locaErrorCheck->throwError(
	   callingFunction,
	   string("Group must implement LOCA::Abstract::TransposeSolveGroup") +
	   string(" to compute initial left null vector"));
      
      Teuchos::RefCountPtr<Teuchos::ParameterList> lsParams =
	parsedParams->getSublist("Linear Solver");
      status = 
	ts_grp->applyJacobianTransposeInverse(*lsParams, (*fdfdp)[1], 
					      *aVecPtr);
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
    }
    else
      *aVecPtr = *bVecPtr;

    // Scale a and b to unit norm
    aVecPtr->scale(std::sqrt(static_cast<double>(aVecPtr->length())) / aVecPtr->norm());
    bVecPtr->scale(std::sqrt(static_cast<double>(bVecPtr->length())) / bVecPtr->norm());
  }

  else {

    // Get initial "a" vector
    if (!turningPointParams->isParameter("Initial A Vector")) {
      globalData->locaErrorCheck->throwError(callingFunction,
					 "\"Initial A Vector\" is not set!");
    }
    aVecPtr = 
      (*turningPointParams).INVALID_TEMPLATE_QUALIFIER 
      get< Teuchos::RefCountPtr<NOX::Abstract::Vector> >("Initial A Vector");

    // Get initial "b" vector
    if (!isSymmetric) {
      if (!turningPointParams->isParameter("Initial B Vector")) {
	globalData->locaErrorCheck->throwError(callingFunction,
					   "\"Initial B Vector\" is not set!");
      }
      bVecPtr = 
	(*turningPointParams).INVALID_TEMPLATE_QUALIFIER 
        get< Teuchos::RefCountPtr<NOX::Abstract::Vector> >("Initial B Vector");
    }
  }
}
