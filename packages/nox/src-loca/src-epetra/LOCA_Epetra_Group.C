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

#include "LOCA_Epetra_Group.H"	          // class definition

#include "LOCA_Epetra_Interface_Required.H"        // class data members
#include "NOX_Parameter_List.H"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Epetra_ShiftInvertOperator.H"
#include "AztecOO.h"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

LOCA::Epetra::Group::Group(
	    NOX::Parameter::List& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess),
  LOCA::Abstract::Group(),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  tmpVectorPtr2(),
  scaleVecPtr()
{
}

LOCA::Epetra::Group::Group(
	    NOX::Parameter::List& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	    const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(),
  params(p),
  userInterface(i),
  userInterfaceTime(),
  tmpVectorPtr2(),
  scaleVecPtr()
{
}

LOCA::Epetra::Group::Group(
	NOX::Parameter::List& printingParams, 
	const Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependent>& i, 
	NOX::Epetra::Vector& initialGuess, 
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	const LOCA::ParameterVector& p) :
  NOX::Epetra::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(),
  params(p),
  userInterface(i), 
  userInterfaceTime(i),
  tmpVectorPtr2(),
  scaleVecPtr()
{
}

LOCA::Epetra::Group::Group(const LOCA::Epetra::Group& source, 
			   NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  params(source.params),
  userInterface(source.userInterface),
  userInterfaceTime(source.userInterfaceTime),
  tmpVectorPtr2(),
  scaleVecPtr()
{
  if (source.scaleVecPtr != Teuchos::null)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
}

LOCA::Epetra::Group::~Group() 
{
}

Teuchos::RefCountPtr<NOX::Abstract::Group>
LOCA::Epetra::Group::clone(NOX::CopyType type) const 
{
  return Teuchos::rcp(new LOCA::Epetra::Group(*this, type));
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Abstract::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Epetra::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Epetra::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Epetra::Group& source)
{
  params = source.params;
  NOX::Epetra::Group::operator=(source);
  LOCA::Abstract::Group::operator=(source);
  params = source.params;
  userInterface = source.userInterface;
  userInterfaceTime = source.userInterfaceTime;
  if (source.scaleVecPtr != Teuchos::null)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
  return *this;
}

void 
LOCA::Epetra::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::Epetra::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::Epetra::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeF() 
{

  if (isF())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface->setParameters(params);
  
  return NOX::Epetra::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::computeJacobian() 
{

  if (isJacobian())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface->setParameters(params);

  return NOX::Epetra::Group::computeJacobian();
}

const LOCA::ParameterVector& 
LOCA::Epetra::Group::getParams() const 
{
  return params;
}

Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>
LOCA::Epetra::Group::getUserInterface()
{
  return userInterface;
}

void
LOCA::Epetra::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Epetra::Vector& x_,
				      const double conParam) const
{
  userInterface->printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Abstract::Vector& x_,
				      const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

double
LOCA::Epetra::Group::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  if (scaleVecPtr == Teuchos::null)
    return a.innerProduct(b) / a.length();
  else {
    Teuchos::RefCountPtr<NOX::Abstract::Vector> as = a.clone(NOX::DeepCopy);
    Teuchos::RefCountPtr<NOX::Abstract::Vector> bs = b.clone(NOX::DeepCopy);
    double d;

    as->scale(*scaleVecPtr);
    bs->scale(*scaleVecPtr);
    d = as->innerProduct(*bs);

    return d;
  }
}
   
void
LOCA::Epetra::Group::setScaleVector(const NOX::Abstract::Vector& s)
{
  scaleVecPtr = s.clone(NOX::DeepCopy);

  return;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeMassMatrix()
{
  // Hardwired for matrix-free mode on mass matrix, so this routine does nothing
  if(userInterfaceTime != Teuchos::null)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyMassMatrix(const NOX::Abstract::Vector& input,
                                        NOX::Abstract::Vector& result) const
{
  // For matrix-free mode, call interface for applying mass matrix
  if(userInterfaceTime != Teuchos::null){
     const NOX::Epetra::Vector& epetraInput = 
       dynamic_cast<const NOX::Epetra::Vector&>(input);
     NOX::Epetra::Vector& epetraResult =
       dynamic_cast<NOX::Epetra::Vector&>(result);
     return userInterfaceTime->applyMassMatrix(epetraInput,epetraResult);
  }
  else
    return NOX::Abstract::Group::BadDependency;
}

bool 
LOCA::Epetra::Group::isMassMatrix()
{
  return true;
}


NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrix(const NOX::Abstract::Vector& input,
                                           NOX::Abstract::Vector& result,
                                           double shift) const
{
  return applyShiftedMatrix(dynamic_cast<const NOX::Epetra::Vector&>(input), 
			    dynamic_cast<NOX::Epetra::Vector&>(result), shift);
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrix(
				       const NOX::Epetra::Vector& epetraInput,
				       NOX::Epetra::Vector& epetraResult,
				       double shift) const
{
  // Hardwired for matrix-free mode, so no shifted matrix is available
  // Create action of shifted matrix in 2 steps

  applyJacobian(epetraInput,epetraResult);

  if(userInterfaceTime != Teuchos::null) {

   // If there is a mass matrix, shifted matrix is J + shift*M

    NOX::Epetra::Vector massResult(epetraResult, NOX::ShapeCopy);

    userInterfaceTime->applyMassMatrix(epetraInput, massResult);

    epetraResult.update(shift, massResult, 1.00);
  }

  else {

    // If no mass matrix, shifted matrix is J + shift*I

     epetraResult.update(shift, epetraInput, 1.00);
  }

  return NOX::Abstract::Group::Ok;
}


NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrixInverse(
					    NOX::Parameter::List& params,
					    const NOX::Abstract::Vector& input,
					    NOX::Abstract::Vector& result,
					    double shift)
                                           
{
  const NOX::Epetra::Vector& epetraInput = 
    dynamic_cast<const NOX::Epetra::Vector&>(input);
  NOX::Epetra::Vector& epetraResult = 
    dynamic_cast<NOX::Epetra::Vector&>(result);

  // If shift is zero, just apply Jacobian inverse

  if(shift == 0.0) {
    applyJacobianInverse(params, epetraInput, epetraResult);
  }
  else {

    // Otherwise, construct a shift and invert operator, and use AztecOO to 
    // solve linear system

    Teuchos::RefCountPtr<LOCA::Epetra::ShiftInvertOperator> A =
      Teuchos::rcp(new LOCA::Epetra::ShiftInvertOperator(
		  Teuchos::rcp(this,false),
		  sharedLinearSystem.getObject(this)->getJacobianOperator(),
		  shift,
		  userInterfaceTime != Teuchos::null));

    NOX::Epetra::Vector dummy(epetraResult, NOX::ShapeCopy);
    Epetra_Vector& epetra_dummy = dummy.getEpetraVector();    
    Teuchos::RefCountPtr<LOCA::Epetra::ShiftInvertInterface> interface = 
      Teuchos::rcp(new LOCA::Epetra::ShiftInvertInterface); 
    NOX::Parameter::List& solveList = params.sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver");

    NOX::Epetra::LinearSystemAztecOO shiftsys(
	params,
	solveList,
	userInterface,
	Teuchos::rcp_dynamic_cast<NOX::Epetra::Interface::Jacobian>(interface),
	A,
	epetra_dummy); 

    shiftsys.setJacobianOperatorForSolve(A);

    shiftsys.applyJacobianInverse(params,epetraInput,epetraResult);

  }

  return NOX::Abstract::Group::Ok;  
}

void
LOCA::Epetra::Group::scaleVector(NOX::Abstract::Vector& x) const
{
  if (scaleVecPtr == Teuchos::null)
    x.scale(1.0 / sqrt(static_cast<double>(x.length())));
  else 
    x.scale(*scaleVecPtr);
}

void
LOCA::Epetra::Group::projectToDraw(const NOX::Abstract::Vector& x,
				   double *px) const
{
  const NOX::Epetra::Vector& ex = 
    dynamic_cast<const NOX::Epetra::Vector&>(x);
  userInterface->projectToDraw(ex, px);
}

int
LOCA::Epetra::Group::projectToDrawDimension() const
{
  return userInterface->projectToDrawDimension();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::augmentJacobianForHomotopy(double conParamValue)
{

  //Allocate temporary vectors if not aready done
  if (tmpVectorPtr == Teuchos::null)
    tmpVectorPtr = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));
  if (tmpVectorPtr2 == Teuchos::null)
    tmpVectorPtr2 = Teuchos::rcp(new Epetra_Vector(xVector.getEpetraVector()));

  tmpVectorPtr2->PutScalar(1.0-conParamValue);

  // See if it is an Epetra_CrsMatrix
  Teuchos::RefCountPtr<const Epetra_CrsMatrix> constTestCrs;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> testCrs;
  constTestCrs = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>
    (sharedLinearSystem.getObject(this)->getJacobianOperator());
  if (constTestCrs != Teuchos::null) {
    testCrs = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(constTestCrs);
    testCrs->Scale(conParamValue);
    testCrs->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testCrs->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;

  }

  // See if it is an Epetra_VbrMatrix
  Teuchos::RefCountPtr<const Epetra_VbrMatrix> constTestVbr;
  Teuchos::RefCountPtr<Epetra_VbrMatrix> testVbr;
  constTestVbr = Teuchos::rcp_dynamic_cast<const Epetra_VbrMatrix>
    (sharedLinearSystem.getObject(this)->getJacobianOperator());
  if (constTestVbr != Teuchos::null) {
    testVbr = Teuchos::rcp_const_cast<Epetra_VbrMatrix>(constTestVbr);
    testVbr->Scale(conParamValue);
    testVbr->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testVbr->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;
  }

  // Otherwise this alg won't work -- return NotDefined

  return LOCA::Abstract::Group::NotDefined;
}

void
LOCA::Epetra::Group::setJacobianOperatorForSolve(
		  const Teuchos::RefCountPtr<const Epetra_Operator>& op) const
{
  // Set Jacobian operator for solve
  sharedLinearSystem.getObject(this)->setJacobianOperatorForSolve(op);
  isValidSolverJacOp = true;
}
