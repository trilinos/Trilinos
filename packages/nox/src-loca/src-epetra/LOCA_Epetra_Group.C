// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Epetra_Group.H"	          // class definition

#include "LOCA_Epetra_Interface.H"        // class data members
#include "NOX_Epetra_SharedOperator.H"
#include "NOX_Parameter_List.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_Operator.H"

// External include files - linking to Aztec00 and Epetra in Trilinos
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Vector.h" 
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"


using namespace LOCA;
using namespace LOCA::Epetra;

Group::Group(const NOX::Parameter::List& params, Interface& i, 
	     const ParameterVector& p, Vector& x, Epetra_Operator& J) :
  NOX::Epetra::Group(params, i, x, J),
  pVectorPtr(new ParameterVector(p)),
  userInterface(i)
{
}

Group::Group(const NOX::Parameter::List& params, Interface& i, 
	     const ParameterVector& p, Vector& x, Epetra_Operator& J, 
	     Epetra_Operator& M) :
  NOX::Epetra::Group(params, i, x, J, M),
  pVectorPtr(new ParameterVector(p)),
  userInterface(i)
{
}

Group::Group(const Group& source, NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  pVectorPtr(new ParameterVector(*(source.pVectorPtr))),
  userInterface(source.userInterface)
{
}

Group::~Group() 
{
  delete pVectorPtr;
}

Abstract::Group* Group::clone(NOX::CopyType type) const 
{
  Group* newgrp = new Group(*this, type);
  return newgrp;
}

Abstract::Group& Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Group& source)
{
  *pVectorPtr = *source.pVectorPtr;
  NOX::Epetra::Group::operator=(source);
  return *this;
}

bool Group::setX(const NOX::Abstract::Vector& y) 
{
  return setX(dynamic_cast<const LOCA::Epetra::Vector&>(y));
}

bool Group::setX(const Vector& y) 
{
  return NOX::Epetra::Group::setX(y);
}

bool Group::computeX(const NOX::Abstract::Group& grp, 
		     const NOX::Abstract::Vector& d, 
		     double step) 
{
  const LOCA::Abstract::Group& epetraGrp = dynamic_cast<const LOCA::Abstract::Group&>(grp);
  const LOCA::Abstract::Vector& epetraD = dynamic_cast<const LOCA::Abstract::Vector&>(d);
  return NOX::Epetra::Group::computeX(epetraGrp, epetraD, step); 
}

bool Group::computeX(const Group& grp, const Vector& d, double step) 
{
  return NOX::Epetra::Group::computeX(grp, d, step);
}

bool Group::setParams(const ParameterVector& p)
{
  *pVectorPtr = p;
  resetIsValid();
  return true;
}

bool Group::computeParams(const ParameterVector& oldParams,
			  const ParameterVector& direction, 
			  double step)
{
  *pVectorPtr = oldParams;
  pVectorPtr->update(step, direction, 1.0);
  resetIsValid();
  return true;
}


bool Group::computeF() 
{

  bool status = false;
  
  // Set the parameters prior to computing F
  status = userInterface.setParameters(*pVectorPtr);

  if (status == false)
    return false;
  
  return NOX::Epetra::Group::computeF();;
}

bool Group::computeJacobian() 
{
  bool status = false;
  
  // Set the parameters prior to computing F
  status = userInterface.setParameters(*pVectorPtr);

  if (status == false)
    return false;
  
  return NOX::Epetra::Group::computeJacobian();
}

bool Group::computeGradient() 
{
  return NOX::Epetra::Group::computeGradient();
}

bool Group::computeNewton(NOX::Parameter::List& p) 
{
  return NOX::Epetra::Group::computeNewton(p);
}

bool Group::applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobian(epetrainput, epetraresult);
}

bool Group::applyJacobian(const Vector& input, Vector& result) const
{
  return 
    NOX::Epetra::Group::applyJacobian(input, result);
}

bool Group::applyJacobianInverse (NOX::Parameter::List &p, const NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const
{
  const Vector& epetraInput = dynamic_cast<const Vector&>(input);
  Vector& epetraResult = dynamic_cast<Vector&>(result);
  return applyJacobianInverse(p, epetraInput, epetraResult);
}

bool Group::applyJacobianInverse (NOX::Parameter::List &p, const Vector &input, Vector &result) const
{
  return NOX::Epetra::Group::applyJacobianInverse(p, input, result);
}

bool Group::applyJacobianDiagonalInverse(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianDiagonalInverse(epetrainput, epetraresult);
}

bool Group::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const
{
  return NOX::Epetra::Group::applyJacobianDiagonalInverse(input, result);
}

bool Group::applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(epetrainput, epetraresult);
}

bool Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  return NOX::Epetra::Group::applyJacobianTranspose(input, result);
}

bool Group::applyRightPreconditioning(NOX::Parameter::List& params,
				      const NOX::Abstract::Vector& input, 
				      NOX::Abstract::Vector& result) const
{
  const Vector& epetraInput = dynamic_cast<const Vector&>(input);
  Vector& epetraResult = dynamic_cast<Vector&>(result);
  
  return applyRightPreconditioning(params, epetraInput, epetraResult);
}

bool Group::applyRightPreconditioning(NOX::Parameter::List& params,
				      const Vector& input, 
				      Vector& result) const
{
  return NOX::Epetra::Group::applyRightPreconditioning(params, input, result);
}

bool Group::isF() const 
{   
  return NOX::Epetra::Group::isF();
}

bool Group::isJacobian() const 
{  
  return NOX::Epetra::Group::isJacobian();
}

bool Group::isGradient() const 
{   
  return NOX::Epetra::Group::isGradient();
}

bool Group::isNewton() const 
{   
  return NOX::Epetra::Group::isNewton();
}

bool Group::isNormNewtonSolveResidual() const 
{   
  return NOX::Epetra::Group::isNormNewtonSolveResidual();
}

bool Group::isPreconditioner() const 
{  
  return NOX::Epetra::Group::isPreconditioner();
}

const Abstract::Vector& Group::getX() const 
{
  const LOCA::Abstract::Vector& locaXVector = dynamic_cast<const LOCA::Abstract::Vector&>(xVector);
  return locaXVector;
}

const ParameterVector& Group::getParams() const 
{
  return *pVectorPtr;
}

const Abstract::Vector& Group::getF() const 
{  
  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::getF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  const LOCA::Abstract::Vector& locaRHSVector = dynamic_cast<const LOCA::Abstract::Vector&>(RHSVector);
  return locaRHSVector;
}

double Group::getNormF() const
{
  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::getNormF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return normRHS;
}

const Abstract::Vector& Group::getGradient() const 
{ 
  if (!isGradient()) {
    cerr << "ERROR: NOX::Epetra::Group::getGradient() - invalid gradient" << endl;
    throw "NOX Error";
  }
    
  const LOCA::Abstract::Vector& locaGradVector = dynamic_cast<const LOCA::Abstract::Vector&>(gradVector);
  return locaGradVector;
}

const Abstract::Vector& Group::getNewton() const 
{
  if (!isNewton()) {
    cerr << "ERROR: NOX::Epetra::Group::getNewton() - invalid Newton vector" << endl;
    throw "NOX Error";
  }
    
  const LOCA::Abstract::Vector& locaNewtonVector = dynamic_cast<const LOCA::Abstract::Vector&>(NewtonVector);
  return locaNewtonVector;
}

double Group::getNormNewtonSolveResidual() const
{
  // Make sure value is not already calculated
  if (isValidNormNewtonSolveResidual) 
    return normNewtonSolveResidual;

  // Otherwise throw an error since a Newton direction has not been calculated
  // wrt this solution group
  /*
  cout << "ERROR: NOX::Epetra::Group::getNormNewtonSolveResidual() - Group has "
       << "not performed a Newton solve corresponding to this solution vector!"
       << endl;
  throw "NOX Error";
  */
  return -1.0;
}  

NOX::Epetra::SharedOperator& Group::getSharedJacobian()
{
  return sharedJacobian;
}

NOX::Epetra::SharedOperator& Group::getSharedPreconditioner()
{
  return sharedPreconditioner;
}

NOX::Epetra::Interface& Group::getUserInterface()
{
  return userInterface;
}
