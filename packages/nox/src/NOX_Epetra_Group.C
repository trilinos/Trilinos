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

#include "NOX_Epetra_Group.H"	// class definition
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_SharedJacobian.H"
#include "NOX_Parameter_List.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_FiniteDifference.H"

// External include files - linking to Aztec00 and Epetra in Trilinos
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Vector.h" 
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

using namespace NOX;
using namespace NOX::Epetra;

Group::Group(const Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_RowMatrix& J):
  xVector(x), // deep copy x     
  RHSVector(x, CopyShape), // new vector of same size
  gradVector(x, CopyShape), // new vector of same size
  NewtonVector(x, CopyShape), // new vector of same size
  tmpVectorPtr(NULL),
  sharedJacobianPtr(new SharedJacobian(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  userInterface(i),
  precType(params.getParameter("Preconditioning Matrix Type","Use Jacobian"))
{
  resetIsValid();

  // Set the Jacobian matrix type
  jacType = getJacobianType();

  // Make sure the preconditioning option is correct for this constructor
  if ((precType != "None") && 
      (precType != "Use Jacobian") && 
      (precType != "User Supplied Action on Vector")) {
    cout << "ERROR:NOX::Epetra::LinearOperator - Preconditioning Matrix Type "
	 << "requested is invalid for this constructor!" << endl;
    throw 1;
  }
}

Group::Group(const Parameter::List& params, Interface& i, 
	     Epetra_Vector& x, Epetra_RowMatrix& J, Epetra_RowMatrix& M):
  xVector(x), // deep copy x     
  RHSVector(x, CopyShape), // new vector of same size
  gradVector(x, CopyShape), // new vector of same size
  NewtonVector(x, CopyShape), // new vector of same size
  tmpVectorPtr(NULL),
  sharedJacobianPtr(new SharedJacobian(J, M)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  userInterface(i),
  precType(params.getParameter("Preconditioning Matrix Type","Use Jacobian"))
{
  resetIsValid();

  // Set the Jacobian matrix type
  jacType = getJacobianType();

  // Make sure the preconditioning option is correct for this constructor
  if ((precType != "None") && 
      (precType != "Use Jacobian") && 
      (precType != "User Supplied Action on Vector")) {
    cout << "ERROR:NOX::Epetra::LinearOperator - Preconditioning Matrix Type "
	 << "requested is invalid for this constructor!" << endl;
    throw 1;
  }
}

Group::Group(const Group& source, CopyType type) :
  xVector(source.xVector.getEpetraVector(), type), 
  RHSVector(source.RHSVector.getEpetraVector(), type), 
  gradVector(source.gradVector.getEpetraVector(), type), 
  NewtonVector(source.NewtonVector.getEpetraVector(), type),
  tmpVectorPtr(NULL),
  sharedJacobianPtr(NULL),
  sharedJacobian(source.sharedJacobian),
  userInterface(source.userInterface),
  jacType(source.jacType),
  precType(source.precType)
{
 
  switch (type) {
    
  case DeepCopy:
    
    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    normRHS = source.normRHS;
    
    // New copy takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedJacobian.getJacobian(this);

    break;

  case CopyShape:
    resetIsValid();
    break;

  default:
    cerr << "ERROR: Invalid ConstructorType for group copy constructor." << endl;
    throw "NOX Error";
  }

}

Group::~Group() 
{
  delete tmpVectorPtr;
  delete sharedJacobianPtr;
}

void Group::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
}


void Group::setAztecOptions(const Parameter::List& p, AztecOO& aztec)
{

  // Preconditioning Matrix Type 
  if (p.isParameterEqual("Preconditioning Matrix Type", "None"))
    aztec.SetAztecOption(AZ_precond, AZ_none);

  // Preconditioning options
  if (p.isParameter("Preconditioning")) {
    
    string prec;
    p.getParameter("Preconditioning",prec);
    
    if (prec == "None") 
      aztec.SetAztecOption(AZ_precond, AZ_none);
    
    else if (prec == "ilut") 
      aztec.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    
    else if (prec == "Polynomial") {
      aztec.SetAztecOption(AZ_precond,AZ_Neumann);
      aztec.SetAztecOption(AZ_poly_ord,3);
    }
    else {
      cout << "ERROR: NOX::Epetra::LinearOperator::setAztecOptions" << endl
	   << "Preconditioning flag is invalid!" << endl;
      throw;
    }
  }

  // Frequency of linear solve residual output
  if (p.isParameter("Iteration Output Frequency")) {
    int i = p.getParameter("Iteration Output Frequency", AZ_last);
    aztec.SetAztecOption(AZ_output, i);
  }
}

Abstract::Group* Group::clone(CopyType type) const 
{
  Group* newgrp = new Group(*this, type);
  return newgrp;
}

Abstract::Group& Group::operator=(const Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Group& source)
{
  // Copy the xVector
  xVector = source.xVector;

  // Copy reference to sharedJacobian
  sharedJacobian = source.sharedJacobian;

  // Update the isValidVectors
  isValidRHS = source.isValidRHS;
  isValidGrad = source.isValidGrad;
  isValidNewton = source.isValidNewton;
  isValidJacobian = source.isValidJacobian;

  // Only copy vectors that are valid
  if (isValidRHS) {
    RHSVector = source.RHSVector;
    normRHS = source.normRHS;
  }

  if (isValidGrad)
    gradVector = source.gradVector;

  if (isValidNewton)
    NewtonVector = source.NewtonVector;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedJacobian.getJacobian(this);
    
  // Copy linear solver options
  jacType = source.jacType;
  precType = source.precType;

  return *this;
}

bool Group::computeX(const Abstract::Group& grp, 
					const Abstract::Vector& d, 
					double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& epetragrp = dynamic_cast<const Group&> (grp);
  const Vector& epetrad = dynamic_cast<const Vector&> (d);
  return computeX(epetragrp, epetrad, step); 
}

bool Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return true;
}

bool Group::computeRHS() 
{
  if (isRHS())
    return true;

  bool status = false;
  
  status = userInterface.computeRHS(xVector.getEpetraVector(), RHSVector.getEpetraVector());

  if (status == false) {
    cout << "ERROR: Epetra::Group::computeRHS() - fill failed!!!"
	 << endl;
    throw "NOX Error";
  } 

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return status;
}

bool Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return true;

  // Take ownership of the Jacobian and get a reference
  Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian(this);

  // Fill the Jacobian
  // Need to add a check here to find out if computing the Jacobian was successful.  
  bool status = false;
  
  if (jacType == "User Supplied") {
    status = userInterface.computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else if (jacType == "Finite Difference") {
    status = (dynamic_cast<FiniteDifference&>(Jacobian)).computeJacobian(xVector.getEpetraVector(), Jacobian);
  }
  else if (jacType == "Matrix-Free") {
    status = (dynamic_cast<MatrixFree&>(Jacobian)).computeJacobian(xVector.getEpetraVector(), Jacobian);
  }

  if (status == false) {
    cout << "ERROR: Epetra::Group::computeJacobian() - fill failed!!!"
	 << endl;
    throw "NOX Error";
  } 

  // Update status
  isValidJacobian = true;

  return status;
}

bool Group::computeGrad() 
{
  if (isGrad())
    return true;
  
  if (!isRHS()) {
    cerr << "ERROR: NOX::Epetra::Group::computeGrad() - RHS is out of date wrt X!" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Epetra::Group::computeGrad() - Jacobian is out of date wrt X!" << endl;
    throw "NOX Error";
  }
  
  // Get a reference to the Jacobian (it's validity was checked above)
  const Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian();

  // Compute grad = Jacobian^T * RHS.
  bool Transpose = true;
  Jacobian.Multiply(Transpose, RHSVector.getEpetraVector(), gradVector.getEpetraVector());

  // Update state
  isValidGrad = true;

  // Return result
  return true;
}

bool Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return true;

  if (!isRHS()) {
    cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid RHS" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Epetra::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }
  
  // Get the Jacobian 
  /* (Have to get non-const version which requires reasserting
     ownership. This is due to a flaw in Epetra_LinearProblem). */
  Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian(this);

  // Create Epetra problem for the linear solve
  Epetra_LinearProblem Problem(&Jacobian, 
  			       &(NewtonVector.getEpetraVector()), 
			       &(RHSVector.getEpetraVector()));

  // Set the default Problem parameters to "hard" (this sets Aztec defaults
  // during the AztecOO instantiation)
  Problem.SetPDL(hard);

  // Create aztec problem
  AztecOO aztec(Problem);  

  // Set specific Aztec parameters requested by NOX
  setAztecOptions(p, aztec);
  
  int maxit = p.getParameter("Max Iterations", 400);
  double tol = p.getParameter("Tolerance", 1.0e-6);

  // Solve Aztec problem
  aztec.Iterate(maxit, tol);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  // Return solution
  return true;
}

bool Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobian(epetrainput, epetraresult);
}

bool Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian (it's validity was check above)
  const Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  bool NoTranspose = false;
  Jacobian.Multiply(NoTranspose, input.getEpetraVector(), result.getEpetraVector());

  return true;
}

bool Group::applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianDiagonalInverse(epetrainput, epetraresult);
}

bool Group::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const
{
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian 
  const Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian();

  // Get epetra reference to the result vector
  Epetra_Vector& r = result.getEpetraVector();

  // Allocate the extra tmpVectorPtr if necessary
  if (tmpVectorPtr == NULL)
    tmpVectorPtr = new Epetra_Vector(r.Map());

  // Get the reference to the temporary vector
  Epetra_Vector& tmpVector = *tmpVectorPtr;

  // Put a copy of the diagonal of the Jacobian into tmpVector
  int retcode = Jacobian.ExtractDiagonalCopy(tmpVector);
  
  // Take element-wise absolute value of diagonal vector
  retcode = tmpVector.Abs(tmpVector);
  
  // Check minimum absolute value of diagonal vector
  double minAbsValue = 0;
  retcode = tmpVector.MinValue(&minAbsValue);

  if(minAbsValue <= 1.e-6) // This minimum threshold can be adjusted
  {
    cout << "Poor scaling on Jacobian diagonal (min abs value: " <<
             minAbsValue << " ) --> NO nonlinear Preconditioning !!" << endl;
    return false; 
  }
  
  // Check if ExtractDiagonalCopy is supported
  if (retcode != 0)
    return false;

  // Calculate r = input ./ tmpVector (./ is element-by-element divide)
  retcode = r.ReciprocalMultiply(1.0, tmpVector, input.getEpetraVector(), 0.0);

  // Check if this worked
  if (retcode != 0)
    return false;

  return true;
}



bool Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& epetrainput = dynamic_cast<const Vector&> (input);
  Vector& epetraresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(epetrainput, epetraresult);
}

bool Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian (it's validity was check above)
  const Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  bool Transpose = true;
  Jacobian.Multiply(Transpose, input.getEpetraVector(), result.getEpetraVector());

  return true;
}



bool Group::isRHS() const 
{   
  return isValidRHS;
}

bool Group::isJacobian() const 
{  
  return ((sharedJacobian.isOwner(this)) && (isValidJacobian));
}

bool Group::isGrad() const 
{   
  return isValidGrad;
}

bool Group::isNewton() const 
{   
  return isValidNewton;
}

const Abstract::Vector& Group::getX() const 
{
  return xVector;
}

const Abstract::Vector& Group::getRHS() const 
{  
  if (!isRHS()) {
    cerr << "ERROR: NOX::Epetra::Group::getRHS() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return RHSVector;
}

double Group::getNormRHS() const
{
  if (!isRHS()) {
    cerr << "ERROR: NOX::Epetra::Group::getNormRHS() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return normRHS;
}

const Abstract::Vector& Group::getGrad() const 
{ 
  if (!isGrad()) {
    cerr << "ERROR: NOX::Epetra::Group::getGrad() - invalid gradient" << endl;
    throw "NOX Error";
  }
    
  return gradVector;
}

const Abstract::Vector& Group::getNewton() const 
{
  if (!isNewton()) {
    cerr << "ERROR: NOX::Epetra::Group::getNewton() - invalid Newton vector" << endl;
    throw "NOX Error";
  }
    
  return NewtonVector;
}


SharedJacobian& Group::getSharedJacobian()
{
  return sharedJacobian;
}

Interface& Group::getUserInterface()
{
  return userInterface;
}


string Group::getJacobianType()
{
  // Get a reference to the Jacobian 
  Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian(this);

  Epetra_RowMatrix* test = 0;
  
  // Is it a "Matrix-Free" Jacobian?
  test = dynamic_cast<MatrixFree*>(&Jacobian);
  if (test != 0) return "Matrix-Free"; 
  
  // Is it a "Finite Difference" Jacobian?
  test = dynamic_cast<FiniteDifference*>(&Jacobian);
  if (test != 0) return "Finite Difference"; 
  
  // Otherwise it must be a "User Supplied" Epetra_RowMatrix.
  return "User Supplied";
}
