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

#include "NOX.H"
#include "NOX_Parameter_List.H"

#include "NOX_Petsc_Group.H"	// class definition
#include "NOX_Petsc_Interface.H"
#include "NOX_Petsc_Vector.H"
#include "NOX_Petsc_SharedJacobian.H"
#include "NOX_Parameter_List.H"

// External include files - linking to Petsc
#include "petscsles.h" 

using namespace NOX;
using namespace NOX::Petsc;

Group::Group(const Parameter::List& params, Interface& i,
             Vec& x, Mat& J, string jType) :
  xVector(x), // deep copy x     
  RHSVector(x, ShapeCopy), // new vector of same size
  gradVector(x, ShapeCopy), // new vector of same size
  NewtonVector(x, ShapeCopy), // new vector of same size
  tmpVectorPtr(NULL),
  sharedJacobianPtr(new SharedJacobian(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  jacType(jType),
  userInterface(i)
{
  resetIsValid();

  // This is material for FD jacobian and is not used for now --> RHooper
  //if(jacType == "Finite Difference") {
  //  int ierr = MatGetColoring(J, MATCOLORING_NATURAL, isColoring);
  //  ierr = MatFDColoringCreate(J, *isColoring, matColoring);
  //  ierr = MatFDColoringSetFunction(matColoring, 
  //              (int (*)(void))ResidualPetscWrapper,
  //              PETSC_NULL);
  //}
}

Group::Group(const Group& source, CopyType type) :
  xVector(source.xVector.getPetscVector(), type), 
  RHSVector(source.RHSVector.getPetscVector(), type), 
  gradVector(source.gradVector.getPetscVector(), type), 
  NewtonVector(source.NewtonVector.getPetscVector(), type),
  tmpVectorPtr(NULL),
  sharedJacobianPtr(NULL),
  sharedJacobian(source.sharedJacobian),
  userInterface(source.userInterface),
  jacType(source.jacType)
{
  switch (type) {
    
  case DeepCopy:
    
    isValidRHS = source.isValidRHS;
    isValidGrad = source.isValidGrad;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    isValidPreconditioner = source.isValidPreconditioner;
    normRHS = source.normRHS;
    
    // New copy takes ownership of the shared Jacobian
    if (isValidJacobian)
      sharedJacobian.getJacobian(this);

    break;

  case ShapeCopy:
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
  isValidPreconditioner = false;
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
  isValidPreconditioner = source.isValidPreconditioner;

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
    
  return *this;
}

bool Group::setX(const Abstract::Vector& y)
{
  return setX(dynamic_cast<const Vector&> (y));
}

bool Group::setX(const Vector& y)
{
  resetIsValid();
  xVector = y;
  return true;
}

bool Group::computeX(const Abstract::Group& grp, 
					const Abstract::Vector& d, 
					double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& petscgrp = dynamic_cast<const Group&> (grp);
  const Vector& petscd = dynamic_cast<const Vector&> (d);
  return computeX(petscgrp, petscd, step); 
}

bool Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return true;
}

bool Group::computeF() 
{
  if (isF())
    return true;

  bool status = false;

  status = userInterface.computeF(xVector.getPetscVector(), RHSVector.getPetscVector());

  if(status == false) {
    cout << "ERROR: Petsc::Group::computeF() - fill failed!!!"
         << endl;
    throw "NOX Error: Fill Failed";
  }

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return true;
}

bool Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return true;

  // Take ownership of the Jacobian and get a reference
  Mat& Jacobian = sharedJacobian.getJacobian(this);

  // Fill the Jacobian

  bool status = false;

  if(jacType == "User Supplied") {

    status = userInterface.computeJacobian(xVector.getPetscVector(), Jacobian);

    if (status == false) {
      cout << "ERROR: Petsc::Group::computeJacobian() - fill failed!!!"
           << endl;
      throw "NOX Error: Fill Failed";
    }
  }
  else if(jacType == "Finite Difference") {
    cout << "Finite Difference evaluation not yet supported !!\n\n";
    throw "NOX Error";
  }

  // Update status
  isValidJacobian = true;

  return true;
}

bool Group::computeGradient() 
{
  if (isGradient())
    return true;
  
  if (!isF()) {
    cerr << "ERROR: NOX::Petsc::Group::computeGradient() - RHS is out of date wrt X!" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Petsc::Group::computeGradient() - Jacobian is out of date wrt X!" << endl;
    throw "NOX Error";
  }
  
  // Get a reference to the Jacobian (it's validity was checked above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Compute grad = Jacobian^T * RHS.
  // Need to add a check on ierr
  int ierr = MatMultTranspose(Jacobian, RHSVector.getPetscVector(), gradVector.getPetscVector());

  // Update state
  isValidGrad = true;

  // Return result
  return true;
}

bool Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return true;

  if (!isF()) {
    cerr << "ERROR: NOX::Petsc::Group::computeNewton() - invalid RHS" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Petsc::Group::computeNewton() - invalid Jacobian" << endl;
    throw "NOX Error";
  }
  
  // Get the Jacobian
  Mat& Jacobian = sharedJacobian.getJacobian(this);

  // Create Petsc SLES problem for the linear solve

  SLES sles;

  int ierr = SLESCreate(PETSC_COMM_WORLD,&sles);CHKERRQ(ierr);

  //ierr = SLESGetKSP(sles, &ksp);CHKERRQ(ierr);
  //ierr = SLESGetPC(sles, &pc);CHKERRQ(ierr);

  //ierr = PCLUSetUseInPlace(pc);CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp, tol, tol, PETSC_DEFAULT, 10);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = SLESSetOperators(sles,Jacobian,Jacobian,
                       DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /*
     Set runtime options (e.g., -ksp_type <type> -pc_type <type>)
  */

  ierr = SLESSetFromOptions(sles);CHKERRQ(ierr);

  /*
     Solve linear system.  Here we explicitly call SLESSetUp() for more
     detailed performance monitoring of certain preconditioners, such
     as ICC and ILU.  This call is optional, as SLESSetUp() will
     automatically be called within SLESSolve() if it hasn't been
     called already.
  */
  ierr = SLESSetUp(sles,RHSVector.getPetscVector(),
                        NewtonVector.getPetscVector());CHKERRQ(ierr);
  int its = 0;
  // Need to put a check on ierr.
  ierr = SLESSolve(sles,RHSVector.getPetscVector(),
                        NewtonVector.getPetscVector(),&its);CHKERRQ(ierr);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  // Return true for now
  return true;
}

bool Group::computePreconditioner()
{
  cout << "NOX::Petsc::Group::computePreconditioner() - Not yet implemented!" << endl;
  exit(0);
  return false;
}

bool Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobian(petscinput, petscresult);
}

bool Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian (it's validity was check above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  MatMult(Jacobian, input.getPetscVector(), result.getPetscVector());

  return true;
}

bool Group::applyJacobianDiagonalInverse(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobianDiagonalInverse(petscinput, petscresult);
}

bool Group::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const
{
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian 
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Get petsc reference to the result vector
  Vec& r = result.getPetscVector();

  // Allocate the extra tmpVectorPtr if necessary
  if (tmpVectorPtr == NULL)
    tmpVectorPtr = new Vec;
    VecDuplicate(r, tmpVectorPtr);

  // Get the reference to the temporary vector
  Vec& tmpVector = *tmpVectorPtr;

  // Put a copy of the diagonal of the Jacobian into tmpVector
  int ierr = MatGetDiagonal(Jacobian, tmpVector); // chk error
  
/*
  // Take element-wise absolute value of diagonal vector
  ierr = VecAbs(tmpVector);  // Check that this works, RH
  
  // Check minimum absolute value of diagonal vector
  int minLocation = 0;
  double minAbsValue = 0;
  ierr = VecMin(tmpVector, &minLocation, &minAbsValue);

  if(minAbsValue <= 1.e-6) // This minimum threshold can be adjusted
  {
    cout << "Poor scaling on Jacobian diagonal (min abs value: " <<
             minAbsValue << " ) --> NO nonlinear Preconditioning !!" << endl;
    return false; 
  }
*/
  
  // Calculate r = input ./ tmpVector (./ is element-by-element divide)
  ierr = VecPointwiseDivide(input.getPetscVector(), tmpVector, r);

  return true;
}


bool Group::applyPreconditionerInverse(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyPreconditionerInverse(petscinput, petscresult);
}

bool Group::applyPreconditionerInverse(const Vector& input, Vector& result) const
{
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian 
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Get petsc reference to the result vector
  Vec& r = result.getPetscVector();

  // Set up preconditioner context
  PC pc;

  int ierr = PCCreate(PETSC_COMM_WORLD,&pc);CHKERRQ(ierr);

  // This sets the type from -pc_type in .petscrc
  // An alternative is to set it explivitly using PCSetType(pc, PCType)

  ierr = PCSetFromOptions(pc);CHKERRQ(ierr);

  /*
     Set operators and vector. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = PCSetOperators(pc,Jacobian,Jacobian,
                       DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = PCSetVector(pc,r);CHKERRQ(ierr);

  // Apply the preconditioner
  ierr = PCApply(pc,input.getPetscVector(),r);CHKERRQ(ierr);

  // Cleanup
  ierr = PCDestroy(pc);CHKERRQ(ierr);

  return true;
}


bool Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(petscinput, petscresult);
}

bool Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return false;

  // Get a reference to the Jacobian (it's validity was check above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  int ierr = MatMultTranspose(Jacobian, input.getPetscVector(), result.getPetscVector());

  return true;
}


// This is intended for FD Jacobian computation --> RHooper
//int Group::ResidualPetscWrapper(SNES, Vec& xx, Vec& fx, void*)
//{
//  bool status = userInterface.computeF(xx, fx);
//
//  return(0);
//}
  

bool Group::isF() const 
{   
  return isValidRHS;
}

bool Group::isJacobian() const 
{  
  return ((sharedJacobian.isOwner(this)) && (isValidJacobian));
}

bool Group::isGradient() const 
{   
  return isValidGrad;
}

bool Group::isNewton() const 
{   
  return isValidNewton;
}

bool Group::isPreconditioner() const 
{   
  return isValidPreconditioner;
}

const Abstract::Vector& Group::getX() const 
{
  return xVector;
}

const Abstract::Vector& Group::getF() const 
{  
  if (!isF()) {
    cerr << "ERROR: NOX::Petsc::Group::getF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return RHSVector;
}

double Group::getNormF() const
{
  if (!isF()) {
    cerr << "ERROR: NOX::Petsc::Group::getNormF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return normRHS;
}

const Abstract::Vector& Group::getGradient() const 
{ 
  if (!isGradient()) {
    cerr << "ERROR: NOX::Petsc::Group::getGradient() - invalid gradient" << endl;
    throw "NOX Error";
  }
    
  return gradVector;
}

const Abstract::Vector& Group::getNewton() const 
{
  if (!isNewton()) {
    cerr << "ERROR: NOX::Petsc::Group::getNewton() - invalid Newton vector" << endl;
    throw "NOX Error";
  }
    
  return NewtonVector;
}


