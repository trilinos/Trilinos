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

Group::Group(Interface& i, Vec& x, Mat& J) :
  xVector(x, "Solution"), // deep copy x     
  RHSVector(x, "RHS", ShapeCopy), // new vector of same size
  gradVector(x, "Grad", ShapeCopy), // new vector of same size
  NewtonVector(x, "Newton", ShapeCopy), // new vector of same size
  sharedJacobianPtr(new SharedJacobian(J)), // pass J to SharedJacobian
  sharedJacobian(*sharedJacobianPtr), // pass J to SharedJacobian
  jacType("User Supplied"), // the only option for now
  userInterface(i)
{
  resetIsValid();
}

Group::Group(const Group& source, CopyType type) :
  xVector(source.xVector, type), 
  RHSVector(source.RHSVector, type), 
  gradVector(source.gradVector, type), 
  NewtonVector(source.NewtonVector, type),
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

void Group::setX(const Abstract::Vector& y)
{
  setX(dynamic_cast<const Vector&> (y));
  return;
}

void Group::setX(const Vector& y)
{
  resetIsValid();
  xVector = y;
  return;
}

void Group::computeX(const Abstract::Group& grp, 
					const Abstract::Vector& d, 
					double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& petscgrp = dynamic_cast<const Group&> (grp);
  const Vector& petscd = dynamic_cast<const Vector&> (d);
  computeX(petscgrp, petscd, step); 
  return;
}

void Group::computeX(const Group& grp, const Vector& d, double step) 
{
  resetIsValid();
  xVector.update(1.0, grp.xVector, step, d);
  return;
}

Abstract::Group::ReturnType Group::computeF() 
{
  if (isF())
    return Abstract::Group::Ok;

  bool status = false;

  status = userInterface.computeF(xVector.getPetscVector(), RHSVector.getPetscVector());

  if(status == false) {
    cout << "ERROR: Petsc::Group::computeF() - fill failed!!!"
         << endl;
    throw "NOX Error: RHS Fill Failed";
  }

  normRHS = RHSVector.norm();

  isValidRHS = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return Abstract::Group::Ok;

  // Take ownership of the Jacobian and get a reference
  Mat& Jacobian = sharedJacobian.getJacobian(this);

  // Fill the Jacobian

  bool status = false;

  if(jacType == "User Supplied") {

    status = userInterface.computeJacobian(xVector.getPetscVector(), Jacobian);

    if (status == false) {
      cout << "ERROR: Petsc::Group::computeJacobian() - fill failed!!!"
           << endl;
      throw "NOX Error: Jacobian Fill Failed";
    }
  }
  else if(jacType == "Finite Difference") {
    cout << "Finite Difference evaluation not yet supported !!\n\n";
    throw "NOX Error";
  }

  // Update status
  isValidJacobian = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeGradient() 
{
  if (isGradient())
    return Abstract::Group::Ok;
  
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
  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType Group::computeNewton(NOX::Parameter::List& p) 
{
  if (isNewton())
    return Abstract::Group::Ok;

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

  int ierr = SLESCreate(PETSC_COMM_WORLD,&sles);//CHKERRQ(ierr);

  //ierr = SLESGetKSP(sles, &ksp);CHKERRQ(ierr);
  //ierr = SLESGetPC(sles, &pc);CHKERRQ(ierr);

  //ierr = PCLUSetUseInPlace(pc);CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp, tol, tol, PETSC_DEFAULT, 10);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = SLESSetOperators(sles,Jacobian,Jacobian,
                       DIFFERENT_NONZERO_PATTERN);//CHKERRQ(ierr);

  /*
     Set runtime options (e.g., -ksp_type <type> -pc_type <type>)
  */

  ierr = SLESSetFromOptions(sles);//CHKERRQ(ierr);

  /*
     Solve linear system.  Here we explicitly call SLESSetUp() for more
     detailed performance monitoring of certain preconditioners, such
     as ICC and ILU.  This call is optional, as SLESSetUp() will
     automatically be called within SLESSolve() if it hasn't been
     called already.
  */
  ierr = SLESSetUp(sles,RHSVector.getPetscVector(),
                        NewtonVector.getPetscVector());//CHKERRQ(ierr);
  int its = 0;
  // Need to put a check on ierr.
  ierr = SLESSolve(sles,RHSVector.getPetscVector(),
                        NewtonVector.getPetscVector(),&its);//CHKERRQ(ierr);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  return Abstract::Group::Ok;
}

Abstract::Group::ReturnType 
Group::applyJacobian(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobian(petscinput, petscresult);
}

Abstract::Group::ReturnType 
Group::applyJacobian(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian (it's validity was checked above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  MatMult(Jacobian, input.getPetscVector(), result.getPetscVector());

  return Abstract::Group::Ok;
}


Abstract::Group::ReturnType 
Group::applyRightPreconditioning(Parameter::List& params, 
                                 const Abstract::Vector& input, 
                                 Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyRightPreconditioning(petscinput, petscresult);
}

Abstract::Group::ReturnType 
Group::applyRightPreconditioning(const Vector& input, Vector& result) const
{
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian 
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Get petsc reference to the result vector
  Vec& r = result.getPetscVector();

  // Set up preconditioner context
  PC pc;

  int ierr = PCCreate(PETSC_COMM_WORLD,&pc);//CHKERRQ(ierr);

  // Here a default to jacobi (jacobian-diagonal-inverse) is established
  // but can be overridden via specification of pc_type in .petscrc
  //ierr = PCSetType(pc, PCJACOBI);CHKERRQ(ierr);

  // This allows more general preconditioning via specification of -pc_type 
  // in .petscrc
  ierr = PCSetFromOptions(pc);//CHKERRQ(ierr);

  /*
     Set operators and vector. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = PCSetOperators(pc,Jacobian,Jacobian,
                       DIFFERENT_NONZERO_PATTERN);//CHKERRQ(ierr);
  ierr = PCSetVector(pc,r);//CHKERRQ(ierr);

  // Apply the preconditioner
  ierr = PCApply(pc,input.getPetscVector(),r, PC_RIGHT);//CHKERRQ(ierr);

  // Cleanup
  ierr = PCDestroy(pc);//CHKERRQ(ierr);

  return Abstract::Group::Ok;
}


Abstract::Group::ReturnType 
Group::applyJacobianTranspose(const Abstract::Vector& input, Abstract::Vector& result) const
{
  const Vector& petscinput = dynamic_cast<const Vector&> (input);
  Vector& petscresult = dynamic_cast<Vector&> (result);
  return applyJacobianTranspose(petscinput, petscresult);
}

Abstract::Group::ReturnType 
Group::applyJacobianTranspose(const Vector& input, Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
    return Abstract::Group::BadDependency;

  // Get a reference to the Jacobian (it's validity was check above)
  const Mat& Jacobian = sharedJacobian.getJacobian();

  // Apply the Jacobian
  int ierr = MatMultTranspose(Jacobian, input.getPetscVector(), result.getPetscVector());

  return Abstract::Group::Ok;
}


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


