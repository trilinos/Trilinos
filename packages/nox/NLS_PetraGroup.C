// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_PetraGroup.H"

NLS_PetraGroup::NLS_PetraGroup(Epetra_Vector& x, NLS_PetraGroupInterface& I) :
  xVector(x, true), // deep copy x     
  RHSVector(x, false), // new vector of same size
  gradVector(x, false), // new vector of same size
  NewtonVector(x, false) // new vector of same size
{
  Jac = NULL;
  Interface = &I;
  reset();
}

NLS_PetraGroup::NLS_PetraGroup(Epetra_Vector& x, Epetra_RowMatrix& J, NLS_PetraGroupInterface& I) :
  xVector(x, true), // deep copy x     
  RHSVector(x, false), // new vector of same size
  gradVector(x, false), // new vector of same size
  NewtonVector(x, false) // new vector of same size
{
  Jac = &J;
  Interface = &I;
  reset();
}

//! NLS_Group deconstructor
NLS_PetraGroup::~NLS_PetraGroup() 
{
  // Calling program is responsible for destructing Jacobian.
  Jac = NULL;
}

//! Copies the values of all vectors in source group to this group
NLS_Group& NLS_PetraGroup::copy(const NLS_Group& source)
{
  cout << "ERROR: NLS_PetraGroup::copy() - requires a PetraGroup object be passed in!" << endl;
  throw;
}

//! Copies the values of all vectors in source group to this group
NLS_Group& NLS_PetraGroup::copy(const NLS_PetraGroup& copyFrom)
{
  // Update the vectors
  xVector.copy(copyFrom.xVector, 1.0);
  
  // IsValid of the Jacobian is necessarily false since we don't copy
  // that.
  isValidJacobian = false;

  // Update the isValidVectors
  isValidRHS = copyFrom.isValidRHS;
  isValidGrad = copyFrom.isValidGrad;
  isValidNewton = copyFrom.isValidNewton;

  // Only copy vectors that are valid
  if (isValidRHS)
    RHSVector.copy(copyFrom.RHSVector);
  if (isValidGrad)
    gradVector.copy(copyFrom.gradVector);
  if (isValidNewton)
    NewtonVector.copy(copyFrom.NewtonVector);
}

//! Compute and return solution vector
const NLS_Vector& NLS_PetraGroup::computeX(const NLS_Group& grp, const NLS_Vector& d, double step) 
{
  cout << "ERROR: NLS_PetraGroup::computeX() - Pass Petra objects in call!" << endl;
  throw;
}

//! Compute and return solution vector
const NLS_Vector& NLS_PetraGroup::computeX(const NLS_PetraGroup& grp, const NLS_PetraVector& d, 
					   double step) 
{
  xVector.update(grp.getX(),  d, step);
  reset();
  return xVector;
}

//! Compute and return RHS
const NLS_Vector& NLS_PetraGroup::computeRHS() 
{
  Interface->computeRHS(xVector.getPetraVector(), RHSVector.getPetraVector());
}

//! Compute Jacobian
void NLS_PetraGroup::computeJacobian() 
{
  if (!isJacobianEnabled())
    throw;

  Interface->computeJacobian(xVector.getPetraVector(), *Jac);
}

//! Compute and return gradient 
/*! Throws an error if RHS and Jacobian have not been computed */
const NLS_Vector& NLS_PetraGroup::computeGrad() 
{
  if (!isValidRHS) {
    cout << "ERROR: NLS_PetraGroup::computeGrad() - RHS is out of date wrt X!" << endl;
    throw;
  }
  if (!isValidJacobian) {
    cout << "ERROR: NLS_PetraGroup::computeGrad() - Jacobian is out of date wrt X!" << endl;
    throw;
  }

  // Compute grad = Jac^T * RHS.
  bool TransJ = true;
  Jac->Multiply(TransJ, RHSVector.getPetraVector(), gradVector.getPetraVector());
  return gradVector;
}

//! Compute and return Newton direction 
/*! Throws an error if RHS and Jacobian have not been computed */
const NLS_Vector& NLS_PetraGroup::computeNewton() 
{
  cout << "ERROR: No direct methods are avialable for matrix inversion yet!\n"
       << "Use the iterative solver call!" << endl;
  throw;
}

//! Compute and return Newton direction, using desired accuracy for nonlinear solve
/*! Throws an error if RHS and Jacobian have not been computed */
const NLS_Vector& NLS_PetraGroup::computeNewton(NLS_ParameterList& parameter) 
{
  if (!isValidRHS) {
    cout << "ERROR: computeNewton() - RHS is out of date wrt X!" << endl;
    throw;
  }
  else if (!isValidJacobian) {
    cout << "ERROR: computeNewton() - Jacobian is out of date wrt X!" << endl;
    throw;
  }

  Epetra_LinearProblem Problem(Jac, &(xVector.getPetraVector()), &(RHSVector.getPetraVector()));
  // For now, set problem level to hard, moderate, or easy
  Problem.SetPDL(hard);
  //AztecOO aztec(Problem);
  //aztec.Iterate(nlParamsPtr->getMaxLinearStep(), eta_k);
  return NewtonVector;
}

/** @name Checks to see if various objects have been computed. 
 *
 * Returns true if the corresponding "compute" function has been
 * called since the last update to the solution vector (via
 * instantiation or computeX). */

bool NLS_PetraGroup::isRHS() const 
{   
  return isValidRHS;
}
bool NLS_PetraGroup::isJacobian() const 
{  
  return isValidJacobian;
}
bool NLS_PetraGroup::isGrad() const 
{   
  return isValidGrad;
}
bool NLS_PetraGroup::isNewton() const 
{   
  return isValidNewton;
}

//! Return solution vector
const NLS_Vector& NLS_PetraGroup::getX() const 
{
  return xVector;
}

//! Return rhs
const NLS_Vector& NLS_PetraGroup::getRHS() const 
{  
  return RHSVector;
}

//! Return gradient 
const NLS_Vector& NLS_PetraGroup::getGrad() const 
{ 
  return gradVector;
}

//! Return Newton direction 
const NLS_Vector& NLS_PetraGroup::getNewton() const 
{
  return NewtonVector;
}

NLS_Group* NLS_PetraGroup::newGroup() const
{
  // Not implemented
  return NULL;
}

void NLS_PetraGroup::reset()
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
}

bool NLS_PetraGroup::isJacobianEnabled() const
{
  return (Jac != NULL);
}
