// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "AztecOO.h"
#include "NLS_Utilities.H"
#include "NLS_PetraGroup.H"

NLS_PetraGroup::NLS_PetraGroup(Epetra_Vector& x, Epetra_RowMatrix& J, NLS_PetraGroupInterface& I) :
  xVector(x, true), // deep copy x     
  RHSVector(x, false), // new vector of same size
  gradVector(x, false), // new vector of same size
  NewtonVector(x, false), // new vector of same size
  Interface(I) // set reference
{
  Jac = &J;
  reset();
}

NLS_PetraGroup::NLS_PetraGroup(const NLS_PetraGroup& copyFrom) :
  xVector(copyFrom.xVector.getPetraVector(), true), 
  RHSVector(copyFrom.RHSVector.getPetraVector(), true), 
  gradVector(copyFrom.gradVector.getPetraVector(), true), 
  NewtonVector(copyFrom.NewtonVector.getPetraVector(), true),
  Interface(copyFrom.Interface)
{
  isValidRHS = copyFrom.isValidRHS;
  isValidGrad = copyFrom.isValidGrad;
  isValidNewton = copyFrom.isValidNewton;

  Jac = NULL;
  isValidJacobian = false;
}

NLS_PetraGroup::~NLS_PetraGroup() 
{
}

void NLS_PetraGroup::reset() //private
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

NLS_Group* NLS_PetraGroup::newCopy(bool isJacobianEnabled) const 
{
  if (isJacobianEnabled)
    return NULL;

  NLS_PetraGroup* newgrp = new NLS_PetraGroup(*this);
  return newgrp;
}

NLS_Group& NLS_PetraGroup::copy(const NLS_Group& source)
{
  copy(dynamic_cast<const NLS_PetraGroup&> (source));
  return *this;
}

NLS_Group& NLS_PetraGroup::copy(const NLS_PetraGroup& copyFrom)
{
  // Update the xVector
  xVector.copy(copyFrom.xVector);
  
  // Don't copy Jacobian
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
  return *this;
}

const NLS_Vector& NLS_PetraGroup::computeX(const NLS_Group& grp, const NLS_Vector& d, double step) 
{
  computeX(dynamic_cast<const NLS_PetraGroup&> (grp),
	   dynamic_cast<const NLS_PetraVector&> (d),
	   step); 
  return xVector;
}

//! Compute and return solution vector
const NLS_Vector& NLS_PetraGroup::computeX(const NLS_PetraGroup& grp, 
					   const NLS_PetraVector& d, 
					   double step) 
{
  xVector.update(dynamic_cast<const NLS_PetraVector&> (grp.getX()),  d, step);
  reset();
  return xVector;
}

const NLS_Vector& NLS_PetraGroup::computeRHS() 
{
  Interface.computeRHS(xVector.getPetraVector(), RHSVector.getPetraVector());
  isValidRHS = true;
  return RHSVector;
}

void NLS_PetraGroup::computeJacobian() 
{
  if (!isJacobianEnabled())
    throw;
  Interface.computeJacobian(xVector.getPetraVector(), *Jac);
  isValidJacobian = true;
}

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
  isValidGrad = true;
  return gradVector;
}

const NLS_Vector& NLS_PetraGroup::computeNewton() 
{
  cout << "ERROR: No direct methods are avialable for linear solves yet!\n"
       << "Use the iterative solver call!" << endl;
  throw;
}

//! Compute and return Newton direction, using desired accuracy for nonlinear solve
/*! Throws an error if RHS and Jacobian have not been computed */
const NLS_Vector& NLS_PetraGroup::computeNewton(NLS_ParameterList& p) 

{
  if (!isValidRHS) {
    cout << "ERROR: computeNewton() - RHS is out of date wrt X!" << endl;
    throw;
  }
  else if (!isValidJacobian) {
    cout << "ERROR: computeNewton() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  
  Epetra_LinearProblem Problem(Jac, 
			       &(NewtonVector.getPetraVector()), 
			       &(RHSVector.getPetraVector()));

  int maxit = p.getParameter("Max Linear Iterations", 400);
  double tol = p.getParameter("Linear Solver Tolerance", 1.0e-6);

  if (NLS_Utilities::doPrint(3)) {
    cout << "      *********************************************" << endl;
    cout << "      NLS_PetraGroup::computeNewton "<< endl;
    cout << "           Max Linear Iterations    = " << maxit << endl;
    cout << "           Linear Solver Tolerance  = " << tol << endl;
    cout << "      *********************************************" << endl 
	 << endl;
  }

  // For now, set problem level to hard, moderate, or easy
  Problem.SetPDL(easy);
  AztecOO aztec(Problem);
  aztec.Iterate(maxit,tol);
  isValidNewton = true;
  return NewtonVector;
}

double NLS_PetraGroup::computeLinearRHSNorm() 
{  
  if (isRHS() && isJacobian() && isNewton()) {
 
    bool TransJ = false;
    NLS_PetraVector* d = dynamic_cast <NLS_PetraVector*> (RHSVector.newcopy());
    Jac->Multiply(TransJ, NewtonVector.getPetraVector(), d->getPetraVector());
    d->update(RHSVector, *d, -1.0);
    double norm = d->norm();
    delete d;
    return norm;
  }
  cout << "ERROR: NLS_PetraGroup::computeLinearRHS - Vectors are out" 
    " of date wrt X!" << endl;
  throw;
}

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

const NLS_Vector& NLS_PetraGroup::getX() const 
{
  return xVector;
}

const NLS_Vector& NLS_PetraGroup::getRHS() const 
{  
  return RHSVector;
}

const NLS_Vector& NLS_PetraGroup::getGrad() const 
{ 
  return gradVector;
}

const NLS_Vector& NLS_PetraGroup::getNewton() const 
{
  return NewtonVector;
}

