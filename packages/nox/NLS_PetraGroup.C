// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_PetraGroup.H"

#include "NLS_Utilities.H"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h" 

NLS_PetraGroup::NLS_PetraGroup(Epetra_Vector& x, NLS_PetraSharedJacobian& sj, 
			       NLS_PetraGroupInterface& i) :
  xVector(x), // deep copy x     
  RHSVector(x, NLS_PetraVector::CopyShape), // new vector of same size
  gradVector(x, NLS_PetraVector::CopyShape), // new vector of same size
  NewtonVector(x, NLS_PetraVector::CopyShape),// new vector of same size
  sharedJacobian(sj), // pass J to SharedJacobian
  interface(i) // set reference
{
  resetIsValid();
}

NLS_PetraGroup::NLS_PetraGroup(const NLS_PetraGroup& copyFrom) :
  xVector(copyFrom.xVector.getPetraVector()), 
  RHSVector(copyFrom.RHSVector.getPetraVector()), 
  gradVector(copyFrom.gradVector.getPetraVector()), 
  NewtonVector(copyFrom.NewtonVector.getPetraVector()),
  sharedJacobian(copyFrom.sharedJacobian),
  interface(copyFrom.interface)
{
  isValidRHS = copyFrom.isValidRHS;
  isValidGrad = copyFrom.isValidGrad;
  isValidNewton = copyFrom.isValidNewton;
  isValidJacobian = copyFrom.isValidJacobian;

  // New copy takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedJacobian.takeOwnership(this);
}

NLS_PetraGroup::~NLS_PetraGroup() 
{
}

void NLS_PetraGroup::resetIsValid() //private
{
  isValidRHS = false;
  isValidJacobian = false;
  isValidGrad = false;
  isValidNewton = false;
}

NLS_Group* NLS_PetraGroup::newCopy() const 
{
  NLS_PetraGroup* newgrp = new NLS_PetraGroup(*this);
  return newgrp;
}

NLS_Group& NLS_PetraGroup::operator=(const NLS_Group& copyFrom)
{
  return operator=(dynamic_cast<const NLS_PetraGroup&> (copyFrom));
}

NLS_Group& NLS_PetraGroup::operator=(const NLS_PetraGroup& copyFrom)
{
  // Copy the xVector
  xVector = copyFrom.xVector;

  // Copy reference to sharedJacobian
  sharedJacobian = copyFrom.sharedJacobian;

  // Update the isValidVectors
  isValidRHS = copyFrom.isValidRHS;
  isValidGrad = copyFrom.isValidGrad;
  isValidNewton = copyFrom.isValidNewton;
  isValidJacobian = copyFrom.isValidJacobian;

  // Only copy vectors that are valid
  if (isValidRHS)
    RHSVector = copyFrom.RHSVector;
  if (isValidGrad)
    gradVector = copyFrom.gradVector;
  if (isValidNewton)
    NewtonVector = copyFrom.NewtonVector;

  // If valid, this takes ownership of the shared Jacobian
  if (isValidJacobian)
    sharedJacobian.takeOwnership(this);
    
  return *this;
}

const NLS_Vector& NLS_PetraGroup::computeX(const NLS_Group& grp, const NLS_Vector& d, double step) 
{
  // Cast to appropriate type, then call other computeX
  const NLS_PetraGroup& petragrp = dynamic_cast<const NLS_PetraGroup&> (grp);
  const NLS_PetraVector& petrad = dynamic_cast<const NLS_PetraVector&> (d);
  return computeX(petragrp, petrad, step); 
}

//! Compute and return solution vector
const NLS_Vector& NLS_PetraGroup::computeX(const NLS_PetraGroup& grp, const NLS_PetraVector& d, 
					   double step) 
{
  resetIsValid();
  xVector.update(grp.xVector, d, step);
  return xVector;
}

const NLS_Vector& NLS_PetraGroup::computeRHS() 
{
  if (isRHS())
    return RHSVector;

  interface.computeRHS(xVector.getPetraVector(), RHSVector.getPetraVector());
  normRHS = RHSVector.norm();
  isValidRHS = true;
  return RHSVector;
}

void NLS_PetraGroup::computeJacobian() 
{
  // Skip if the Jacobian is already valid
  if (isJacobian())
    return;

  // Take ownership of the Jacobian and get a reference
  Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian(this);

  // Fill the Jacobian
  interface.computeJacobian(xVector.getPetraVector(), Jacobian);

  // Update status
  isValidJacobian = true;
}

const NLS_Vector& NLS_PetraGroup::computeGrad() 
{
  if (isGrad())
    return gradVector;
  
  if (!isRHS()) {
    cout << "ERROR: NLS_PetraGroup::computeGrad() - RHS is out of date wrt X!" << endl;
    throw;
  }

  if (!isJacobian()) {
    cout << "ERROR: NLS_PetraGroup::computeGrad() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  
  // Get a reference to the Jacobian (it's validity was check above)
  const Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian();

  // Compute grad = Jacobian^T * RHS.
  bool Transpose = true;
  Jacobian.Multiply(Transpose, RHSVector.getPetraVector(), gradVector.getPetraVector());

  // Update state
  isValidGrad = true;

  // Return result
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
  if (isNewton())
    return NewtonVector;

  if (!isRHS()) {
    cout << "ERROR: computeNewton() - RHS is out of date wrt X!" << endl;
    throw;
  }

  // Check the we own the shared Jacobian AND that the Jacobian is valid
  if (!isJacobian()) {
    cout << "ERROR: computeNewton() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  
  // Get the Jacobian 
  /* (Have to get non-const version which requires reasserting
     ownership. This is due to a flaw in Epetra_LinearProblem). */
  Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian(this);

  // Create Epetra problem to solve for NewtonVector 
  // (Jacobian * NewtonVector = RHSVector)
  Epetra_LinearProblem Problem(&Jacobian, 
			       &(NewtonVector.getPetraVector()), 
			       &(RHSVector.getPetraVector()));

  // For now, set problem level to hard, moderate, or easy
  Problem.SetPDL(easy);

  // Create aztec problem
  AztecOO aztec(Problem);

  int maxit = p.getParameter("Max Linear Iterations", 400);
  double tol = p.getParameter("Linear Solver Tolerance", 1.0e-6);

  if (NLS_Utilities::doPrint(3)) {
    cout << "      *********************************************" << "\n";
    cout << "      NLS_PetraGroup::computeNewton "<< "\n";
    cout << "           Max Linear Iterations    = " << maxit << "\n";
    cout << "           Linear Solver Tolerance  = " << tol << "\n";
    cout << "      *********************************************" << "\n" << endl;
  }

  // Solve Aztex problem
  aztec.Iterate(maxit, tol);

  // Scale soln by -1
  NewtonVector.scale(-1.0);

  // Update state
  isValidNewton = true;

  // Return solution
  return NewtonVector;
}

double NLS_PetraGroup::computeNormPredictedRHS(const NLS_Vector& s) 
{  
  if (!isRHS()) {
    cout << "ERROR: computePredictedRHSNorm() - RHS is out of date wrt X!" << endl;
    throw;
  }
  if (!isJacobian()) {
    cout << "ERROR: computePredictedRHSNorm() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  
  // Cast s to a const NLS_PetraVector
  const NLS_PetraVector& petras = dynamic_cast <const NLS_PetraVector&> (s);

  // Create a new vector to be the predicted RHS
  NLS_PetraVector v(RHSVector.getPetraVector(), NLS_PetraVector::CopyShape);

  // Get a reference to the Jacobian (it's validity was check above)
  const Epetra_RowMatrix& Jacobian = sharedJacobian.getJacobian();

  // Compute v = Jacobian * petras
  bool NoTranspose = false;
  Jacobian.Multiply(NoTranspose, petras.getPetraVector(), v.getPetraVector());

  // Compute v = RHSVector + v (this is the predicted RHS)
  v.update(RHSVector, v, 1.0);

  // Return norm of predicted RHS
  return v.norm();
}

bool NLS_PetraGroup::isRHS() const 
{   
  return isValidRHS;
}

bool NLS_PetraGroup::isJacobian() const 
{  
  return ((sharedJacobian.isOwner(this)) && (isValidJacobian));
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

double NLS_PetraGroup::getNormRHS() const
{
  return normRHS;
}

const NLS_Vector& NLS_PetraGroup::getGrad() const 
{ 
  return gradVector;
}

const NLS_Vector& NLS_PetraGroup::getNewton() const 
{
  return NewtonVector;
}

