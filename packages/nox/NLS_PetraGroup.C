// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_PetraGroup.H"

//! Constructor
NLS_PetraGroup::NLS_PetraGroup() :
  xVector(NULL),
  RHSVector(NULL),
  gradVector(NULL),
  NewtonVector(NULL)
{
  resetVectorStatus();
}

//! Copy constructor
NLS_PetraGroup::NLS_PetraGroup(NLS_Group& copyFrom) {

}

//! Create a new group where the new solution vector is grp.x() + step * d
NLS_PetraGroup::NLS_PetraGroup(NLS_Group& grp, 
				     NLS_Vector& d, double step) {

}
 
//! NLS_Group deconstructor
NLS_PetraGroup::~NLS_PetraGroup() {

}

//! Copy Constructor
NLS_PetraGroup& NLS_PetraGroup::operator=(NLS_PetraGroup& 
							copyFrom) {

}

//! Compute and return solution vector
NLS_Vector& NLS_PetraGroup::computeX(NLS_Group& x, NLS_Vector& d, double step) {
  xVector.update(step, d, 1.0);
  resetVectorStatus();
}

//! Compute and return RHS
NLS_Vector& NLS_PetraGroup::computeRHS() {
  interface.loadFunction(xVector);
}

//! Compute RHS
void NLS_PetraGroup::computeJacobian() {
  interface.loadJacobian(xVector);
}

//! Compute and return gradient 
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_PetraGroup::computeGrad() {
  if (!isRHS) {
    cout << "ERROR: computeGrad() - RHS is out of date wrt X!" << endl;
    throw;
  }
  else if (!isJacobian) {
    cout << "ERROR: computeGrad() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  bool TransJ = true;
  double scaleFactor = 1.0/RHSVector.norm2();
  interface.Jacobian.Multiply(TransJ, RHSVector.petraVec, gradVector.petraVec);
  gradVector.scaleVector(scaleFactor);
  return gradVector;
}

//! Compute and return Newton direction 
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_PetraGroup::computeNewton() {
  cout << "ERROR: No direct methods are avialable for matrix inversion yet!\n"
    "Use the iterative solver!" << endl;
  throw;
}

//! Compute and return Newton direction, using desired accuracy for nonlinear solve
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_PetraGroup::computeNewton(string& name, NLS_Parameter& parameter) {
  if (!isRHS) {
    cout << "ERROR: computeNewton() - RHS is out of date wrt X!" << endl;
    throw;
  }
  else if (!isJacobian) {
    cout << "ERROR: computeNewton() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  Petra_RDP_LinearProblem Problem(interface.Jacobian,
				  xVector.petraVec,
				  RHSVector.petraVec);
  // For now, set problem level to hard, moderate, or easy
  Problem.SetPDL(hard);
  Aztec_OO aztec(Problem);
  aztec.Iterate(nlParamsPtr->getMaxLinearStep(), eta_k);
  return NewtonVector;
}

/** @name Checks to see if various objects have been computed. 
 *
 * Returns true if the corresponding "compute" function has been
 * called since the last update to the solution vector (via
 * instantiation or computeX). */

bool NLS_PetraGroup::isRHS() {   
  if (statusRHS) return true; 
  else return false; 
}
bool NLS_PetraGroup::isJacobian() {   
  if (statusJacobian) return true; 
  else return false; 
}
bool NLS_PetraGroup::isGrad() {   
  if (statusGrad) return true; 
  else return false; 
}
bool NLS_PetraGroup::isNewton() {   
  if (statusNewton) return true; 
  else return false;
}

//! Return solution vector
NLS_Vector& NLS_PetraGroup::getX() {return *xVector;}

//! Return rhs (throws an error if RHS has not been computed)
NLS_Vector& NLS_PetraGroup::getRHS() {  
  if (isRHS) return *RHSVector;
  else {
    cout << "ERROR: RHS Vector does NOT correspond to current Group "
      "solution vector!" << endl;
    throw;      
  }
}

//! Return gradient (throws an error if gradient has not been computed)
NLS_Vector& NLS_PetraGroup::getGrad() { 
  if (isGrad) return *gradVector;
  else {
    cout << "ERROR: Grad Vector does NOT correspond to current Group "
      "solution vector!" << endl;
    throw;      
  }
}

//! Return Newton direction (throws an error if newton direction has not been computed)
NLS_Vector& NLS_PetraGroup::getNewton() {
  if (isNewton) return *NewtonVector;
  else {
    cout << "ERROR: Newton Vector does NOT correspond to current Group "
      "solution vector!" << endl;
    throw;      
  }
}

void NLS_PetraGroup::resetVectorStatus() {
  statusRHS = false;
  statusJacobian = false;
  statusGrad = false;
  statusNewton = false;
}
