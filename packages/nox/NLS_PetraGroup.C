// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_PetraGroup.H"

NLS_PetraGroup::NLS_PetraGroup(Epetra_Vector& x) :
  xVector(NULL),
  RHSVector(NULL),
  gradVector(NULL),
  NewtonVector(NULL),
  Jac(NULL)
{
  xVector = new NLS_PetraVector(x);
  RHSVector = new NLS_PetraVector(x,false);
  gradVector = new NLS_PetraVector(x,false);
  NewtonVector = new NLS_PetraVector(x,false);
  resetVectorStatus();
  doDeleteJacobian = false;
}

NLS_PetraGroup::NLS_PetraGroup(Epetra_Vector& x, Epetra_RowMatrix& J) :
  xVector(NULL),
  RHSVector(NULL),
  gradVector(NULL),
  NewtonVector(NULL),
  Jac(&J)
{
  xVector = new NLS_PetraVector(x);
  RHSVector = new NLS_PetraVector(x,false);
  gradVector = new NLS_PetraVector(x,false);
  NewtonVector = new NLS_PetraVector(x,false);
  resetVectorStatus();
  doDeleteJacobian = false;
}

//! Copy constructor
NLS_PetraGroup::NLS_PetraGroup(NLS_Group& copyFrom, bool newJacobian) 
{
  // Copy vectors
  xVector = dynamic_cast<NLS_PetraVector*>((copyFrom.getX()).newcopy());
  RHSVector = dynamic_cast<NLS_PetraVector*>((copyFrom.getRHS()).newcopy());
  gradVector = dynamic_cast<NLS_PetraVector*>((copyFrom.getGrad()).newcopy());
  NewtonVector = dynamic_cast<NLS_PetraVector*>((copyFrom.getNewton()).newcopy());
  // Copy Status Flags
  statusRHS = copyFrom.isRHS();
  statusJacobian = copyFrom.isJacobian();
  statusGrad = copyFrom.isGrad();
  statusNewton = copyFrom.isNewton();
  // Copy Jacobian Matrix if required
  NLS_PetraGroup& tmpGroup = dynamic_cast<NLS_PetraGroup&> (copyFrom);
  if (newJacobian) {
    Jac = new Epetra_CrsMatrix(*(dynamic_cast<Epetra_CrsMatrix*> (tmpGroup.Jac)));
    doDeleteJacobian = true;
  } else {
    Jac = tmpGroup.Jac;
    doDeleteJacobian = false;
  }
}

//! Create a new group where the new solution vector is grp.x() + step * d
NLS_PetraGroup::NLS_PetraGroup(NLS_Group& grp, NLS_Vector& d, double step,
			       bool newJacobian) 
{
  // Copy vectors
  xVector = dynamic_cast<NLS_PetraVector*>((grp.getX()).newcopy());
  xVector->update(1.0,d,step);
  RHSVector = dynamic_cast<NLS_PetraVector*>((grp.getRHS()).newcopy());
  gradVector = dynamic_cast<NLS_PetraVector*>((grp.getGrad()).newcopy());
  NewtonVector = dynamic_cast<NLS_PetraVector*>((grp.getNewton()).newcopy());
  // Copy Status Flags
  statusRHS = grp.isRHS();
  statusJacobian = grp.isJacobian();
  statusGrad = grp.isGrad();
  statusNewton = grp.isNewton();
  // Copy Jacobian Matrix if required
  NLS_PetraGroup& tmpGroup = dynamic_cast<NLS_PetraGroup&> (grp);
  if (newJacobian) {
    Jac = new Epetra_CrsMatrix(*(dynamic_cast<Epetra_CrsMatrix*> (tmpGroup.Jac)));
    doDeleteJacobian = true;
  } else {
    Jac = tmpGroup.Jac;
    doDeleteJacobian = false;
  }
}
 
//! NLS_Group deconstructor
NLS_PetraGroup::~NLS_PetraGroup() 
{
  delete xVector;
  delete RHSVector;
  delete gradVector;
  delete NewtonVector;
  if (doDeleteJacobian) delete Jac;

}

//! Copies the values of all vectors in source group to this group
NLS_Group& NLS_PetraGroup::copy(const NLS_Group& source)
{
    cout << "ERROR: NLS_PetraGroup::copy() - requires a PetraGroup object be passed in!" << endl;
  throw;
}

//! Copies the values of all vectors in source group to this group
NLS_Group& NLS_PetraGroup::copy(const NLS_PetraGroup& source)
{
  // Update the vectors
  xVector->copy(*source.xVector, 1.0);
  RHSVector->copy(*source.RHSVector, 1.0);
  gradVector->copy(*source.gradVector, 1.0);
  NewtonVector->copy(*source.NewtonVector, 1.0);
  // Update the statusVectors
  statusRHS = source.statusRHS;
  statusJacobian = source.statusJacobian;
  statusGrad = source.statusGrad;
  statusNewton = source.statusNewton;
}

//! Compute and return solution vector
NLS_Vector& NLS_PetraGroup::computeX(const NLS_Group& x, const NLS_Vector& d, double step) {
  cout << "ERROR: NLS_PetraGroup::computeX() - Pass Petra objects in call!" << endl;
  throw;
}

//! Compute and return solution vector
NLS_Vector& NLS_PetraGroup::computeX(const NLS_PetraGroup& x, 
				     const NLS_PetraVector& d, 
				     double step) {
  xVector->copy(x.getX());
  xVector->update(1.0, d, step);
  resetVectorStatus();
  return *xVector;
}

//! Compute and return RHS
NLS_Vector& NLS_PetraGroup::computeRHS() {
  cout << "ERROR: NLS_PetraGroup::computeRHS() - This function must be overloaded to\n"
    "interface with user's code!" << endl;
  exit(-1);
}

//! Compute RHS
void NLS_PetraGroup::computeJacobian() {
  cout << "ERROR: NLS_PetraGroup::computeJacobian() - This function must be overloaded to\n"
    "interface with user's code!" << endl;
  exit(-1);
}

//! Compute and return gradient 
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_PetraGroup::computeGrad() {
  if (!isRHS()) {
    cout << "ERROR: NLS_PetraGroup::computeGrad() - RHS is out of date wrt X!" << endl;
    throw;
  }
  if (!isJacobian()) {
    cout << "ERROR: NLS_PetraGroup::computeGrad() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  bool TransJ = true;
  double scaleFactor = 1.0/RHSVector->norm();
  //interface.Jacobian.Multiply(TransJ, RHSVector.petraVec, gradVector.petraVec);
  //gradVector->scale(scaleFactor);
  return *gradVector;
}

//! Compute and return Newton direction 
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_PetraGroup::computeNewton() {
  cout << "ERROR: No direct methods are avialable for matrix inversion yet!\n"
    "Use the iterative solver call!" << endl;
  throw;
}

//! Compute and return Newton direction, using desired accuracy for nonlinear solve
/*! Throws an error if RHS and Jacobian have not been computed */
NLS_Vector& NLS_PetraGroup::computeNewton(NLS_ParameterList& parameter) {
  if (!isRHS()) {
    cout << "ERROR: computeNewton() - RHS is out of date wrt X!" << endl;
    throw;
  }
  else if (!isJacobian()) {
    cout << "ERROR: computeNewton() - Jacobian is out of date wrt X!" << endl;
    throw;
  }
  Epetra_LinearProblem Problem(Jac, xVector->petraVec, RHSVector->petraVec);
  // For now, set problem level to hard, moderate, or easy
  Problem.SetPDL(hard);
  //AztecOO aztec(Problem);
  //aztec.Iterate(nlParamsPtr->getMaxLinearStep(), eta_k);
  return *NewtonVector;
}

/** @name Checks to see if various objects have been computed. 
 *
 * Returns true if the corresponding "compute" function has been
 * called since the last update to the solution vector (via
 * instantiation or computeX). */

bool NLS_PetraGroup::isRHS() const {   
  if (statusRHS) return true; 
  else return false; 
}
bool NLS_PetraGroup::isJacobian() const {   
  if (statusJacobian) return true; 
  else return false; 
}
bool NLS_PetraGroup::isGrad() const {   
  if (statusGrad) return true; 
  else return false; 
}
bool NLS_PetraGroup::isNewton() const {   
  if (statusNewton) return true; 
  else return false;
}

//! Return solution vector
NLS_Vector& NLS_PetraGroup::getX() const {return *xVector;}

//! Return rhs (throws an error if RHS has not been computed)
NLS_Vector& NLS_PetraGroup::getRHS() const {  
  if (isRHS()) return *RHSVector;
  else {
    cout << "ERROR: RHS Vector does NOT correspond to current Group "
      "solution vector!" << endl;
    throw;      
  }
}

//! Return gradient (throws an error if gradient has not been computed)
NLS_Vector& NLS_PetraGroup::getGrad() const { 
  if (isGrad()) return *gradVector;
  else {
    cout << "ERROR: Grad Vector does NOT correspond to current Group "
      "solution vector!" << endl;
    throw;      
  }
}

//! Return Newton direction (throws an error if newton direction has not been computed)
NLS_Vector& NLS_PetraGroup::getNewton() const {
  if (isNewton()) return *NewtonVector;
  else {
    cout << "ERROR: Newton Vector does NOT correspond to current Group "
      "solution vector!" << endl;
    throw;      
  }
}

NLS_Group* NLS_PetraGroup::newcopy(bool newJacobian) 
{
  // **************** Not Finished yet!! ****************
}

void NLS_PetraGroup::resetVectorStatus() 
{
  statusRHS = false;
  statusJacobian = false;
  statusGrad = false;
  statusNewton = false;
}
