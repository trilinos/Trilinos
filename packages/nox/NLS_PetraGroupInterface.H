
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NLS_PETRAGROUPINTERFACE_H
#define _NLS_PETRAGROUPINTERFACE_H
  
//! NLSPACK pure abstract interface to provide a link to the user's code 
//! for RHS and Jacobian fills called by NLS_PetraGroup.

#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"

class NLS_PetraGroupInterface {

public:
  
  virtual void computeRHS(Epetra_Vector& x, Epetra_Vector& RHS) = 0;
  
  virtual void computeJacobian(Epetra_Vector& x, Epetra_RowMatrix& Jac) = 0;

};

#endif
