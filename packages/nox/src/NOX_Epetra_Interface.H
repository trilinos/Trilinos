// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_EPETRA_INTERFACE_H
#define _NOX_EPETRA_INTERFACE_H
  
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"

namespace NOX {
namespace Epetra {

  /*! \brief Used by NOX Epetra Group to provide a link to the
   *  external code for RHS and Jacobian fills. */

class Interface {

public:
 
  //! Compute RHS given the specified input vector, x.
  virtual void computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS) = 0;
  
  //! Compute Jacobian given the specified input vector, x.
  virtual void computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac) = 0;

};
}
}

#endif
