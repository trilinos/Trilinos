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

#ifndef NOX_EPETRA_INTERFACE_H
#define NOX_EPETRA_INTERFACE_H
  
// External Includes from Trilinos
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"

namespace NOX {
namespace Epetra {

  /*! \brief Used by NOX::Epetra::Group to provide a link to the
   *  external code for RHS and Jacobian fills. */
  /*! At a minimum the user should only supply a RHS fill.  Jacobian fills 
   * and preconditioning options need only be supplied if certain methods 
   * require it. Therefore the implemetations below are set to throw an
   * error if not supplied.*/

class Interface {

public:
 
  //! Constructor
  Interface() {};

  //! Destructor
  virtual ~Interface() {};

  //! Compute RHS given the specified input vector, x.
  virtual void computeRHS(const Epetra_Vector& x, Epetra_Vector& RHS) = 0;
  
  //! Compute Jacobian given the specified input vector, x.
  virtual void computeJacobian(const Epetra_Vector& x, Epetra_RowMatrix& Jac) = 0;

  //! Compute the matrix M that will be used as the preconditioner.
  virtual void computePreconditioner(Epetra_RowMatrix& M) = 0;

  //! Return the action of the preconditioner on a vector.
  virtual void preconditionVector(Epetra_Vector& y) = 0;

};
} // namespace Epetra
} // namespace NOX

#endif
