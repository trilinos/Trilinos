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

#include "NOX_Common.H"
  
// Forward declarations
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_Operator;

namespace NOX {
namespace Epetra {

  /*! \brief Used by NOX::Epetra::Group to provide a link to the
   *  external code for RHS and Jacobian fills ()if needed. */
  /*! At a minimum the user must supply a computeF().  
   * If the user is supplying their own Epetra_RowMatrix for the Jacobian,
   * then they must supply a computeJacobian().  If computeJacobian is
   * not supplied, a Matrix-Free or Finite Differencing Epetra_RowMatrix
   * can be used in NOX.  Preconditioning functions are only  needed 
   * if the user wishes to supply their own preconditioner.  Otherwise
   * Aztec has built-in preconditioners to use.  */
class Interface {

public:
 
  //! Constructor
  Interface() {};

  //! Destructor
  virtual ~Interface() {};

  //! Compute RHS given the specified input vector x.  Returns true if computation was successful.
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& RHS) = 0;
  
  //! Compute Jacobian given the specified input vector x.  Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) = 0;

  //! Compute the Epetra_RowMatrix M, that will be used by the Aztec preconditioner instead of the Jacobian.  This is used when there is no explicit Jacobian present (i.e. Matrix-Free Newton-Krylov).  This MUST BE an Epetra_RowMatrix since the Aztec preconditioners need to know the sparsity pattern of the matrix.  Returns true if computation was successful.
  virtual bool computePreconditioner(const Epetra_Vector& x, Epetra_RowMatrix& M) = 0;

  //! Computes a user supplied preconditioner based on input vector x.  Returns true if computation was successful.
  virtual bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M) = 0;

};
} // namespace Epetra
} // namespace NOX

#endif
