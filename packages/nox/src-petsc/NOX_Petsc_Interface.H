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

#ifndef NOX_PETSC_INTERFACE_H
#define NOX_PETSC_INTERFACE_H
  
// Forward declarations

namespace NOX {
namespace Petsc {

  /*! \brief Used by NOX::Petsc::Group to provide a link to the
   *  external code for F and Jacobian evaluations (if needed). */
  /*! At a minimum the user must supply a computeF().  
   * If the user is supplying their own Petsc Mat for the Jacobian,
   * then they must supply a computeJacobian().  Preconditioning 
   * functions are only needed if the user wishes to supply their own 
   * preconditioner.  Otherwise Petsc SLES has built-in preconditioners 
   * to use.  */
class Interface {

public:
 
  //! Constructor
  Interface() {};

  //! Destructor
  virtual ~Interface() {};

  //! Compute F given the specified input vector, x.  Returns true if computation was successful.
  virtual bool computeF(const Vec& x, Vec& RHS) = 0;
  
  //! Compute Jacobian given the specified input vector, x.  Returns true if computation was successful.
  virtual bool computeJacobian(const Vec& x, Mat& Jac) = 0;

  //! Compute the matrix M that will be used as the preconditioner.  Returns true if computation was successful.
  virtual bool computePreconditioner(Mat& M) = 0;

  //! Return the action of the preconditioner on a vector.  Returns true if computation was successful.
  virtual bool preconditionVector(Vec& y) = 0;

};
} // namespace Petsc
} // namespace NOX

#endif
