// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_INTERFACE_REQUIRED_H
#define NOX_EPETRA_INTERFACE_REQUIRED_H

#include "NOX_Common.H"
  
// Forward declarations
class Epetra_Vector;

namespace NOX {
namespace Epetra {

  /*! 
    \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.  

    Contains interfaces for the user to supply (1) the evaluation of the nonlinear equations, (2) the Jacobian, and (3) any preconditioning if required.
  */
namespace Interface {

  /*!
    \brief Supplies NOX with the set nonlinear equations.
 
    This is the minimum required information to solve a nonlinear
    problem using the NOX::Epetra objects for the linear algebra
    implementation.  Used by NOX::Epetra::Group to provide a link
    to the external code for residual fills.
  */
class Required {

public:
 
  //! Type of fill that a computeF() method is used for.  
  /*! computeF() can be called for a variety of reasons:

  - To evaluate the function residuals.
  - To be used in an approximation to the Jacobian (finite difference or directional derivative).
  - To be used in an approximation to the preconditioner.
  
  This flag tells computeF() what the evaluation is used for.  This allows the user to change the fill process to eliminate costly terms.  For example, sometimes, terms in the function are very expensive and can be ignored in a Jacobian calculation.  The user can query this flag and determine not to recompute such terms if the computeF() is used in a Jacobian calculation.
   */
  enum FillType {
    //! The exact residual (F) is being calculated.
    Residual, 
    //! The Jacobian matrix is being estimated.
    Jac, 
    //! The preconditioner matrix is being estimated.
    Prec, 
    //! The fill context is from a FD approximation (includes FDC)
    FD_Res,
    //! The fill context is from a MF approximation
    MF_Res, 
    //! The fill context is from a MF computeJacobian() approximation
    MF_Jac, 
    //! A user defined estimation is being performed.
    User
  };

  //! Constructor
  Required() {};

  //! Destructor
  virtual ~Required() {};

  //! Compute the function, F, given the specified input vector x.  Returns true if computation was successful.
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			const FillType fillFlag) = 0;
  
};
} // namespace Interface
} // namespace Epetra
} // namespace NOX

#endif
