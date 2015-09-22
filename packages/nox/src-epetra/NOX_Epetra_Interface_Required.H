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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
