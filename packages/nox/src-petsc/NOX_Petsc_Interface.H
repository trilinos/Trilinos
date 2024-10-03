// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
