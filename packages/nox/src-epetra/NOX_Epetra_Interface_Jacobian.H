// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_INTERFACE_JACOBIAN_H
#define NOX_EPETRA_INTERFACE_JACOBIAN_H

#include "NOX_Common.H"

// Forward declarations
class Epetra_Vector;
class Epetra_Operator;

namespace NOX {
namespace Epetra {
namespace Interface {

  /*! \brief Used by NOX::Epetra to provide a link to the
    external code for Jacobian fills.

    This is only required if the user wishes to supply their own Jacobian
    operator.
  */
class Jacobian {

public:

  //! Constructor.
  Jacobian() {};

  //! Destructor.
  virtual ~Jacobian() {};

  /*! Compute Jacobian given the specified input vector x.  Returns
    true if computation was successful.
   */
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) = 0;

};
} // namespace Interface
} // namespace Epetra
} // namespace NOX

#endif
