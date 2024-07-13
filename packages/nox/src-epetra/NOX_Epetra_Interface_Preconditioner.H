// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_INTERFACE_PRECONDITIONER_H
#define NOX_EPETRA_INTERFACE_PRECONDITIONER_H

#include "NOX_Common.H"

// Forward declarations
class Epetra_Vector;
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}

namespace NOX {
namespace Epetra {
namespace Interface {

  /*! \brief Used by NOX::Epetra to provide a link to the
    external code for Precondtioner fills.

    This is only required if the user wishes to supply their own
    preconditioner operator.
  */
class Preconditioner {

public:

  //! Constructor
  Preconditioner() {};

  //! Destructor
  virtual ~Preconditioner() {};

  //! Computes a user defined preconditioner.
  virtual bool computePreconditioner(const Epetra_Vector& x,
                     Epetra_Operator& M,
                     Teuchos::ParameterList* precParams = 0) = 0;

};
} // namespace Interface
} // namespace Epetra
} // namespace NOX

#endif
