// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_SCHURINTERFACE_H
#define NOX_EPETRA_SCHURINTERFACE_H

#include "NOX_Epetra.H"

//-----------------------------------------------------------------------------

// This is a new class that may eventually get moved into NOX::Epetra.  For now,
// this is simply used as a testbed for driving new coupling algorithms - RWH 10/23/2006

namespace NOX {

namespace Epetra {

class SchurInterface
{
public:

  SchurInterface() { };

  virtual ~SchurInterface() { };

  //! Apply the action of a matrix block on a vector
  virtual bool applyBlockAction( int rowBlock, int colBlock, const Epetra_Vector & x, Epetra_Vector & result ) = 0;

  //! Apply the action of a matrix block inverse on a vector
  virtual bool applyBlockInverseAction( int rowBlock, int colBlock, const Epetra_Vector & x, Epetra_Vector & result ) = 0;

  //! Get the communicator
  virtual const Epetra_Comm & getComm() = 0;

  //! Query availability of explicit operator
  virtual bool hasExplicitOperator( int rowBlock, int colBlock )
    { return false; }

  //! Access explicit operator
  virtual Teuchos::RCP<Epetra_Operator> getExplicitOperator( int rowBlock, int colBlock )
    { return Teuchos::null; }

};
} // Epetra
} // NOX

#endif
