// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PROBLEM_NOX_PROTOTYPE_H
#define PROBLEM_NOX_PROTOTYPE_H

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

class ProblemNOXPrototype {

  public:

    // Constructor
    ProblemNOXPrototype();

    // Destructor
    virtual ~ProblemNOXPrototype();

    // Compute F
    virtual bool ComputeF(const Epetra_Vector & x,
    Epetra_Vector & f)=0;

    // Compute the Jacobian of F
    virtual bool ComputeJacF(const Epetra_Vector & x)=0;

    // Returns a pointer to the Jacobian of F
    virtual Teuchos::RCP <Epetra_CrsMatrix> GetJacF() const =0;

    // Returns a pointer to the initial guess
    virtual Teuchos::RCP <Epetra_Vector> GetInitialGuess() const =0;

};
#endif
