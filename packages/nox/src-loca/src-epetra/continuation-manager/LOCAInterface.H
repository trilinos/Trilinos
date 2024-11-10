// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_INTERFACE_H
#define LOCA_INTERFACE_H

// Trilinos headers
#include <Epetra_Vector.h>
#include <LOCA_Epetra.H>

// ParaCont headers
#include "ProblemLOCAPrototype.H"
#include "ContinuationManager.H"

class  LOCAInterface :
//  public LOCA::Epetra::Interface::Required,
  public LOCA::Epetra::Interface::TimeDependent,
  public NOX::Epetra::Interface::Jacobian {

    public:

      // Constructor
      LOCAInterface( Teuchos::RCP <ProblemLOCAPrototype> & aProblem ,
      Teuchos::RCP <ContinuationManager> aContinuationManager);

      // Destructor
      ~LOCAInterface();

      // Compute and return F
      bool computeF(const Epetra_Vector& x, Epetra_Vector& f,
            const NOX::Epetra::Interface::Required::FillType F);

      // Compute an explicit Jacobian
      bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

      // Set a parameter in the user's code.
      void setParameters(const LOCA::ParameterVector& params);

      // Print the solution
      void printSolution (const Epetra_Vector &x, const double conParam);

      // Compute the Shifted matrix needed for linear stability
      virtual bool computeShiftedMatrix (double alpha, double beta,
                                 const Epetra_Vector &x, Epetra_Operator &A);

      virtual void setXdot(const Epetra_Vector& xdot, const double time);

    private:

      // The continuation manager
      Teuchos::RCP <ContinuationManager> continuationManager;

      // The Problem
      Teuchos::RCP <ProblemLOCAPrototype> problem;

      double t; // time, for non-autonomous problems

};

#endif
