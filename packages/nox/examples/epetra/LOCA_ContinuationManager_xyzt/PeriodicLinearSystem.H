// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

// Trilinos headers
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "LOCA_Parameter_Vector.H"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// ParaCont headers
#include "ProblemLOCAPrototype.H"

class PeriodicLinearSystem:
  public ProblemLOCAPrototype{

    public:

      // Constructor
      PeriodicLinearSystem( const Teuchos::RCP <Epetra_Comm> aComm );

      // Destructor
      ~PeriodicLinearSystem();

      // Implementation of ProblemNOXPrototype
      // virtual methods  ---------------------------------------
      //
      // Compute F
      virtual bool ComputeF(const Epetra_Vector & x,
          Epetra_Vector & f);

      // Compute the Jacobian of F
      virtual bool ComputeJacF(const Epetra_Vector & x);

      // Returns a pointer to the Jacobian of F
      virtual Teuchos::RCP <Epetra_CrsMatrix>
    GetJacF() const;

      // Returns a pointer to the initial guess
      virtual Teuchos::RCP <Epetra_Vector>
    GetInitialGuess() const;

      // Implementation of ProblemLOCAPrototype
      // virtual methods  --------------------------------------
      //
      // Pointer to the continuable parameter list
      virtual LOCA::ParameterVector
    GetContinuableParams() const;

      // Setting one continuation parameter
      virtual bool SetContinuableParameter(std::string label,
      double value);

      // Updates the continuation file
      virtual bool UpdateContinuationFile(
      const std::string & fileName,
          const int & idStep,
          const Teuchos::ParameterList & continuationFileParams);

      // Setting the continuation file parameters
      virtual bool SetContinuationFileParameters(const Epetra_Vector & x);

      // Getting the continuation file parameters
      virtual Teuchos::RCP <Teuchos::ParameterList>
    GetContinuationFileParameters();

      // Printing the solution at each step
      virtual bool PrintSolutionFile(const std::string & fileName,
      const Epetra_Vector & x,
      const Teuchos::ParameterList & xParams);

    private:

      // Internally stored variables ----------------------
      //
      // Communicator
      Teuchos::RCP <Epetra_Comm> comm;

      // Continuable parameter list
      LOCA::ParameterVector continuableParams;

      // The continuation file parameters
      Teuchos::RCP <Teuchos::ParameterList>
    continuationFileParams;

      // Initial Guess
      Teuchos::RCP <Epetra_Vector> initialGuess;

      // Jacobian
      Teuchos::RCP <Epetra_CrsMatrix>  jacobian;
      Teuchos::RCP <Epetra_CrsMatrix>  savedJacobian;

      // From Local to global indices (vector map)
      int * myGlobalElements;

      // Local number of elements (vector map)
      int numMyElements;

      // Epetra Map of the vector's components
      Teuchos::RCP <Epetra_Map>  vectorMap;
  };

#endif
