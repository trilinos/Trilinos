// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Trilinos headers
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// ParaCont headers
#include "ContinuationManager.H"
#include "LinearSystem.H"

// Main driver
int main( int argc, char **argv )
{

  // Initialise MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  bool success = false;
  bool verbose = false;
  try {

#ifdef HAVE_MPI
    // Create a communicator
    Teuchos::RCP <Epetra_MpiComm> comm =
      Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    // Create a communicator
    Teuchos::RCP <Epetra_SerialComm> comm =
      Teuchos::rcp(new Epetra_SerialComm);
#endif

    std::string fileName = "task.xml";
    if (argc>1)
      fileName = argv[1];

    // Instantiate the continuation manager
    Teuchos::RCP <ContinuationManager> contManager =
      Teuchos::rcp(new ContinuationManager(comm,fileName));

    // Instantiate the problem
    Teuchos::RCP <LinearSystem> problem =
      Teuchos::rcp(new LinearSystem(comm));

    // Set the problem in the continuation manager
    contManager->SetLOCAProblem(problem);

    // Prepare to run LOCA
    contManager->BuildLOCAStepper();

    // Run LOCA
    success = contManager->RunLOCAStepper();

    if (success)
      std::cout << "\nAll tests passed" << std::endl;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
