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
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// ParaCont headers
#include "ContinuationManager.H"
#include "PeriodicLinearSystem.H"

// Main driver
int main( int argc, char **argv )
{
  // Initialise MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  bool success = false;
  bool verbose = false;
  try {
#ifdef HAVE_MPI
    // Create a two-level communicator for space&time parallelism
    int numSpatialProcs = 1;
    int numTimeSteps = 2;
    Teuchos::RCP <EpetraExt::MultiMpiComm> globalComm =
      Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD,
            numSpatialProcs,
            numTimeSteps));

    // Get the communicator for the spatial sub-problem
    //Teuchos::RCP <Epetra_MpiComm> comm =
    //   Teuchos::rcp(&(globalComm->SubDomainComm()), false);
#else
    // Create a two-level communicator for space&time parallelism
    int numTimeSteps = 2;
    Teuchos::RCP <EpetraExt::MultiSerialComm> globalComm =
      Teuchos::rcp(new EpetraExt::MultiSerialComm(numTimeSteps));

    // Get the communicator for the spatial sub-problem
    // Teuchos::RCP <Epetra_SerialComm> comm =
    //   Teuchos::rcp(&(globalComm->SubDomainComm()), false);
#endif
    Teuchos::RCP <Epetra_Comm> comm =
      Teuchos::rcp(&(globalComm->SubDomainComm()), false);

    std::string fileName = "task.xml";
    if (argc>1)
      fileName = argv[1];

    // Instantiate the continuation manager
    Teuchos::RCP <ContinuationManager> contManager =
      Teuchos::rcp(new ContinuationManager(comm,fileName));

    // Instantiate the problem
    Teuchos::RCP <PeriodicLinearSystem> problem =
      Teuchos::rcp(new PeriodicLinearSystem(comm));

    // Set the problem in the continuation manager
    contManager->SetLOCAProblem(problem);

    // Prepare to run LOCA
    contManager->BuildLOCAPeriodicStepper(globalComm);

    // Run LOCA
    success = contManager->RunLOCAStepper();

    if (success)
      std::cout << "\nAll tests passed" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
