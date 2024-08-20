// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// NOX headers
#include "NOX.H"  // Required headers
#include "NOX_Epetra.H" // Epetra Interface headers
#include "NOX_TestError.H" // Test Suite headers

// Trilinos headers
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[]) {
  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  bool success = true;
  try {
    // Get the process ID and the total number of processors
    int MyPID = Comm.MyPID();
#ifdef HAVE_MPI
    int NumProc = Comm.NumProc();
#endif

    // Set up the printing utilities
    Teuchos::RCP<Teuchos::ParameterList> noxParamsPtr =
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& noxParams = *(noxParamsPtr.get());
    // Only print output if the "-v" flag is set on the command line
    Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
    printParams.set("MyPID", MyPID);
    printParams.set("Output Precision", 5);
    printParams.set("Output Processor", 0);
    if( verbose )
      printParams.set("Output Information",
          NOX::Utils::OuterIteration +
          NOX::Utils::OuterIterationStatusTest +
          NOX::Utils::InnerIteration +
          NOX::Utils::Parameters +
          NOX::Utils::Details +
          NOX::Utils::Warning +
          NOX::Utils::TestDetails);
    else
      printParams.set("Output Information", NOX::Utils::Error +
          NOX::Utils::TestDetails);

    NOX::Utils printing(printParams);

    // Identify the test problem
    if (printing.isPrintType(NOX::Utils::TestDetails))
      printing.out() << "Starting epetra/NOX_NewTest/NOX_NewTest.exe" << std::endl;

    // Identify processor information
#ifdef HAVE_MPI
    if (printing.isPrintType(NOX::Utils::TestDetails)) {
      printing.out() << "Parallel Run" << std::endl;
      printing.out() << "Number of processors = " << NumProc << std::endl;
      printing.out() << "Print Process = " << MyPID << std::endl;
    }
    Comm.Barrier();
    if (printing.isPrintType(NOX::Utils::TestDetails))
      printing.out() << "Process " << MyPID << " is alive!" << std::endl;
    Comm.Barrier();
#else
    if (printing.isPrintType(NOX::Utils::TestDetails))
      printing.out() << "Serial Run" << std::endl;
#endif

    // *** Insert your testing here! ***

    // Final return value (0 = successfull, non-zero = failure)
    int status = 0;

    // Summarize test results
    if (status == 0)
    {
      success = true;
      printing.out() << "Test passed!" << std::endl;
    }
    else
    {
      success = false;
      printing.out() << "Test failed!" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

/*
  end of file test.C
*/
