//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
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
    printing.out() << "Test passed!" << std::endl;
  else 
    printing.out() << "Test failed!" << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  // Final return value (0 = successfull, non-zero = failure)
  return status;
}

/*
  end of file test.C
*/
