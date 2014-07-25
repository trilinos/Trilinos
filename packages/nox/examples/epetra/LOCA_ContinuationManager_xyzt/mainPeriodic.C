/*
//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
*/

// Trilinos headers
#ifdef HAVE_MPI
#include "mpi.h"
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
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
    Teuchos::RefCountPtr <EpetraExt::MultiMpiComm> globalComm =
      Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD,
            numSpatialProcs,
            numTimeSteps));

    // Get the communicator for the spatial sub-problem
    //Teuchos::RefCountPtr <Epetra_MpiComm> comm =
    //   Teuchos::rcp(&(globalComm->SubDomainComm()), false);
#else
    // Create a two-level communicator for space&time parallelism
    int numTimeSteps = 2;
    Teuchos::RefCountPtr <EpetraExt::MultiSerialComm> globalComm =
      Teuchos::rcp(new EpetraExt::MultiSerialComm(numTimeSteps));

    // Get the communicator for the spatial sub-problem
    // Teuchos::RefCountPtr <Epetra_SerialComm> comm =
    //   Teuchos::rcp(&(globalComm->SubDomainComm()), false);
#endif
    Teuchos::RefCountPtr <Epetra_Comm> comm =
      Teuchos::rcp(&(globalComm->SubDomainComm()), false);

    std::string fileName = "task.xml";
    if (argc>1)
      fileName = argv[1];

    // Instantiate the continuation manager
    Teuchos::RefCountPtr <ContinuationManager> contManager =
      Teuchos::rcp(new ContinuationManager(comm,fileName));

    // Instantiate the problem
    Teuchos::RefCountPtr <PeriodicLinearSystem> problem =
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
