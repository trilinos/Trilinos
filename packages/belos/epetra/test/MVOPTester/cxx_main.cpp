//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//
//  This test uses the MVOPTester.hpp functions to test the Belos adapters
//  to Epetra.
//

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosOutputManager.hpp"

#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[])
{
  int i;
  bool ierr, gerr;
  gerr = true;

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  // If we aren't using MPI, then setup a serial communicator.
  Teuchos::RCP<Epetra_SerialComm> Comm = Teuchos::rcp( new Epetra_SerialComm() );
#endif

  bool verbose = false;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  bool success = true;
  try {
    // number of global elements
    int dim = 100;
    int blockSize = 5;

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );

    // Get update list and number of local equations from newly created Map.
    int NumMyElements = Map->NumMyElements();
    std::vector<int> MyGlobalElements(NumMyElements);
    Map->MyGlobalElements(&MyGlobalElements[0]);

    // Create an integer std::vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
    // on this processor
    std::vector<int> NumNz(NumMyElements);

    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    for (i=0; i<NumMyElements; i++) {
      if (MyGlobalElements[i]==0 || MyGlobalElements[i] == dim-1) {
        NumNz[i] = 2;
      }
      else {
        NumNz[i] = 3;
      }
    }

    // Create an Epetra_Matrix
    Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]) );

    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    std::vector<double> Values(2);
    Values[0] = -1.0; Values[1] = -1.0;
    std::vector<int> Indices(2);
    double two = 2.0;
    int NumEntries;
    for (i=0; i<NumMyElements; i++) {
      if (MyGlobalElements[i]==0) {
        Indices[0] = 1;
        NumEntries = 1;
      }
      else if (MyGlobalElements[i] == dim-1) {
        Indices[0] = dim-2;
        NumEntries = 1;
      }
      else {
        Indices[0] = MyGlobalElements[i]-1;
        Indices[1] = MyGlobalElements[i]+1;
        NumEntries = 2;
      }
      ierr = A->InsertGlobalValues(MyGlobalElements[i],NumEntries,&Values[0],&Indices[0]);
      assert(ierr==0);
      // Put in the diagonal entry
      ierr = A->InsertGlobalValues(MyGlobalElements[i],1,&two,&MyGlobalElements[i]);
      assert(ierr==0);
    }

    // Finish building the epetra matrix A
    ierr = A->FillComplete();
    assert(ierr==0);

    // Issue several useful typedefs;
    typedef Belos::MultiVec<double> EMV;
    //typedef Belos::Operator<double> EOP; // unused

    // Create an Epetra_MultiVector for an initial std::vector to start the solver.
    // Note that this needs to have the same number of columns as the blocksize.
    Teuchos::RCP<Belos::EpetraMultiVec> ivec = Teuchos::rcp( new Belos::EpetraMultiVec(*Map, blockSize) );
    ivec->Random();

    // Create an output manager to handle the I/O from the solver
    Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>() );
    if (verbose) {
      MyOM->setVerbosity( Belos::Warnings );
    }

    // test the Epetra adapter multivector
    ierr = Belos::TestMultiVecTraits<double,EMV>(MyOM,ivec);
    gerr &= ierr;
    if (ierr) {
      MyOM->print(Belos::Warnings,"*** EpetraAdapter PASSED TestMultiVecTraits()\n");
    }
    else {
      MyOM->print(Belos::Warnings,"*** EpetraAdapter FAILED TestMultiVecTraits() ***\n\n");
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if (!gerr) {
      success = false;
      MyOM->print(Belos::Warnings,"End Result: TEST FAILED\n");
    } else {
      success = true;
      MyOM->print(Belos::Warnings,"End Result: TEST PASSED\n");
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
