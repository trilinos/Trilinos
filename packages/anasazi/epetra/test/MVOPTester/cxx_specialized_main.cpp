// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
// @HEADER
//
//  This test uses the MVOPTester.hpp functions to test the Anasazi adapters
//  to Epetra and Thyra.
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

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziSpecializedEpetraAdapter.hpp"
#include "AnasaziBasicOutputManager.hpp"

int main(int argc, char *argv[])
{
  int i, epetra_ierr;
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

   // number of global elements
  int dim = 100;
  int blockSize = 5;

  bool verbose = false;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );
  
  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map->NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map->MyGlobalElements(&MyGlobalElements[0]);

  // Create an Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, 1) );
   
  // Add  rows one-at-a-time
  std::vector<double> Values(1, 1.0);
  for (i=0; i<NumMyElements; i++) {
    // Put in the diagonal entry
    epetra_ierr = A->InsertGlobalValues(MyGlobalElements[i],1,&Values[0],&MyGlobalElements[i]);
    assert(epetra_ierr==0);
  }
   
  // Finish building the epetra matrix A
  epetra_ierr = A->FillComplete();
  assert(epetra_ierr==0);

  // Issue several useful typedefs;
  typedef Anasazi::MultiVec<double> EMV;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note that this needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Anasazi::EpetraOpMultiVec> ivec = Teuchos::rcp( new Anasazi::EpetraOpMultiVec(A, *Map, blockSize) );
  ivec->MvRandom();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::BasicOutputManager<double>() );
  if (verbose) {
    MyOM->setVerbosity( Anasazi::Warnings );
  }

  // test the Epetra adapter multivector
  ierr = Anasazi::TestMultiVecTraits<double,EMV>(MyOM,ivec);
  gerr &= ierr;
  if (ierr) {
    MyOM->print(Anasazi::Warnings,"*** EpetraSpecializedAdapter PASSED TestMultiVecTraits()\n");
  }
  else {
    MyOM->print(Anasazi::Warnings,"*** EpetraSpecializedAdapter FAILED TestMultiVecTraits() ***\n\n");
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (gerr == false) {
    MyOM->print(Anasazi::Warnings,"End Result: TEST FAILED\n");
    return -1;
  }
  //
  // Default return value
  //
  MyOM->print(Anasazi::Warnings,"End Result: TEST PASSED\n");
  return 0;

}
