// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    A->InsertGlobalValues(MyGlobalElements[i],1,&Values[0],&MyGlobalElements[i]);
  }
   
  // Finish building the epetra matrix A
  A->FillComplete();

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
