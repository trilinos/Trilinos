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
#include "Teuchos_StandardCatchMacros.hpp"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicOutputManager.hpp"

#ifdef HAVE_EPETRA_THYRA
#include "AnasaziThyraAdapter.hpp"
#include "AnasaziThyraDebugAdapter.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif

int runTest(
int argc,
char* argv[],
#ifdef HAVE_MPI
Teuchos::RCP<Epetra_MpiComm>& Comm
#else
Teuchos::RCP<Epetra_SerialComm>& Comm
#endif
){
  using std::cout;
  using std::endl;
  using Teuchos::rcp_implicit_cast;

  int i;
  bool ierr, gerr;
  gerr = true;

   // number of global elements
  int dim = 100;
  int blockSize = 5;

  bool verbose = false;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::BasicOutputManager<double>() );
  if (verbose) {
    MyOM->setVerbosity( Anasazi::Warnings );
  }

#ifndef HAVE_EPETRA_THYRA
  MyOM->stream(Anasazi::Warnings)
    << "Please configure Anasazi with:" << endl
    << "--enable-epetra-thyra" << endl
    << "--enable-anasazi-thyra" << endl;
  return -1;
#endif

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map->NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map->MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
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

  // Create an Anasazi::EpetraSymOp from this Epetra_CrsMatrix
  Teuchos::RCP<Anasazi::EpetraSymOp> op = Teuchos::rcp(new Anasazi::EpetraSymOp(A));

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note that this needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Anasazi::EpetraMultiVec> ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec(*Map, blockSize) );
  ivec->Random();

#ifdef HAVE_EPETRA_THYRA
  typedef Thyra::MultiVectorBase<double> TMVB;
  typedef Thyra::LinearOpBase<double>    TLOB;
  typedef Anasazi::MultiVec<double>      AMV;
  typedef Anasazi::Operator<double>      AOP;
  // create thyra objects from the epetra objects

  // first, a Thyra::VectorSpaceBase
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > epetra_vs =
    Thyra::create_VectorSpace(Map);

  // then, a MultiVectorBase (from the Epetra_MultiVector)
  Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_ivec =
    Thyra::create_MultiVector(rcp_implicit_cast<Epetra_MultiVector>(ivec),epetra_vs);

  // then, a LinearOpBase (from the Epetra_CrsMatrix)
  Teuchos::RCP<const Thyra::LinearOpBase<double> > thyra_op =
    Thyra::epetraLinearOp(A);

  // test the Thyra adapter multivector
  ierr = Anasazi::TestMultiVecTraits<double,TMVB>(MyOM,thyra_ivec);
  gerr |= ierr;
  if (ierr) {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter PASSED TestMultiVecTraits()" << endl;
  }
  else {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter FAILED TestMultiVecTraits() ***" << endl << endl;
  }

  // test the Thyra adapter operator
  ierr = Anasazi::TestOperatorTraits<double,TMVB,TLOB>(MyOM,thyra_ivec,thyra_op);
  gerr |= ierr;
  if (ierr) {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter PASSED TestOperatorTraits()" << endl;
  }
  else {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter FAILED TestOperatorTraits() ***" << endl << endl;
  }

  // test the Thyra debug adapter, ThyraMultiVec
  Teuchos::RCP<AMV> thyraMultiVec = Teuchos::rcp( new Anasazi::ThyraMultiVec<double>( thyra_ivec ) );
  ierr = Anasazi::TestMultiVecTraits<double,AMV>(MyOM,thyraMultiVec);
  gerr |= ierr;
  if (ierr) {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraDebugAdapter PASSED TestMultiVecTraits()" << endl;
  }
  else {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraDebugAdapter FAILED TestMultiVecTraits() ***" << endl << endl;
  }

  // test the Thyra debug adapter, ThyraOp
  Teuchos::RCP<AOP> thyraOp = Teuchos::rcp( new Anasazi::ThyraOp<double>( thyra_op ) );
  ierr = Anasazi::TestOperatorTraits<double,AMV,AOP>(MyOM,thyraMultiVec,thyraOp);
  gerr |= ierr;
  if (ierr) {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraDebugAdapter PASSED TestOperatorTraits()" << endl;
  }
  else {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraDebugAdapter FAILED TestOperatorTraits() ***" << endl << endl;
  }

  Teuchos::TimeMonitor::summarize(MyOM->stream(Anasazi::Warnings));
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


int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  // If we aren't using MPI, then setup a serial communicator.
  Teuchos::RCP<Epetra_SerialComm> Comm = Teuchos::rcp( new Epetra_SerialComm() );
#endif

  bool success = false;
  bool verbose = false;
  try {
    int toReturn = runTest(argc, argv, Comm);
    success = toReturn==0;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  Comm = Teuchos::null;
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
